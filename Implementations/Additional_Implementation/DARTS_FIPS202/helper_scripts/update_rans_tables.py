#!/usr/bin/env python3
"""Validate rANS statistics and synchronize all checked-in table copies.

The script is check-only unless --write is supplied.  It consumes one
test_rANS report for each security mode, updates the frequency inputs in all
prec_encoding.c copies in memory, runs the real C precomputation program, and
then compares or updates the six generated tables in every encoding.c copy.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import shlex
import shutil
import stat
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


SCALE = 1024
UINT64_MAX = (1 << 64) - 1
MODES = (128, 256, 512)
EXPECTED_IMPLEMENTATIONS = 15
TABLE_NAMES = (
    "esyms_h",
    "dsyms_h",
    "symbol_h",
    "esyms_hb_z1",
    "dsyms_hb_z1",
    "symbol_hb_z1",
)


@dataclass(frozen=True)
class ModeSpec:
    params: Tuple[int, int, int, int, int]
    bounds: Tuple[str, str, str]
    m_h: int
    m_hb_z1: int


MODE_SPECS = {
    128: ModeSpec((512, 130817, 10, 2, 1), ("5866.64", "6054.88", "6057.90"), 23, 11),
    256: ModeSpec(
        (512, 130817, 9, 3, 2), ("18440.32", "18910.70", "18916.70"), 27, 13
    ),
    512: ModeSpec(
        (1024, 260609, 10, 3, 2), ("73440.24", "74216.86", "74220.97"), 71, 35
    ),
}


@dataclass(frozen=True)
class Report:
    mode: int
    path: Path
    samples: int
    mapped_candidates: int
    f_h: Tuple[int, ...]
    f_hb_z1: Tuple[int, ...]
    sha256: str
    contents: bytes


@dataclass(frozen=True)
class Definition:
    start: int
    end: int
    statement: str
    ctype: str
    initializer: str


@dataclass(frozen=True)
class GeneratedTable:
    statement: str
    ctype: str
    values: Tuple[int, ...]


class UpdateError(RuntimeError):
    pass


def read_text(path: Path) -> str:
    try:
        return path.read_bytes().decode("utf-8")
    except OSError as exc:
        raise UpdateError(f"cannot read {path}: {exc}") from exc
    except UnicodeDecodeError as exc:
        raise UpdateError(f"{path} is not valid UTF-8: {exc}") from exc


def one_match(pattern: str, text: str, label: str) -> re.Match[str]:
    matches = list(re.finditer(pattern, text, re.MULTILINE))
    if len(matches) != 1:
        raise UpdateError(f"expected exactly one {label}, found {len(matches)}")
    return matches[0]


def parse_frequency(text: str, name: str, expected_len: int, path: Path) -> Tuple[int, ...]:
    match = one_match(
        rf"^{re.escape(name)}\s*=\s*\{{([^}}\r\n]*)\}}\s*$",
        text,
        f"{name} line in {path}",
    )
    pieces = [piece.strip() for piece in match.group(1).split(",")]
    if not pieces or any(not re.fullmatch(r"[0-9]+", piece) for piece in pieces):
        raise UpdateError(f"{path}: {name} must be a comma-separated decimal integer list")
    values = tuple(int(piece) for piece in pieces)
    if len(values) != expected_len:
        raise UpdateError(
            f"{path}: {name} has {len(values)} entries, expected {expected_len}"
        )
    if any(value <= 0 for value in values):
        raise UpdateError(f"{path}: every {name} entry must be positive")
    if sum(values) != SCALE:
        raise UpdateError(f"{path}: {name} sums to {sum(values)}, expected {SCALE}")
    return values


def parse_decimal_field(text: str, name: str, path: Path) -> int:
    match = one_match(
        rf"^{re.escape(name)}\s*=\s*(\d+)\s*$",
        text,
        f"{name} field in {path}",
    )
    value = int(match.group(1))
    if value > UINT64_MAX:
        raise UpdateError(f"{path}: {name} exceeds uint64_t")
    return value


def parse_count_array(
    text: str, name: str, expected_len: int, path: Path
) -> Tuple[int, ...]:
    match = one_match(
        rf"^{re.escape(name)}\s*=\s*\{{([^}}\r\n]*)\}}\s*$",
        text,
        f"{name} line in {path}",
    )
    pieces = [piece.strip() for piece in match.group(1).split(",")]
    if not pieces or any(not re.fullmatch(r"[0-9]+", piece) for piece in pieces):
        raise UpdateError(
            f"{path}: {name} must be a comma-separated non-negative integer list"
        )
    values = tuple(int(piece) for piece in pieces)
    if len(values) != expected_len:
        raise UpdateError(
            f"{path}: {name} has {len(values)} entries, expected {expected_len}"
        )
    if any(value > UINT64_MAX for value in values):
        raise UpdateError(f"{path}: {name} contains a value exceeding uint64_t")
    return values


def normalize_counts(counts: Sequence[int], name: str, path: Path) -> Tuple[int, ...]:
    """Reproduce test_rANS.c normalize_frequencies exactly."""
    total = sum(counts)
    if total == 0:
        raise UpdateError(f"{path}: cannot normalize empty {name} counts")

    frequencies: List[int] = []
    max_frequency = 0
    max_index = 0
    for index, count in enumerate(counts):
        if count > UINT64_MAX // SCALE:
            raise UpdateError(f"{path}: {name}[{index}] overflows count * {SCALE}")
        frequency = (count * SCALE) // total
        if frequency == 0:
            frequency = 1
        frequencies.append(frequency)
        # The C implementation uses a strict comparison, selecting the first
        # bin if multiple normalized frequencies share the maximum.
        if frequency > max_frequency:
            max_frequency = frequency
            max_index = index

    adjusted = frequencies[max_index] + SCALE - sum(frequencies)
    if adjusted <= 0 or adjusted > SCALE:
        raise UpdateError(
            f"{path}: {name} normalization cannot adjust its first maximum bin"
        )
    frequencies[max_index] = adjusted
    if sum(frequencies) != SCALE:
        raise UpdateError(f"{path}: normalized {name} frequencies do not sum to {SCALE}")
    return tuple(frequencies)


def parse_report(mode: int, path: Path, min_samples: int) -> Report:
    try:
        raw = path.read_bytes()
        text = raw.decode("utf-8")
    except OSError as exc:
        raise UpdateError(f"cannot read {path}: {exc}") from exc
    except UnicodeDecodeError as exc:
        raise UpdateError(f"{path} is not valid UTF-8: {exc}") from exc
    spec = MODE_SPECS[mode]

    params = one_match(
        r"^参数:\s*N=(\d+),\s*Q=(\d+),\s*D=(\d+),\s*L=(\d+),\s*K=(\d+)\s*$",
        text,
        f"parameter line in {path}",
    )
    actual_params = tuple(int(value) for value in params.groups())
    if actual_params != spec.params:
        raise UpdateError(
            f"{path}: parameters {actual_params} do not identify mode {mode}; "
            f"expected {spec.params}"
        )

    bounds = one_match(
        r"^B0=(\d+\.\d{2}),\s*B1=(\d+\.\d{2}),\s*B=(\d+\.\d{2})\s*$",
        text,
        f"B0/B1/B line in {path}",
    )
    actual_bounds = bounds.groups()
    if actual_bounds != spec.bounds:
        raise UpdateError(
            f"{path}: B0/B1/B={actual_bounds} are stale for mode {mode}; "
            f"expected {spec.bounds}"
        )

    sample_match = one_match(r"^采样次数:\s*(\d+)\s*$", text, f"sample count in {path}")
    samples = int(sample_match.group(1))
    if samples < min_samples:
        raise UpdateError(
            f"{path}: only {samples} samples; at least {min_samples} are required "
            "(override deliberately with --min-samples)"
        )

    encoded = one_match(
        r"^可映射候选/成功编码:\s*(\d+)\s*/\s*(\d+)\s*$",
        text,
        f"mapped/encoded result in {path}",
    )
    mapped_candidates, encoded_success = (int(value) for value in encoded.groups())
    if mapped_candidates != encoded_success:
        raise UpdateError(
            f"{path}: mapped candidates / successful encodings are "
            f"{mapped_candidates}/{encoded_success}; both values must be equal"
        )

    machine_mapped = parse_decimal_field(text, "PREPACK_MAPPED_CANDIDATES", path)
    size_rejections = parse_decimal_field(text, "PREPACK_PACK_SIZE_REJECTIONS", path)
    if machine_mapped != mapped_candidates:
        raise UpdateError(
            f"{path}: PREPACK_MAPPED_CANDIDATES={machine_mapped}, but the "
            f"mapped/encoded result is {mapped_candidates}"
        )
    expected_mapped = samples + size_rejections
    if mapped_candidates != expected_mapped:
        raise UpdateError(
            f"{path}: mapped candidates={mapped_candidates}, but samples "
            f"({samples}) + PREPACK_PACK_SIZE_REJECTIONS ({size_rejections}) "
            f"is {expected_mapped}"
        )

    for name in ("PREPACK_H_COUNT", "PREPACK_HB_Z1_COUNT", "PREPACK_TOTAL_COUNT"):
        count = parse_decimal_field(text, name, path)
        if count != mapped_candidates:
            raise UpdateError(
                f"{path}: {name}={count}, expected mapped candidate count "
                f"{mapped_candidates}"
            )

    # Alphabet misses are legitimate signing retries, so their values need not
    # be zero.  Requiring both fields still catches truncated/incomplete reports.
    parse_decimal_field(text, "UNMAPPED_H", path)
    parse_decimal_field(text, "UNMAPPED_HB_Z1", path)

    for label in ("h", "hb_z1"):
        decoded = one_match(
            rf"^{label}\s+解码:\s*成功\s+(\d+),\s*失败\s+(\d+)\s*$",
            text,
            f"{label} decoding result in {path}",
        )
        successes, failures = (int(value) for value in decoded.groups())
        if successes != mapped_candidates or failures != 0:
            raise UpdateError(
                f"{path}: {label} decoding is success={successes}, failure={failures}; "
                f"expected success={mapped_candidates}, failure=0"
            )

    mapping_h_total = parse_decimal_field(text, "MAPPING_H_TOTAL", path)
    mapping_hb_z1_total = parse_decimal_field(text, "MAPPING_HB_Z1_TOTAL", path)
    counts_h = parse_count_array(text, "COUNTS_H", spec.m_h, path)
    counts_hb_z1 = parse_count_array(text, "COUNTS_HB_Z1", spec.m_hb_z1, path)

    n, _q, _d, l, k = spec.params
    expected_h_total = mapped_candidates * k * n
    expected_hb_z1_total = mapped_candidates * l * n
    if sum(counts_h) != mapping_h_total:
        raise UpdateError(
            f"{path}: COUNTS_H sums to {sum(counts_h)}, but "
            f"MAPPING_H_TOTAL={mapping_h_total}"
        )
    if mapping_h_total != expected_h_total:
        raise UpdateError(
            f"{path}: MAPPING_H_TOTAL={mapping_h_total}, expected "
            f"mapped*K*N={expected_h_total}"
        )
    if sum(counts_hb_z1) != mapping_hb_z1_total:
        raise UpdateError(
            f"{path}: COUNTS_HB_Z1 sums to {sum(counts_hb_z1)}, but "
            f"MAPPING_HB_Z1_TOTAL={mapping_hb_z1_total}"
        )
    if mapping_hb_z1_total != expected_hb_z1_total:
        raise UpdateError(
            f"{path}: MAPPING_HB_Z1_TOTAL={mapping_hb_z1_total}, expected "
            f"mapped*L*N={expected_hb_z1_total}"
        )

    f_h = parse_frequency(text, "F_H", spec.m_h, path)
    f_hb_z1 = parse_frequency(text, "F_HB_Z1", spec.m_hb_z1, path)
    recomputed_f_h = normalize_counts(counts_h, "COUNTS_H", path)
    recomputed_f_hb_z1 = normalize_counts(counts_hb_z1, "COUNTS_HB_Z1", path)
    if f_h != recomputed_f_h:
        mismatch = next(index for index in range(spec.m_h) if f_h[index] != recomputed_f_h[index])
        raise UpdateError(
            f"{path}: F_H does not match COUNTS_H normalization at index {mismatch}: "
            f"reported {f_h[mismatch]}, recomputed {recomputed_f_h[mismatch]}"
        )
    if f_hb_z1 != recomputed_f_hb_z1:
        mismatch = next(
            index
            for index in range(spec.m_hb_z1)
            if f_hb_z1[index] != recomputed_f_hb_z1[index]
        )
        raise UpdateError(
            f"{path}: F_HB_Z1 does not match COUNTS_HB_Z1 normalization at "
            f"index {mismatch}: reported {f_hb_z1[mismatch]}, "
            f"recomputed {recomputed_f_hb_z1[mismatch]}"
        )

    return Report(
        mode=mode,
        path=path,
        samples=samples,
        mapped_candidates=mapped_candidates,
        f_h=f_h,
        f_hb_z1=f_hb_z1,
        sha256=hashlib.sha256(raw).hexdigest(),
        contents=raw,
    )


def definition_pattern(name: str) -> re.Pattern[str]:
    return re.compile(
        rf"(?P<statement>^[ \t]*static[ \t]+(?:const[ \t]+)?"
        rf"(?P<ctype>RansEncSymbol|RansDecSymbol|uint16_t|uint32_t)[ \t]+"
        rf"{re.escape(name)}[ \t]*\[[^\]\r\n]+\][ \t]*=[ \t]*"
        rf"(?P<initializer>\{{.*?\}})[ \t]*;)",
        re.MULTILINE | re.DOTALL,
    )


def find_definition(text: str, name: str, where: str) -> Definition:
    matches = list(definition_pattern(name).finditer(text))
    if len(matches) != 1:
        raise UpdateError(
            f"{where}: expected exactly one definition of {name}, found {len(matches)}"
        )
    match = matches[0]
    return Definition(
        start=match.start("statement"),
        end=match.end("statement"),
        statement=match.group("statement"),
        ctype=match.group("ctype"),
        initializer=match.group("initializer"),
    )


def initializer_values(initializer: str, where: str) -> Tuple[int, ...]:
    without_comments = re.sub(r"/\*.*?\*/", "", initializer, flags=re.DOTALL)
    without_comments = re.sub(r"//[^\r\n]*", "", without_comments)
    token = r"(?<![A-Za-z0-9_])[-+]?\d+(?![A-Za-z0-9_])"
    values = tuple(int(value) for value in re.findall(token, without_comments))
    leftover = re.sub(token, "", without_comments)
    if re.sub(r"[{},\s]", "", leftover):
        raise UpdateError(f"{where}: initializer contains unsupported non-integer syntax")
    return values


def mode_branch(text: str, macro: str, mode: int, where: str) -> Tuple[int, int]:
    start_pattern = re.compile(
        rf"^[ \t]*#(?:if|elif)[ \t]+{re.escape(macro)}[ \t]*==[ \t]*{mode}[ \t\r]*$",
        re.MULTILINE,
    )
    starts = list(start_pattern.finditer(text))
    if len(starts) != 1:
        raise UpdateError(
            f"{where}: expected one {macro} == {mode} branch, found {len(starts)}"
        )
    start = starts[0].end()
    end_pattern = re.compile(
        rf"^[ \t]*#(?:elif[ \t]+{re.escape(macro)}[ \t]*==|endif\b)", re.MULTILINE
    )
    end = end_pattern.search(text, start)
    if end is None:
        raise UpdateError(f"{where}: unterminated {macro} == {mode} branch")
    return start, end.start()


def newline_for(text: str) -> str:
    return "\r\n" if text.count("\r\n") > text.count("\n") / 2 else "\n"


def format_frequency(name: str, values: Sequence[int], newline: str) -> str:
    size_macro = "M_H" if name == "f_h" else "M_HB_Z1"
    lines = [f"static const uint32_t {name}[{size_macro}] = {{"]
    for offset in range(0, len(values), 12):
        chunk = ", ".join(str(value) for value in values[offset : offset + 12])
        suffix = "," if offset + 12 < len(values) else ""
        lines.append(f"    {chunk}{suffix}")
    lines.append("};")
    return newline.join(lines)


def update_frequency_source(
    text: str, reports: Mapping[int, Report], where: str
) -> Tuple[str, List[str]]:
    changes: List[str] = []
    for mode in MODES:
        begin, end = mode_branch(text, "DARTS_MODE", mode, where)
        branch = text[begin:end]
        replacements: List[Tuple[int, int, str]] = []
        wanted = {
            "f_h": reports[mode].f_h,
            "f_hb_z1": reports[mode].f_hb_z1,
        }
        for name, values in wanted.items():
            definition = find_definition(branch, name, f"{where}, mode {mode}")
            if definition.ctype != "uint32_t":
                raise UpdateError(
                    f"{where}, mode {mode}: {name} has type {definition.ctype}, expected uint32_t"
                )
            current = initializer_values(
                definition.initializer, f"{where}, mode {mode}, {name}"
            )
            if current != values:
                replacement = format_frequency(name, values, newline_for(text))
                replacements.append((definition.start, definition.end, replacement))
                changes.append(f"mode {mode} {name}")
        for start, stop, replacement in sorted(replacements, reverse=True):
            branch = branch[:start] + replacement + branch[stop:]
        text = text[:begin] + branch + text[end:]
    return text, changes


def expected_type(name: str) -> str:
    if name.startswith("esyms_"):
        return "RansEncSymbol"
    if name.startswith("dsyms_"):
        return "RansDecSymbol"
    return "uint16_t"


def normalized_table_values(
    name: str, values: Sequence[int], mode: int, where: str
) -> Tuple[int, ...]:
    spec = MODE_SPECS[mode]
    alphabet_size = spec.m_hb_z1 if name.endswith("hb_z1") else spec.m_h
    if name.startswith("esyms_"):
        expected = alphabet_size * 5
        if len(values) != expected:
            raise UpdateError(f"{where}: {name} has {len(values)} fields, expected {expected}")
        normalized: List[int] = []
        for offset in range(0, len(values), 5):
            normalized.extend(value & 0xFFFFFFFF for value in values[offset : offset + 3])
            normalized.extend(value & 0xFFFF for value in values[offset + 3 : offset + 5])
        return tuple(normalized)
    if name.startswith("dsyms_"):
        expected = alphabet_size * 2
        if len(values) != expected:
            raise UpdateError(f"{where}: {name} has {len(values)} fields, expected {expected}")
        return tuple(value & 0xFFFF for value in values)
    if len(values) != SCALE:
        raise UpdateError(f"{where}: {name} has {len(values)} entries, expected {SCALE}")
    return tuple(value & 0xFFFF for value in values)


def parse_generated_tables(output: str, mode: int) -> Dict[str, GeneratedTable]:
    tables: Dict[str, GeneratedTable] = {}
    for name in TABLE_NAMES:
        definition = find_definition(output, name, f"generator output for mode {mode}")
        wanted_type = expected_type(name)
        if definition.ctype != wanted_type:
            raise UpdateError(
                f"generator output for mode {mode}: {name} has type "
                f"{definition.ctype}, expected {wanted_type}"
            )
        values = initializer_values(
            definition.initializer, f"generator output for mode {mode}, {name}"
        )
        normalized_table_values(name, values, mode, f"generator output for mode {mode}")
        tables[name] = GeneratedTable(definition.statement, definition.ctype, values)
    return tables


def run_generator(
    canonical_source: str, header: Path, compiler: Sequence[str]
) -> Dict[int, Dict[str, GeneratedTable]]:
    if not header.is_file():
        raise UpdateError(f"missing rANS header: {header}")
    generated: Dict[int, Dict[str, GeneratedTable]] = {}
    with tempfile.TemporaryDirectory(prefix="darts-rans-tables-") as temp_name:
        temp_root = Path(temp_name)
        source_dir = temp_root / "precomputations_rANS"
        source_dir.mkdir()
        source = source_dir / "prec_encoding.c"
        source.write_bytes(canonical_source.encode("utf-8"))
        shutil.copy2(header, temp_root / "rans_byte.h")
        for mode in MODES:
            executable = temp_root / f"prec_encoding_{mode}"
            command = list(compiler) + [
                "-O2",
                "-std=gnu99",
                f"-DDARTS_MODE={mode}",
                str(source),
                "-o",
                str(executable),
            ]
            try:
                compiled = subprocess.run(command, capture_output=True, text=True, check=False)
            except OSError as exc:
                raise UpdateError(f"cannot run compiler {compiler[0]!r}: {exc}") from exc
            if compiled.returncode != 0:
                raise UpdateError(
                    f"generator compile failed for mode {mode}:\n{compiled.stderr.strip()}"
                )
            ran = subprocess.run([str(executable)], capture_output=True, text=True, check=False)
            if ran.returncode != 0:
                raise UpdateError(
                    f"generator failed for mode {mode}:\n{ran.stderr.strip()}"
                )
            generated[mode] = parse_generated_tables(ran.stdout, mode)
    return generated


def encoding_macro(text: str, where: str) -> str:
    found = [
        macro
        for macro in ("DARTS_MODE", "SIGN_MODE")
        if re.search(rf"^[ \t]*#if[ \t]+{macro}[ \t]*==", text, re.MULTILINE)
    ]
    if len(found) != 1:
        raise UpdateError(
            f"{where}: expected exactly one mode macro family, found {', '.join(found) or 'none'}"
        )
    return found[0]


def update_encoding_source(
    text: str, generated: Mapping[int, Mapping[str, GeneratedTable]], where: str
) -> Tuple[str, List[str]]:
    macro = encoding_macro(text, where)
    changes: List[str] = []
    for mode in MODES:
        begin, end = mode_branch(text, macro, mode, where)
        branch = text[begin:end]
        replacements: List[Tuple[int, int, str]] = []
        for name in TABLE_NAMES:
            current_def = find_definition(branch, name, f"{where}, mode {mode}")
            wanted = generated[mode][name]
            if current_def.ctype != wanted.ctype:
                raise UpdateError(
                    f"{where}, mode {mode}: {name} has type {current_def.ctype}, "
                    f"expected {wanted.ctype}"
                )
            current_values = initializer_values(
                current_def.initializer, f"{where}, mode {mode}, {name}"
            )
            current_normalized = normalized_table_values(
                name, current_values, mode, f"{where}, mode {mode}"
            )
            wanted_normalized = normalized_table_values(
                name, wanted.values, mode, f"generated mode {mode}"
            )
            if current_normalized != wanted_normalized:
                replacements.append((current_def.start, current_def.end, wanted.statement))
                changes.append(f"mode {mode} {name}")
        for start, stop, replacement in sorted(replacements, reverse=True):
            branch = branch[:start] + replacement + branch[stop:]
        text = text[:begin] + branch + text[end:]
    return text, changes


def discover_sources(root: Path) -> Tuple[List[Path], List[Path]]:
    if not root.is_dir():
        raise UpdateError(f"implementations root is not a directory: {root}")
    prec_paths = sorted(root.rglob("precomputations_rANS/prec_encoding.c"))
    encoding_paths = sorted(root.rglob("encoding.c"))
    if len(prec_paths) != EXPECTED_IMPLEMENTATIONS:
        raise UpdateError(
            f"found {len(prec_paths)} prec_encoding.c files under {root}, "
            f"expected {EXPECTED_IMPLEMENTATIONS}"
        )
    if len(encoding_paths) != EXPECTED_IMPLEMENTATIONS:
        raise UpdateError(
            f"found {len(encoding_paths)} encoding.c files under {root}, "
            f"expected {EXPECTED_IMPLEMENTATIONS}"
        )
    prec_roots = {path.parent.parent.resolve() for path in prec_paths}
    encoding_roots = {path.parent.resolve() for path in encoding_paths}
    if prec_roots != encoding_roots:
        missing_encoding = sorted(str(path) for path in prec_roots - encoding_roots)
        missing_prec = sorted(str(path) for path in encoding_roots - prec_roots)
        raise UpdateError(
            "prec_encoding.c/encoding.c implementation roots do not match; "
            f"missing encoding roots={missing_encoding}, missing prec roots={missing_prec}"
        )
    return prec_paths, encoding_paths


def parse_report_overrides(items: Iterable[str]) -> Dict[int, Path]:
    overrides: Dict[int, Path] = {}
    for item in items:
        if "=" not in item:
            raise UpdateError(f"invalid --report {item!r}; expected MODE=PATH")
        mode_text, path_text = item.split("=", 1)
        try:
            mode = int(mode_text)
        except ValueError as exc:
            raise UpdateError(f"invalid report mode {mode_text!r}") from exc
        if mode not in MODES:
            raise UpdateError(f"unsupported report mode {mode}; choose 128, 256, or 512")
        if mode in overrides:
            raise UpdateError(f"mode {mode} was supplied more than once with --report")
        if not path_text:
            raise UpdateError(f"mode {mode} has an empty report path")
        overrides[mode] = Path(path_text).expanduser().resolve()
    return overrides


def atomic_write_bytes(path: Path, data: bytes) -> None:
    mode = stat.S_IMODE(path.stat().st_mode) if path.exists() else 0o644
    descriptor, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=str(path.parent))
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(temp_name, mode)
        os.replace(temp_name, path)
    except Exception:
        try:
            os.unlink(temp_name)
        except FileNotFoundError:
            pass
        raise


def atomic_write(path: Path, text: str) -> None:
    atomic_write_bytes(path, text.encode("utf-8"))


def build_parser(default_root: Path) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Check or synchronize every DARTS rANS frequency/generated-table copy."
    )
    parser.add_argument(
        "--report",
        action="append",
        default=[],
        metavar="MODE=PATH",
        help="override one report path; repeat for modes 128, 256, and 512",
    )
    parser.add_argument(
        "--reports-dir",
        type=Path,
        help="directory containing rANS_128.txt, rANS_256.txt, and rANS_512.txt",
    )
    parser.add_argument(
        "--implementations-root",
        type=Path,
        default=default_root,
        help=f"DARTS Implementations directory (default: {default_root})",
    )
    parser.add_argument(
        "--min-samples",
        type=int,
        default=200000,
        help="minimum declared successful samples per report (default: 200000)",
    )
    parser.add_argument(
        "--cc",
        default=os.environ.get("CC", "cc"),
        help="C compiler command used for prec_encoding.c (default: CC or cc)",
    )
    action = parser.add_mutually_exclusive_group()
    action.add_argument(
        "--check",
        action="store_true",
        help="validate and report drift without writing (the default)",
    )
    action.add_argument(
        "--write",
        action="store_true",
        help="atomically replace source files whose rANS tables are stale",
    )
    parser.add_argument(
        "--sync-reports",
        action="store_true",
        help=(
            "include all 45 report replicas in the check; with --write, copy "
            "each complete input report to all 15 precomputations_rANS directories"
        ),
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    default_root = Path(__file__).resolve().parents[3]
    parser = build_parser(default_root)
    args = parser.parse_args(argv)
    try:
        if args.min_samples <= 0:
            raise UpdateError("--min-samples must be positive")
        root = args.implementations_root.expanduser().resolve()
        canonical_dir = root / "Reference_Implementation" / "DARTS128"
        canonical_prec = canonical_dir / "precomputations_rANS" / "prec_encoding.c"
        header = canonical_dir / "rans_byte.h"

        report_dir = (
            args.reports_dir.expanduser().resolve()
            if args.reports_dir
            else canonical_prec.parent
        )
        report_paths = {mode: report_dir / f"rANS_{mode}.txt" for mode in MODES}
        report_paths.update(parse_report_overrides(args.report))
        missing = [f"mode {mode}: {path}" for mode, path in report_paths.items() if not path.is_file()]
        if missing:
            raise UpdateError("missing report file(s): " + "; ".join(missing))

        reports = {
            mode: parse_report(mode, report_paths[mode], args.min_samples) for mode in MODES
        }
        for mode in MODES:
            report = reports[mode]
            print(
                f"report {mode}: samples={report.samples}, "
                f"mapped={report.mapped_candidates}, "
                f"F_H={len(report.f_h)}/{sum(report.f_h)}, "
                f"F_HB_Z1={len(report.f_hb_z1)}/{sum(report.f_hb_z1)}, "
                f"sha256={report.sha256}"
            )

        prec_paths, encoding_paths = discover_sources(root)
        if canonical_prec not in prec_paths:
            raise UpdateError(f"canonical generator is outside the discovered source set: {canonical_prec}")

        updated: Dict[Path, str] = {}
        details: Dict[Path, List[str]] = {}
        prec_texts: Dict[Path, str] = {}
        for path in prec_paths:
            original = read_text(path)
            transformed, changes = update_frequency_source(original, reports, str(path))
            prec_texts[path] = transformed
            if transformed != original:
                updated[path] = transformed
                details[path] = changes

        compiler = shlex.split(args.cc)
        if not compiler:
            raise UpdateError("--cc resolves to an empty compiler command")
        generated = run_generator(prec_texts[canonical_prec], header, compiler)

        for path in encoding_paths:
            original = read_text(path)
            transformed, changes = update_encoding_source(original, generated, str(path))
            if transformed != original:
                updated[path] = transformed
                details[path] = changes

        report_updates: Dict[Path, bytes] = {}
        if args.sync_reports:
            report_bytes = {mode: reports[mode].contents for mode in MODES}
            report_targets = [
                path.parent / f"rANS_{mode}.txt"
                for path in prec_paths
                for mode in MODES
            ]
            if len(report_targets) != EXPECTED_IMPLEMENTATIONS * len(MODES):
                raise UpdateError(
                    f"computed {len(report_targets)} report targets, expected "
                    f"{EXPECTED_IMPLEMENTATIONS * len(MODES)}"
                )
            for target in report_targets:
                mode_match = re.fullmatch(r"rANS_(128|256|512)\.txt", target.name)
                if mode_match is None:
                    raise UpdateError(f"internal error: unexpected report target {target}")
                wanted_bytes = report_bytes[int(mode_match.group(1))]
                current_bytes = target.read_bytes() if target.exists() else None
                if current_bytes != wanted_bytes:
                    report_updates[target] = wanted_bytes

        print(
            f"validated {len(prec_paths)} prec_encoding.c and "
            f"{len(encoding_paths)} encoding.c files; generated six tables for each of 3 modes"
        )
        if not updated and not report_updates:
            print("OK: all rANS frequency and generated tables match the reports")
            if args.sync_reports:
                print("OK: all 45 checked-in rANS reports match the three input reports")
            return 0

        for path in sorted(updated):
            relative = path.relative_to(root)
            print(f"STALE: {relative}: {', '.join(details[path])}")
        for path in sorted(report_updates):
            print(f"STALE REPORT: {path.relative_to(root)}")
        if not args.write:
            print(
                f"CHECK FAILED: {len(updated)} source files and "
                f"{len(report_updates)} report files need synchronization; "
                "rerun with --write"
            )
            return 1

        # Every report, source, generated table, and replacement has been validated
        # before the first checked-in file is touched.
        for path in sorted(updated):
            atomic_write(path, updated[path])
        for path in sorted(report_updates):
            atomic_write_bytes(path, report_updates[path])
        print(
            f"UPDATED: synchronized {len(updated)} source files and "
            f"{len(report_updates)} report files"
        )
        return 0
    except UpdateError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    except OSError as exc:
        print(f"ERROR: filesystem operation failed: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
