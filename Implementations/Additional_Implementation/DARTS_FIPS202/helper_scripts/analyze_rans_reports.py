#!/usr/bin/env python3
"""Validate and summarize machine-readable ``test_rANS`` reports.

Example::

    python3 analyze_rans_reports.py \
        --report 128=/tmp/rans-128.txt \
        --report 256=/tmp/rans-256.txt \
        --report 512=/tmp/rans-512.txt

The script is deliberately read-only.  The current rANS bases, size cap, and
fixed part of the signature are read from the reference ``params.h`` so that
the recommendations track parameter updates made elsewhere in the tree.
"""

from __future__ import annotations

import argparse
import ast
import math
import operator
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


MODES = (128, 256, 512)
KINDS = ("H", "HB_Z1", "TOTAL")
QUANTILES = (
    ("Q99", 99, 100),
    ("Q99.9", 999, 1000),
    ("Q99.99", 9999, 10000),
)


class ReportError(ValueError):
    """An input report or local parameter file is inconsistent."""


@dataclass(frozen=True)
class ModeParameters:
    base_h: int
    base_hb_z1: int
    cap: int
    fixed_part: int


@dataclass
class ParsedReport:
    path: Path
    counts: Dict[str, int] = field(default_factory=dict)
    histograms: Dict[str, Dict[int, int]] = field(
        default_factory=lambda: {kind: {} for kind in KINDS}
    )
    reported_quantiles: Dict[Tuple[str, str], int] = field(default_factory=dict)
    mapped_candidates: int | None = None
    pack_size_rejections: int | None = None
    step9: Dict[int, Tuple[int, int]] = field(default_factory=dict)
    final_b: Tuple[int, int] | None = None


_BINOPS = {
    ast.Add: operator.add,
    ast.Sub: operator.sub,
    ast.Mult: operator.mul,
    ast.Div: operator.floordiv,
    ast.FloorDiv: operator.floordiv,
    ast.LShift: operator.lshift,
    ast.RShift: operator.rshift,
}
_UNARYOPS = {ast.UAdd: operator.pos, ast.USub: operator.neg}


def _eval_integer_expression(expression: str, macros: Mapping[str, str]) -> int:
    """Evaluate the small integer-only subset used by the relevant macros."""

    try:
        root = ast.parse(expression.strip(), mode="eval").body
    except SyntaxError as exc:
        raise ReportError(f"cannot parse parameter expression {expression!r}") from exc

    resolving: set[str] = set()

    def visit(node: ast.AST) -> int:
        if isinstance(node, ast.Constant) and isinstance(node.value, int):
            return node.value
        if isinstance(node, ast.Name):
            if node.id not in macros:
                raise ReportError(f"unknown macro {node.id!r} in {expression!r}")
            if node.id in resolving:
                raise ReportError(f"recursive macro {node.id!r}")
            resolving.add(node.id)
            try:
                return _eval_integer_expression(macros[node.id], macros)
            finally:
                resolving.remove(node.id)
        if isinstance(node, ast.BinOp) and type(node.op) in _BINOPS:
            return int(_BINOPS[type(node.op)](visit(node.left), visit(node.right)))
        if isinstance(node, ast.UnaryOp) and type(node.op) in _UNARYOPS:
            return int(_UNARYOPS[type(node.op)](visit(node.operand)))
        raise ReportError(f"unsupported parameter expression {expression!r}")

    return visit(root)


def _mode_block(params_text: str, mode: int) -> str:
    start_re = re.compile(
        rf"^\s*#if\s+DARTS_MODE\s*==\s*{mode}\s*$", re.MULTILINE
    )
    found_mode_block = False
    for start in start_re.finditer(params_text):
        found_mode_block = True
        end = re.search(
            r"^\s*#endif\s*$", params_text[start.end() :], re.MULTILINE
        )
        if end is None:
            raise ReportError(f"unterminated DARTS_MODE == {mode} block in params.h")
        block = params_text[start.end() : start.end() + end.start()]
        # Mode 512 also has a short MUBYTES-only conditional near the top.
        if re.search(r"^\s*#define\s+BASE_ENC_H\b", block, re.MULTILINE):
            return block
    if not found_mode_block:
        raise ReportError(f"params.h has no DARTS_MODE == {mode} block")
    raise ReportError(f"DARTS_MODE == {mode} has no rANS parameter block")


def load_parameters(script_path: Path) -> Dict[int, ModeParameters]:
    darts_root = script_path.resolve().parents[4]
    params_path = (
        darts_root
        / "Implementations"
        / "Reference_Implementation"
        / "DARTS128"
        / "params.h"
    )
    try:
        params_text = params_path.read_text(encoding="utf-8")
    except OSError as exc:
        raise ReportError(f"cannot read current parameters: {params_path}: {exc}") from exc

    result: Dict[int, ModeParameters] = {}
    for mode in MODES:
        block = _mode_block(params_text, mode)
        macros = {
            match.group(1): match.group(2).strip()
            for match in re.finditer(
                r"^\s*#define\s+([A-Za-z_]\w*)\s+(.+?)\s*(?://.*)?$",
                block,
                re.MULTILINE,
            )
        }

        needed = (
            "BASE_ENC_H",
            "BASE_ENC_HB_Z1",
            "ENCODE_MAXSIZE",
            "POLY_C_PACKEDBYTES",
            "POLY_LOWBITS_PACKEDBYTES",
            "L",
        )
        missing = [name for name in needed if name not in macros]
        if missing:
            raise ReportError(
                f"DARTS_MODE {mode} is missing macros: {', '.join(missing)}"
            )

        values = {
            name: _eval_integer_expression(macros[name], macros) for name in needed
        }
        fixed_part = (
            values["POLY_C_PACKEDBYTES"]
            + values["L"] * values["POLY_LOWBITS_PACKEDBYTES"]
            + 2
        )
        result[mode] = ModeParameters(
            base_h=values["BASE_ENC_H"],
            base_hb_z1=values["BASE_ENC_HB_Z1"],
            cap=values["ENCODE_MAXSIZE"],
            fixed_part=fixed_part,
        )
    return result


def _set_once(
    current: int | Tuple[int, int] | None,
    value: int | Tuple[int, int],
    label: str,
    path: Path,
) -> int | Tuple[int, int]:
    if current is not None:
        raise ReportError(f"{path}: duplicate {label}")
    return value


def parse_report(path: Path) -> ParsedReport:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise ReportError(f"cannot read report {path}: {exc}") from exc

    report = ParsedReport(path=path)
    hist_re = re.compile(r"^PREPACK_(H|HB_Z1|TOTAL)_HIST=(\d+),(\d+)$")
    count_re = re.compile(r"^PREPACK_(H|HB_Z1|TOTAL)_COUNT=(\d+)$")
    quantile_re = re.compile(
        r"^PREPACK_(H|HB_Z1|TOTAL)_(Q50|Q90|Q95|Q99|Q999|Q9999|Q100)=(\d+)$"
    )
    step9_re = re.compile(r"^\s*b=([01]):\s*(\d+)/(\d+)(?:\s*\([^)]*\))?\s*$")
    final_re = re.compile(r"^最终 b=0/b=1:\s*(\d+)/(\d+)\s*$")

    for line_number, line in enumerate(lines, start=1):
        match = hist_re.fullmatch(line)
        if match:
            kind, size_text, frequency_text = match.groups()
            size, frequency = int(size_text), int(frequency_text)
            if frequency <= 0:
                raise ReportError(
                    f"{path}:{line_number}: histogram frequency must be positive"
                )
            if size in report.histograms[kind]:
                raise ReportError(
                    f"{path}:{line_number}: duplicate PREPACK_{kind}_HIST size {size}"
                )
            report.histograms[kind][size] = frequency
            continue

        match = count_re.fullmatch(line)
        if match:
            kind, count_text = match.groups()
            if kind in report.counts:
                raise ReportError(f"{path}:{line_number}: duplicate PREPACK_{kind}_COUNT")
            report.counts[kind] = int(count_text)
            continue

        match = quantile_re.fullmatch(line)
        if match:
            kind, label, value_text = match.groups()
            key = (kind, label)
            if key in report.reported_quantiles:
                raise ReportError(
                    f"{path}:{line_number}: duplicate PREPACK_{kind}_{label}"
                )
            report.reported_quantiles[key] = int(value_text)
            continue

        if line.startswith("PREPACK_MAPPED_CANDIDATES="):
            value = int(line.partition("=")[2])
            report.mapped_candidates = int(
                _set_once(
                    report.mapped_candidates, value, "PREPACK_MAPPED_CANDIDATES", path
                )
            )
            continue

        if line.startswith("PREPACK_PACK_SIZE_REJECTIONS="):
            value = int(line.partition("=")[2])
            report.pack_size_rejections = int(
                _set_once(
                    report.pack_size_rejections,
                    value,
                    "PREPACK_PACK_SIZE_REJECTIONS",
                    path,
                )
            )
            continue

        match = step9_re.fullmatch(line)
        if match:
            bit, rejects, attempts = (int(value) for value in match.groups())
            if bit in report.step9:
                raise ReportError(f"{path}:{line_number}: duplicate Step9 b={bit}")
            report.step9[bit] = (rejects, attempts)
            continue

        match = final_re.fullmatch(line)
        if match:
            value = (int(match.group(1)), int(match.group(2)))
            report.final_b = tuple(
                _set_once(report.final_b, value, "final b=0/b=1", path)
            )  # type: ignore[arg-type]

    return report


def histogram_count(histogram: Mapping[int, int]) -> int:
    return sum(histogram.values())


def histogram_sum(histogram: Mapping[int, int]) -> int:
    return sum(size * frequency for size, frequency in histogram.items())


def histogram_quantile(
    histogram: Mapping[int, int], numerator: int, denominator: int
) -> int:
    total = histogram_count(histogram)
    if total <= 0:
        raise ReportError("cannot compute a quantile of an empty histogram")
    target = (total * numerator + denominator - 1) // denominator
    cumulative = 0
    for size, frequency in sorted(histogram.items()):
        cumulative += frequency
        if cumulative >= target:
            return size
    raise AssertionError("histogram count changed while computing a quantile")


def validate_report(report: ParsedReport, params: ModeParameters) -> List[str]:
    errors: List[str] = []
    for kind in KINDS:
        if kind not in report.counts:
            errors.append(f"missing PREPACK_{kind}_COUNT")
        if not report.histograms[kind]:
            errors.append(f"missing PREPACK_{kind}_HIST")

    if errors:
        return errors

    common_count = report.counts["TOTAL"]
    for kind in KINDS:
        declared = report.counts[kind]
        observed = histogram_count(report.histograms[kind])
        if declared != observed:
            errors.append(
                f"PREPACK_{kind}_COUNT={declared}, but its HIST sum is {observed}"
            )
        if declared != common_count:
            errors.append(
                f"PREPACK_{kind}_COUNT={declared}, expected TOTAL_COUNT={common_count}"
            )

    if report.mapped_candidates is None:
        errors.append("missing PREPACK_MAPPED_CANDIDATES")
    elif report.mapped_candidates != common_count:
        errors.append(
            f"PREPACK_MAPPED_CANDIDATES={report.mapped_candidates}, "
            f"expected {common_count}"
        )

    component_sum = histogram_sum(report.histograms["H"]) + histogram_sum(
        report.histograms["HB_Z1"]
    )
    total_sum = histogram_sum(report.histograms["TOTAL"])
    if total_sum != component_sum:
        errors.append(
            f"TOTAL HIST weighted sum is {total_sum}, but H + HB_Z1 is {component_sum}"
        )

    min_component_total = min(report.histograms["H"]) + min(
        report.histograms["HB_Z1"]
    )
    max_component_total = max(report.histograms["H"]) + max(
        report.histograms["HB_Z1"]
    )
    outside_total = sum(
        frequency
        for size, frequency in report.histograms["TOTAL"].items()
        if size < min_component_total or size > max_component_total
    )
    if outside_total:
        errors.append(
            f"TOTAL HIST has {outside_total} samples outside the possible marginal "
            f"range [{min_component_total}, {max_component_total}]"
        )

    printed_labels = (("Q99", 99, 100), ("Q999", 999, 1000), ("Q9999", 9999, 10000))
    for kind in KINDS:
        for label, numerator, denominator in printed_labels:
            key = (kind, label)
            if key not in report.reported_quantiles:
                errors.append(f"missing PREPACK_{kind}_{label}")
                continue
            computed = histogram_quantile(
                report.histograms[kind], numerator, denominator
            )
            if report.reported_quantiles[key] != computed:
                errors.append(
                    f"PREPACK_{kind}_{label}={report.reported_quantiles[key]}, "
                    f"but HIST gives {computed}"
                )

    tail_count = sum(
        frequency
        for size, frequency in report.histograms["TOTAL"].items()
        if size > params.cap
    )
    if report.pack_size_rejections is None:
        errors.append("missing PREPACK_PACK_SIZE_REJECTIONS")
    elif report.pack_size_rejections != tail_count:
        errors.append(
            f"PREPACK_PACK_SIZE_REJECTIONS={report.pack_size_rejections}, "
            f"but TOTAL HIST has {tail_count} samples above cap {params.cap}"
        )

    for bit in (0, 1):
        if bit not in report.step9:
            errors.append(f"missing Step9 b={bit} attempts/rejects")
        else:
            rejects, attempts = report.step9[bit]
            if attempts <= 0:
                errors.append(f"Step9 b={bit} attempts must be positive")
            if rejects > attempts:
                errors.append(
                    f"Step9 b={bit} rejects {rejects} exceed attempts {attempts}"
                )
    if report.final_b is None:
        errors.append("missing final b=0/b=1 counts")
    elif sum(report.final_b) <= 0:
        errors.append("final b=0/b=1 count is zero")

    return errors


def _rate(numerator: int, denominator: int) -> float:
    return numerator / denominator if denominator else math.nan


def _format_z(value: float) -> str:
    return "undefined" if not math.isfinite(value) else f"{value:+.4f}"


def print_summary(mode: int, report: ParsedReport, params: ModeParameters) -> None:
    print(f"\nDARTS-{mode}: {report.path}")
    print("  validation: OK")
    for kind in KINDS:
        histogram = report.histograms[kind]
        count = histogram_count(histogram)
        mean = histogram_sum(histogram) / count
        quantile_values = [
            histogram_quantile(histogram, numerator, denominator)
            for _, numerator, denominator in QUANTILES
        ]
        print(
            f"  {kind:5s}: n={count}, min={min(histogram)}, max={max(histogram)}, "
            f"mean={mean:.4f}, Q99={quantile_values[0]}, "
            f"Q99.9={quantile_values[1]}, Q99.99={quantile_values[2]}"
        )

    total_hist = report.histograms["TOTAL"]
    total_count = histogram_count(total_hist)
    tail_count = sum(
        frequency
        for size, frequency in total_hist.items()
        if size > params.cap
    )
    suggested_cap = histogram_quantile(total_hist, 999, 1000)
    print(
        f"  current cap={params.cap}: tail={tail_count}/{total_count} "
        f"({_rate(tail_count, total_count):.6%})"
    )
    print(
        f"  empirical minimum cap for tail <= 0.1%: {suggested_cap} "
        f"(empirical Q99.9)"
    )

    for kind, base in (("H", params.base_h), ("HB_Z1", params.base_hb_z1)):
        histogram = report.histograms[kind]
        upper = base + 255
        below = sum(freq for size, freq in histogram.items() if size < base)
        above = sum(freq for size, freq in histogram.items() if size > upper)
        print(
            f"  {kind:5s} byte range [{base}, {upper}]: "
            f"out-of-range={below + above} (below={below}, above={above})"
        )

    rejects0, attempts0 = report.step9[0]
    rejects1, attempts1 = report.step9[1]
    rate0, rate1 = _rate(rejects0, attempts0), _rate(rejects1, attempts1)
    pooled = _rate(rejects0 + rejects1, attempts0 + attempts1)
    variance = pooled * (1.0 - pooled) * (1.0 / attempts0 + 1.0 / attempts1)
    step9_z = (rate0 - rate1) / math.sqrt(variance) if variance > 0.0 else math.nan
    print(
        f"  Step9: b=0 {rejects0}/{attempts0} ({rate0:.6%}), "
        f"b=1 {rejects1}/{attempts1} ({rate1:.6%})"
    )
    print(
        f"         rejection-rate difference (b0-b1)={rate0 - rate1:+.6%}, "
        f"two-proportion z={_format_z(step9_z)}"
    )

    final0, final1 = report.final_b  # type: ignore[misc]
    final_total = final0 + final1
    final_rate0 = final0 / final_total
    balance_z = (final0 - final1) / math.sqrt(final_total)
    print(
        f"  final b balance: b=0 {final0}/{final_total} ({final_rate0:.6%}), "
        f"b=1 {final1}/{final_total} ({1.0 - final_rate0:.6%}), "
        f"delta-from-50%={final_rate0 - 0.5:+.6%}, z={balance_z:+.4f}"
    )
    print(
        f"  suggested CRYPTO_SIGNATUREBYTES={params.fixed_part + suggested_cap} "
        f"(fixed {params.fixed_part} + cap {suggested_cap})"
    )


def parse_report_arguments(values: Iterable[str]) -> Dict[int, Path]:
    reports: Dict[int, Path] = {}
    for value in values:
        mode_text, separator, path_text = value.partition("=")
        if not separator or not mode_text or not path_text:
            raise ReportError(f"invalid --report {value!r}; expected MODE=PATH")
        try:
            mode = int(mode_text)
        except ValueError as exc:
            raise ReportError(f"invalid report mode {mode_text!r}") from exc
        if mode not in MODES:
            raise ReportError(f"unsupported report mode {mode}; expected 128, 256, or 512")
        if mode in reports:
            raise ReportError(f"duplicate report for mode {mode}")
        reports[mode] = Path(path_text)

    missing = [str(mode) for mode in MODES if mode not in reports]
    if missing:
        raise ReportError(f"missing --report for mode(s): {', '.join(missing)}")
    return reports


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate and summarize DARTS test_rANS machine-readable reports."
    )
    parser.add_argument(
        "--report",
        action="append",
        required=True,
        metavar="MODE=PATH",
        help="report for DARTS mode 128, 256, or 512 (repeat once per mode)",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        paths = parse_report_arguments(args.report)
        parameters = load_parameters(Path(__file__))
        reports = {mode: parse_report(paths[mode]) for mode in MODES}
    except (ReportError, ValueError) as exc:
        parser.error(str(exc))

    validation_failed = False
    for mode in MODES:
        errors = validate_report(reports[mode], parameters[mode])
        if errors:
            validation_failed = True
            print(f"DARTS-{mode}: {reports[mode].path}", file=sys.stderr)
            for error in errors:
                print(f"  ERROR: {error}", file=sys.stderr)

    if validation_failed:
        return 1

    print("All report consistency checks passed.")
    for mode in MODES:
        print_summary(mode, reports[mode], parameters[mode])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
