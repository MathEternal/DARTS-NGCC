#!/usr/bin/env python3
"""Generate and audit the DARTS-Q80 Gaussian-sampler constants.

The base sampler is the half Gaussian D_{Z+,16}, whose probability mass is

    p(k) = exp(-k^2 / 512) / sum_{j >= 0} exp(-j^2 / 512).

For k >= 1 the 112-bit implementation probability is rounded down.  All
rounding residue (including the unrepresented tail) is assigned to k = 0.
This is the one-sided base-sampler construction used by the Renyi-divergence
argument; it is deliberately not a rounded CDF.

The exp(-z) approximation uses the degree-12 Chebyshev truncation on
[0, 391/512], lifted by 2^-62.  The upper endpoint is (2*195+1)/512,
corresponding to the largest base-CDT sample and 72-bit tail.  Coefficients and every Horner multiplication
are rounded upward at Q63 precision.  The script prints a conservative
analytic relative-error bound as well as an optional dense implementation
scan.  The analytic proof relies on the documented Bessel-tail bound below;
the dense scan is only a regression check.

The production input is not exact: truncation in fixpoint_square followed by
Q63 rounding contributes at most 0.500061 Q63 ulp.  Both the one-sided
ceiling check and the reported relative-error bound include that perturbation.

Requires only Python 3 and mpmath.  Constants are computed twice at different
precisions and must agree before they are emitted.
"""

from __future__ import annotations

import argparse
import sys
from typing import Iterable, List, Sequence, Tuple

import mpmath as mp


CDT_BITS = 112
EXP_BITS = 63
EXP_DEGREE = 12
EXP_ZMAX_NUM = 391
EXP_ZMAX_DEN = 512
EXP_LIFT_BITS = 62
EXP_INPUT_ERROR_ULPS = mp.mpf("0.500061")
RENYI_ORDERS = (255, 511, 1023)
CDT_LOG2_RENYI_TARGET = mp.mpf("-105.5")
EXP_K_TARGET = mp.mpf("58.0")


def _half_gaussian_probabilities(dps: int) -> List[mp.mpf]:
    """Return probabilities through the first term rounding to zero at Q112."""
    with mp.workdps(dps):
        # k=2047 is vastly beyond what is needed at this precision.  The fixed
        # upper bound makes generation deterministic across mpmath versions.
        weights = [mp.exp(-mp.mpf(k * k) / 512) for k in range(2048)]
        z_half = mp.fsum(weights)
        probs = [w / z_half for w in weights]
        scale = 1 << CDT_BITS
        end = 1
        while int(mp.floor(probs[end] * scale)) != 0:
            end += 1
        return [+p for p in probs[: end + 1]]


def generate_cdt(dps: int) -> Tuple[List[int], List[int], List[mp.mpf]]:
    """Return thresholds, exact integer counts, and ideal probabilities."""
    probs = _half_gaussian_probabilities(dps)
    scale = 1 << CDT_BITS
    positive_counts: List[int] = []
    for p in probs[1:]:
        count = int(mp.floor(p * scale))
        if count == 0:
            break
        positive_counts.append(count)

    counts = [scale - sum(positive_counts), *positive_counts]
    assert counts[-1] > 0
    assert sum(counts) == scale

    # With thresholds T_i and uniform U in [0,2^112), the sampler returns the
    # number of T_i strictly smaller than U.  Hence T_i = cumulative_count-1.
    thresholds: List[int] = []
    cumulative = 0
    for count in counts[:-1]:
        cumulative += count
        thresholds.append(cumulative - 1)
    assert len(thresholds) == len(counts) - 1
    assert thresholds[-1] < scale - 1
    return thresholds, counts, probs


def cdt_renyi_log2(
    counts: Sequence[int], probs: Sequence[mp.mpf], alpha: int
) -> mp.mpf:
    """Compute log2(R_alpha(q || p)-1) at high precision."""
    scale = mp.mpf(2) ** CDT_BITS
    terms = []
    for count, ideal in zip(counts, probs):
        if count:
            q = mp.mpf(count) / scale
            terms.append(q**alpha * ideal ** (1 - alpha))
    r_alpha = mp.fsum(terms) ** (mp.mpf(1) / (alpha - 1))
    return mp.log(r_alpha - 1, 2)


def _poly_add(a: Sequence[mp.mpf], b: Sequence[mp.mpf]) -> List[mp.mpf]:
    n = max(len(a), len(b))
    return [
        (a[i] if i < len(a) else mp.mpf(0))
        + (b[i] if i < len(b) else mp.mpf(0))
        for i in range(n)
    ]


def _poly_scale(a: Sequence[mp.mpf], scale: mp.mpf) -> List[mp.mpf]:
    return [scale * x for x in a]


def _poly_mul_t(a: Sequence[mp.mpf], zmax: mp.mpf) -> List[mp.mpf]:
    """Multiply a power-basis polynomial by t = 2*z/zmax - 1."""
    result = [mp.mpf(0)] * (len(a) + 1)
    for i, value in enumerate(a):
        result[i] -= value
        result[i + 1] += (2 / zmax) * value
    return result


def generate_exp_coefficients(
    dps: int,
) -> Tuple[List[int], mp.mpf, mp.mpf, mp.mpf]:
    """Return coefficients, tail bound, relative bound, and ceiling margin."""
    with mp.workdps(dps):
        degree = EXP_DEGREE
        zmax = mp.mpf(EXP_ZMAX_NUM) / EXP_ZMAX_DEN
        center = zmax / 2

        # T_k(2*z/zmax-1), represented in the power basis in z.
        cheb = [[mp.mpf(1)], [mp.mpf(-1), 2 / zmax]]
        for _ in range(2, degree + 1):
            cheb.append(
                _poly_add(
                    _poly_scale(_poly_mul_t(cheb[-1], zmax), 2),
                    _poly_scale(cheb[-2], -1),
                )
            )

        power = [mp.mpf(0)] * (degree + 1)
        for k in range(degree + 1):
            if k == 0:
                cheb_coefficient = mp.exp(-center) * mp.besseli(0, center)
            else:
                cheb_coefficient = (
                    2
                    * mp.exp(-center)
                    * ((-1) ** k)
                    * mp.besseli(k, center)
                )
            for j, value in enumerate(cheb[k]):
                power[j] += cheb_coefficient * value

        # For k0=d+1, use
        # I_k(c) <= (c/2)^k/k! * exp(c^2/(4(k+1))).
        # Later leading terms have ratio <= c/(2(k+1)), giving a geometric
        # bound for the entire omitted Chebyshev tail.
        k0 = degree + 1
        bessel_tail_bound = (
            2
            * mp.exp(-center)
            * (center / 2) ** k0
            / mp.factorial(k0)
            * mp.exp(center * center / (4 * (k0 + 1)))
            / (1 - center / (2 * (k0 + 1)))
        )
        lift = mp.mpf(2) ** (-EXP_LIFT_BITS)
        if not bessel_tail_bound < lift:
            raise RuntimeError("the configured lift does not cover the Chebyshev tail")
        power[0] += lift

        scale = 1 << EXP_BITS
        coefficients = [int(mp.ceil(value * scale)) for value in power]

        # Bound coefficient-ceiling and Horner-ceiling errors.  Each error is
        # below one Q63 ulp and z is non-negative.
        coefficient_error = mp.fsum(zmax**k for k in range(degree + 1)) / scale
        horner_error = mp.fsum(zmax**k for k in range(degree)) / scale
        absolute_error_bound = (
            lift + bessel_tail_bound + coefficient_error + horner_error
        )
        input_error = EXP_INPUT_ERROR_ULPS / scale

        # If w is the implemented input and z the exact input, |w-z| is at
        # most input_error.  The lifted real polynomial remains above the
        # exact target provided lift-tail >= expm1(input_error); using 1 as an
        # upper bound for exp(-w) is conservative over the whole domain.
        ceiling_margin = lift - bessel_tail_bound - mp.expm1(input_error)
        if ceiling_margin < 0:
            raise RuntimeError("the lift does not cover the exp-input error")

        # For the other direction, the input perturbation contributes at most
        # expm1(input_error) to the relative error.  The polynomial/rounding
        # absolute error is maximized relative to exp(-z) at z=zmax.
        relative_error_bound = (
            absolute_error_bound / mp.exp(-zmax) + mp.expm1(input_error)
        )
        return (
            coefficients,
            bessel_tail_bound,
            relative_error_bound,
            ceiling_margin,
        )


def _ceil_div_pow2(value: int, bits: int) -> int:
    return (value + (1 << bits) - 1) // (1 << bits)


def eval_exp_q63(coefficients: Sequence[int], xi: int) -> int:
    """Bit-exact scalar Q63 Horner evaluation, including probability clamp."""
    result = coefficients[-1]
    for coefficient in reversed(coefficients[:-1]):
        result = _ceil_div_pow2(result * xi, EXP_BITS) + coefficient
    # The lifted constant is 2^63+1.  A rejection threshold cannot exceed 1.
    return min(result, 1 << EXP_BITS)


def dense_exp_scan(coefficients: Sequence[int], points: int) -> Tuple[mp.mpf, mp.mpf]:
    if points <= 0:
        return mp.nan, mp.nan
    scale = 1 << EXP_BITS
    xmax = (EXP_ZMAX_NUM * scale) // EXP_ZMAX_DEN
    maximum = mp.mpf(0)
    minimum = mp.inf
    for i in range(points + 1):
        xi = (xmax * i) // points
        implemented_z = mp.mpf(xi) / scale
        perturbation = EXP_INPUT_ERROR_ULPS / scale
        exact_z_lo = max(mp.mpf(0), implemented_z - perturbation)
        exact_z_hi = min(mp.mpf(EXP_ZMAX_NUM) / EXP_ZMAX_DEN, implemented_z + perturbation)
        probability = mp.mpf(eval_exp_q63(coefficients, xi)) / scale
        relative_lo = (probability - mp.exp(-exact_z_lo)) / mp.exp(-exact_z_lo)
        relative_hi = (probability - mp.exp(-exact_z_hi)) / mp.exp(-exact_z_hi)
        maximum = max(maximum, relative_hi)
        minimum = min(minimum, relative_lo)
    return minimum, maximum


def _format_c_array(name: str, c_type: str, values: Iterable[int], width: int) -> str:
    values = list(values)
    macro = "UINT64_C" if c_type == "uint64_t" else "UINT64_C"
    lines = [f"static const {c_type} {name}[CDTLEN] = {{"]
    row: List[str] = []
    for value in values:
        row.append(f"{macro}(0x{value:0{width}x})")
        if len(row) == 3:
            lines.append("    " + ", ".join(row) + ",")
            row = []
    if row:
        lines.append("    " + ", ".join(row) + ",")
    lines.append("};")
    return "\n".join(lines)


def emit_cdt(thresholds: Sequence[int]) -> None:
    high = [value >> 64 for value in thresholds]
    low = [value & ((1 << 64) - 1) for value in thresholds]
    print(f"#define CDTLEN {len(thresholds)}")
    print(_format_c_array("CDT112_HI", "uint64_t", high, 12))
    print()
    print(_format_c_array("CDT112_LO", "uint64_t", low, 16))


def emit_exp(coefficients: Sequence[int]) -> None:
    print("/* Includes a 0.500061-Q63-ulp bound for exp-input quantization. */")
    print(f"#define EXP_FRAC_BITS {EXP_BITS}")
    print(f"#define EXP_DEGREE {EXP_DEGREE}")
    print(f"#define EXP_ZMAX_NUM {EXP_ZMAX_NUM}")
    print(f"#define EXP_ZMAX_DEN {EXP_ZMAX_DEN}")
    print("/* Coefficients are c_0 through c_12; c_0 is intentionally unsigned. */")
    for i, value in enumerate(coefficients):
        if i == 0:
            literal = f"UINT64_C({value})"
        elif value < 0:
            literal = f"(-INT64_C({-value}))"
        else:
            literal = f"INT64_C({value})"
        print(f"#define EXP_C{i} {literal}")


def emit_header(thresholds: Sequence[int], coefficients: Sequence[int]) -> None:
    print("#ifndef DARTS_GAUSS_Q80_CONSTANTS_H")
    print("#define DARTS_GAUSS_Q80_CONSTANTS_H")
    print()
    print("#include <stdint.h>")
    print()
    print("/* Generated by helper_scripts/gen_gauss_q80.py --emit all. */")
    print("#define CDT_BITS 112")
    emit_cdt(thresholds)
    print()
    emit_exp(coefficients)
    print()
    print("#endif")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--emit", choices=("none", "cdt", "exp", "all"), default="none")
    parser.add_argument("--scan-points", type=int, default=200_000)
    parser.add_argument("--quiet", action="store_true", help="suppress the audit report on stderr")
    args = parser.parse_args()

    mp.mp.dps = 180
    thresholds_180, counts_180, probs_180 = generate_cdt(180)
    thresholds_300, counts_300, probs_300 = generate_cdt(300)
    if thresholds_180 != thresholds_300 or counts_180 != counts_300:
        raise RuntimeError("CDT generation is precision-dependent")
    max_base_sample = len(counts_300) - 1
    if EXP_ZMAX_NUM != 2 * max_base_sample + 1:
        raise RuntimeError(
            "exp domain is inconsistent with the CDT support: "
            f"{EXP_ZMAX_NUM}/512 != (2*{max_base_sample}+1)/512"
        )

    exp_180, tail_180, rel_180, margin_180 = generate_exp_coefficients(180)
    exp_300, tail_300, rel_300, margin_300 = generate_exp_coefficients(300)
    if exp_180 != exp_300:
        raise RuntimeError("exp coefficient generation is precision-dependent")
    exp_scale = 1 << EXP_BITS
    if exp_300[0] != exp_scale + 2:
        raise RuntimeError("unexpected lifted Q63 constant coefficient")
    if any(not (-(1 << 63) <= value < (1 << 63)) for value in exp_300[1:]):
        raise RuntimeError("a signed Q63 coefficient does not fit int64_t")
    exp_xmax = (EXP_ZMAX_NUM * exp_scale) // EXP_ZMAX_DEN
    if eval_exp_q63(exp_300, 0) != exp_scale:
        raise RuntimeError("clamped Q63 approximation is not exactly one at zero")
    if not (0 < eval_exp_q63(exp_300, exp_xmax) <= exp_scale):
        raise RuntimeError("Q63 approximation is outside the probability range")

    renyi_results = {
        alpha: cdt_renyi_log2(counts_300, probs_300, alpha)
        for alpha in RENYI_ORDERS
    }
    for alpha, value in renyi_results.items():
        if value > CDT_LOG2_RENYI_TARGET:
            raise RuntimeError(
                f"CDT112 fails the alpha={alpha} Renyi target: "
                f"{value} > {CDT_LOG2_RENYI_TARGET}"
            )
    exp_k_bound = -mp.log(rel_300, 2)
    if exp_k_bound < EXP_K_TARGET:
        raise RuntimeError(
            f"Q63 degree-12 approximation fails the K target: "
            f"{exp_k_bound} < {EXP_K_TARGET}"
        )

    if not args.quiet:
        print(
            f"CDT112: {len(thresholds_300)} thresholds, outputs 0..{len(counts_300)-1}",
            file=sys.stderr,
        )
        for alpha in RENYI_ORDERS:
            value = renyi_results[alpha]
            print(f"  alpha={alpha:4d}: log2(R_alpha-1)={mp.nstr(value, 18)}", file=sys.stderr)

        print(
            f"EXP Q63 degree 12: Chebyshev tail bound < 2^-{EXP_LIFT_BITS}",
            file=sys.stderr,
        )
        print(f"  computed tail upper bound: 2^{mp.nstr(mp.log(tail_300, 2), 18)}", file=sys.stderr)
        print(
            f"  analytic relative-error bound: K >= {mp.nstr(exp_k_bound, 18)}",
            file=sys.stderr,
        )
        print(
            "  exact-input ceiling margin: "
            f"{mp.nstr(margin_300 * (1 << EXP_BITS), 18)} Q63 ulp",
            file=sys.stderr,
        )
    scan_min, scan_max = dense_exp_scan(exp_300, args.scan_points)
    if args.scan_points > 0 and not args.quiet:
        print(
            f"  dense scan ({args.scan_points} intervals): min_rel={mp.nstr(scan_min, 8)}, "
            f"K_observed={mp.nstr(-mp.log(scan_max, 2), 18)}",
            file=sys.stderr,
        )
    if args.scan_points > 0:
        if scan_min < 0:
            raise RuntimeError("dense scan found an under-accepting input")

    if args.emit == "all":
        emit_header(thresholds_300, exp_300)
    elif args.emit == "cdt":
        emit_cdt(thresholds_300)
    elif args.emit == "exp":
        emit_exp(exp_300)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
