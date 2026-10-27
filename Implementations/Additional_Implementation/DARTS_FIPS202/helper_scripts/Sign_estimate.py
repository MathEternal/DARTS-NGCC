from math import isfinite, log

from security_estimate.LWE_estimate import *
from security_estimate.SIS_estimate import *
from helper_pycode.repetitions_compute import (
    from_B11_to_B0,
    from_B11_to_B1,
    from_B11_to_BSIS,
)
from helper_pycode.entropy_coordinate_refined_hyperball import *


# Implementation-layout data and empirical h entropy are intentionally tied to
# the complete, confirmed parameter fingerprint.  Reusing either value for a
# new candidate would silently produce a misleading size estimate.
#
# h_entropy is computed from COUNTS_H in the corrected 200000-signature rANS
# fixed-point reports (SHA-256 prefixes: 128=f142579a, 256=0449dfac,
# 512=5a7b6a5d).  Refresh it after any parameter or signing-path change.
SIZE_PROFILES = {
    (512, 1, 2, 130817, 0.15, 30, 10, 7500.00, 191.26): {
        "seed_bytes": 32,
        "lowbits_packed_bytes": 456,
        "h_entropy": 2.408300028561897,
    },
    (512, 2, 3, 130817, 0.3125, 59, 9, 23000.00, 476.38): {
        "seed_bytes": 32,
        "lowbits_packed_bytes": 520,
        "h_entropy": 2.6713478956882128,
    },
    (1024, 2, 3, 260609, 0.20, 115, 10, 80000.00, 780.72): {
        "seed_bytes": 64,
        "lowbits_packed_bytes": 1024,
        "h_entropy": 4.11038329585329,
    },
}


def _size_profile_key(n, k, l, q, p, tau, d, B11, t):
    return (
        n,
        k,
        l,
        q,
        round(float(p), 12),
        tau,
        d,
        round(float(B11), 2),
        round(float(t), 2),
    )


def _positive_integer_or_none(value, name):
    if value is None:
        return None
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a positive integer")
    try:
        converted = int(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a positive integer") from exc
    if converted <= 0 or value != converted:
        raise ValueError(f"{name} must be a positive integer")
    return converted


class SignParameters:
    def __init__(
        self,
        n,
        k,
        l,
        q,
        p,
        tau,
        d,
        B11,
        t,
        seed_bytes=None,
        *,
        lowbits_packed_bytes=None,
        h_entropy=None,
    ):
        self.n = n      # 多项式维数
        self.k = k      # 行数
        self.l = l      # 列数
        self.q = q      # 模数
        self.p = p      # 三元分布取 1 (-1) 的概率
        self.tau = tau  # 挑战 c 中非零系数的数量
        self.d = d      # 压缩位数
        self.B11 = B11  # 重复采样参数 B''
        self.t = t      # sc 的上界 t = gamma * sqrt(tau)

        profile = SIZE_PROFILES.get(
            _size_profile_key(n, k, l, q, p, tau, d, B11, t)
        )
        if profile is not None:
            if seed_bytes is None:
                seed_bytes = profile["seed_bytes"]
            if lowbits_packed_bytes is None:
                lowbits_packed_bytes = profile["lowbits_packed_bytes"]
            if h_entropy is None:
                h_entropy = profile["h_entropy"]

        seed_bytes = _positive_integer_or_none(seed_bytes, "seed_bytes")
        lowbits_packed_bytes = _positive_integer_or_none(
            lowbits_packed_bytes, "lowbits_packed_bytes"
        )
        if h_entropy is not None:
            try:
                h_entropy = float(h_entropy)
            except (TypeError, ValueError, OverflowError) as exc:
                raise ValueError(
                    "h_entropy must be a finite positive number"
                ) from exc
            if not isfinite(h_entropy) or h_entropy <= 0:
                raise ValueError(
                    "h_entropy must be a finite positive number"
                )

        self.seed_bytes = seed_bytes
        self.lowbits_packed_bytes = lowbits_packed_bytes
        self.h_entropy = h_entropy


def Sign_to_MLWE(ps):
    """
    将签名参数转换为 MLWE 参数集。

    输入:
        ps: 签名参数
    输出:
        n: 环维数
        d: MLWE 的维数
        m: MLWE 的样本数
    """
    d = ps.l             # MLWE 维数
    m = ps.k             # MLWE 样本数
    return MLWEParameterSet(ps.n, d, m, ps.q, ps.p)


def Sign_to_MSIS(ps):
    """
    将签名参数转换为 MSIS 参数集。

    输入:
        ps: 签名参数
    输出:
        n: 环维数
        h: MSIS 的行数，即方程数
        w: MSIS 的列数，即变量数
        B: MSIS 的范数上界（采用2范数估计安全性）
    """
    B = from_B11_to_BSIS(ps.B11, ps.n, ps.k, ps.l, ps.q, ps.d)
    h = ps.k
    w = ps.k + ps.l
    return MSISParameterSet(ps.n, w, h, B, ps.q, method="refined_l2")


def compute_sign_size(ps):
    """
    估计 rANS 编码后的平均签名大小（字节）。

    h 使用修正统计程序得到的独立经验熵；解析超球模型只用于
    HighBits(z1)，不能把它复用于 h。低位尺寸采用实现的 BAT 打包
    字节数，并计入两路 rANS 各自 4 字节的状态刷新开销。

    该值仍是理想熵估计，不包含 1024 频率量化的交叉熵损失和逐字节
    取整，也不是由 CRYPTO_SIGNATUREBYTES 定义的固定 wire-format 长度。

    输入:
        ps: 签名参数
    输出:
        sign_size_byte: 签名压缩尺寸的熵估计值
    """
    if ps.lowbits_packed_bytes is None:
        raise ValueError(
            "No BAT lowbits layout is known for this parameter set; pass "
            "the per-polynomial lowbits_packed_bytes explicitly."
        )
    if ps.h_entropy is None:
        raise ValueError(
            "No empirical h entropy is known for this parameter set; run the "
            "rANS statistics driver and pass h_entropy explicitly."
        )

    m = ps.k + ps.l

    B1 = from_B11_to_B1(ps.B11, ps.n, ps.k, ps.l, ps.q, ps.d)
    B0 = from_B11_to_B0(ps.B11, ps.n, ps.k, ps.l, ps.q, ps.d, ps.t)

    # 解析模型描述的是 HighBits(z1)，不是 h。
    (hb_prob, _) = refined_compute_f_table(ps.n, m, B1, B0, ps.d, ps.q)
    hb_entropy = -sum(p * log(p, 2) for p in hb_prob if p > 0)

    challenge_bytes = (ps.n + 7) // 8
    lowbits_bytes = ps.l * ps.lowbits_packed_bytes
    length_bytes = 2
    rans_state_bytes = 2 * 4
    h_bytes = ps.h_entropy * ps.n * ps.k / 8
    hb_z1_bytes = hb_entropy * ps.n * ps.l / 8

    return (
        challenge_bytes
        + lowbits_bytes
        + length_bytes
        + rans_state_bytes
        + h_bytes
        + hb_z1_bytes
    )


def compute_pk_size(ps):
    """
    计算公钥大小（字节）。

    公钥包括 seed 和矩阵 A_0。每个多项式按实现单独向上取整打包，
    系数位宽为 ceil(log2(q)) = (q-1).bit_length()。

    输入:
        ps: 签名参数
    输出:
        公钥大小（字节）
    """
    if ps.seed_bytes is None:
        raise ValueError(
            "No seed layout is known for this parameter set; pass seed_bytes "
            "explicitly."
        )
    # int() keeps this compatible with Sage Integer objects used by the
    # notebook kernel (which do not consistently expose int.bit_length()).
    q_bits = (int(ps.q) - 1).bit_length()
    poly_q_packed_bytes = (int(ps.n) * q_bits + 7) // 8
    return int(ps.seed_bytes) + int(ps.k) * poly_q_packed_bytes


def summarize_sign(ps):
    pk_size = compute_pk_size(ps)
    estimated_sign_size = compute_sign_size(ps)

    print("Sign parameters:", ps.__dict__)
    print(
        "Bound for SIS:",
        from_B11_to_BSIS(ps.B11, ps.n, ps.k, ps.l, ps.q, ps.d),
    )
    print("Public key size (bytes): ", pk_size)
    print(
        "Estimated mean compressed signature size (bytes): ",
        f"{estimated_sign_size:.2f}",
    )
    print(
        "Estimated public key + mean signature size (bytes):",
        f"{pk_size + estimated_sign_size:.2f}",
    )
    print("MLWE security:")
    MLWE_summarize_attacks(Sign_to_MLWE(ps))
    print("MSIS security:")
    MSIS_summarize_attacks(Sign_to_MSIS(ps))
