from math import sqrt


# 估计离散化系数 N 的值（下界）

# 默认 M0 和 c 的值
M0 = 1.001
c = 1.001


# 半径参数与 security_and_size_estimate.ipynb 中确认的候选值保持一致。
# m 不单独存储，计算时始终使用 k + l，避免参数集更新时不同步。
PARAMETER_SETS = (
    {
        "security": 128,
        "n": 512,
        "k": 1,
        "l": 2,
        "B": 6057.90,
        "B0": 5866.64,
        "B1": 6054.88,
    },
    {
        "security": 256,
        "n": 512,
        "k": 2,
        "l": 3,
        "B": 18916.70,
        "B0": 18440.32,
        "B1": 18910.70,
    },
    {
        "security": 512,
        "n": 1024,
        "k": 2,
        "l": 3,
        "B": 74220.97,
        "B0": 73440.24,
        "B1": 74216.86,
    },
)


def N_1_compute(m, n, B, M0=None):
    """
    利用连续的超球采样实现离散超球采样时对 N 的要求。

    :param m: m = k + l
    :param n: 环的维数
    :param B: 采样 y 时的超球半径
    :param M0: 成功概率为 1/M0
    """
    if M0 is None:
        M0 = 1.001
    return (sqrt(m * n) / (2 * B)) * (
        (M0 ** (1 / (m * n)) + 1) / (M0 ** (1 / (m * n)) - 1)
    )


def N_2_compute(m, n, B, B_candidate, c=None):
    """
    双峰拒绝采样时对 N 的要求。

    :param m: m = k + l
    :param n: 环的维数
    :param B: 采样 y 时的超球半径
    :param B_candidate: 候选球半径（B0 或 B1）
    :param c: M' = cM，其中 M' 为离散超球的拒绝采样次数，
              M 为连续超球的拒绝采样次数
    """
    if c is None:
        c = 1.001
    return (
        1 / (c ** (1 / (m * n)) - 1)
        * sqrt(m * n / 4)
        * (c ** (1 / (m * n)) / B + 1 / B_candidate)
    )


def main():
    for params in PARAMETER_SETS:
        m = params["k"] + params["l"]
        n = params["n"]
        bound1 = N_1_compute(m, n, params["B"], M0)
        bound2_b0 = N_2_compute(m, n, params["B"], params["B0"], c)
        bound2_b1 = N_2_compute(m, n, params["B"], params["B1"], c)

        print(f"Security {params['security']} (m = k + l = {m}):")
        print("N lower bound (continuous hyperball):", bound1)
        print("N lower bound (B0 candidate):", bound2_b0)
        print("N lower bound (B1 candidate):", bound2_b1)


if __name__ == "__main__":
    main()
