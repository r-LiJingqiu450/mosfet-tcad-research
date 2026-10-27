# -*- coding: utf-8 -*-
"""
系列统一三对角求解器。

三个入口：
    thomas_solve        标量三对角（纯 Python 追赶法，S1/S2/S3 原始实现）
    thomas_solve_banded 同签名，scipy LAPACK 带状求解优先（快 ~18x），
                        异常时自动回退纯 Python
    block_thomas        块三对角（块大小 3x3，S2 全 Newton 系统）

统一约定：a 下对角（a[0] 置 0）、b 主对角、c 上对角（c[-1] 置 0）、
d 右端项；返回解向量。
"""

import numpy as np

try:                                   # 可选加速：LAPACK 带状求解
    import scipy.linalg as _sla
    _HAVE_SCIPY = True
except ImportError:                    # pragma: no cover - 无 scipy 时回退
    _HAVE_SCIPY = False


def thomas_solve(a, b, c, d):
    """
    追赶法求解三对角方程组 A x = d。

    矩阵结构（下标 0..n-1）：
        b[0] c[0]
        a[1] b[1] c[1]
             ...
                 a[n-1] b[n-1]

    参数：a 下对角(a[0] 置 0)，b 主对角，c 上对角(c[-1] 置 0)，d 右端项
    返回：解向量 x
    """
    n = len(d)
    cp = np.zeros(n)   # 消元后的上对角
    dp = np.zeros(n)   # 消元后的右端项

    # ---- 前向消元 ----
    cp[0] = c[0] / b[0]
    dp[0] = d[0] / b[0]
    for i in range(1, n):
        denom = b[i] - a[i] * cp[i - 1]
        cp[i] = c[i] / denom
        dp[i] = (d[i] - a[i] * dp[i - 1]) / denom

    # ---- 回代 ----
    x = np.zeros(n)
    x[-1] = dp[-1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]

    return x


def thomas_solve_banded(a, b, c, d):
    """三对角求解：scipy LAPACK 带状优先，纯 Python 回退（与
    thomas_solve 同签名、同结果）。"""
    if _HAVE_SCIPY:
        n = len(d)
        ab = np.zeros((3, n))
        ab[0, 1:] = c[:-1]
        ab[1, :] = b
        ab[2, :-1] = a[1:]
        try:
            return _sla.solve_banded((1, 1), ab, d)
        except Exception:               # pragma: no cover - 奇异等异常时回退
            pass
    return thomas_solve(a, b, c, d)


def block_thomas(A_diag, A_low, A_up, rhs):
    """
    块三对角求解（块大小 3x3，块追赶法）。

    推导（块 LU / Schur 消元，与标量 Thomas 逐行对应）：
        行消元：row_i -= F_i @ row_{i-1}^{elim},  F_i = A_low[i] @ S_{i-1}^{-1}
        Schur 块：S_i = A_diag[i] - F_i @ A_up[i-1]
                  = A_diag[i] - A_low[i] @ S_{i-1}^{-1} @ A_up[i-1]
        右端：    g~_i = rhs[i] - F_i @ g~_{i-1}
        回代：    u[i] = S_i^{-1} (g~_i - A_up[i] @ u[i+1])
    代码中用 W[i] = S_i^{-1}、g[i] = S_i^{-1} @ g~_i 折叠表示。
    注意 A_up[i-1] 因子不可省（标量对应 cp[i] = c[i]/(b[i]-a[i]cp[i-1])
    中的 cp[i-1] 而非 1/cp[i-1]）；曾遗漏该因子导致 Newton 方向
    完全失真（DEF-S2-02，回归测试见 stage2_pn/tests/）。

    参数
    ----
    A_diag : (M,3,3) 主对角块
    A_low  : (M,3,3) 次对角块，A_low[i] 耦合节点 i 与 i-1（A_low[0] 须为 0）
    A_up   : (M,3,3) 超对角块，A_up[i] 耦合节点 i 与 i+1（A_up[-1] 须为 0）
    rhs    : (M,3)   右端项

    返回 (M,3) 解。
    """
    M = len(rhs)
    W = [None] * M
    g = [None] * M
    W[0] = np.linalg.inv(A_diag[0])
    g[0] = W[0] @ rhs[0]
    for i in range(1, M):
        Mi = A_diag[i] - A_low[i] @ W[i - 1] @ A_up[i - 1]
        W[i] = np.linalg.inv(Mi)
        g[i] = W[i] @ (rhs[i] - A_low[i] @ g[i - 1])
    u = [None] * M
    u[M - 1] = g[M - 1]
    for i in range(M - 2, -1, -1):
        u[i] = g[i] - W[i] @ (A_up[i] @ u[i + 1])
    return np.array(u)
