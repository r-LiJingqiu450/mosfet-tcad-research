# -*- coding: utf-8 -*-
"""stage2_pn.postprocess —— 电流密度 / 耗尽层宽度 / 能带图数据。

    current_density   SG 格式截面电流（默认取结中心截面，数值最稳健）
    depletion_width   电荷法（与解析定义严格对应）与 rho 衰减判据（对照）
    band_diagram      平衡/偏压能带数据（相对 p 侧深部，eV）

依赖：common.constants / common.bernoulli。"""

import numpy as np

from common.constants import D_N, D_P, EG, NI, Q_E, VT
from common.bernoulli import bernoulli


# ================================================================================
# 9. 后处理：电流 / 耗尽宽度 / 能带
# ================================================================================

def current_density(phi, n, p, x, idx_interface=None):
    """
    在指定截面计算总电流密度 J = J_n + J_p (A/m^2)。

    默认取结中心截面（x 最接近 0 的界面）：
    该截面载流子浓度低，SG 两项无大数相消，数值最稳健。

    同时返回各截面电流列表用于验证电流守恒（J_n + J_p 应处处相等）。
    """
    hw = np.diff(x)
    delta = (phi[:-1] - phi[1:]) / VT
    Bp = bernoulli(delta)
    Bm = bernoulli(-delta)

    # 电子电流（界面 i+1/2）
    Jn = (Q_E * D_N / hw) * (Bp * n[1:] - Bm * n[:-1])
    # 空穴电流（界面 i+1/2）
    Jp = -(Q_E * D_P / hw) * (Bm * p[1:] - Bp * p[:-1])
    Jt = Jn + Jp

    if idx_interface is None:
        idx_interface = int(np.argmin(np.abs(x[:-1] + 0.5 * hw)))
    J = float(Jt[idx_interface])
    return J, Jt, Jn, Jp


def depletion_width(phi, n, p, Na, Nd, x):
    """
    从数值解提取耗尽层宽度 W (m)，两种定义：

    定义 1（电荷法，与解析定义严格对应）：
        W = |Q_p|/(q N_a) + |Q_n|/(q N_d)
        其中 Q_p = ∫_{p区} rho dx（负），Q_n = ∫_{n区} rho dx（正）。
        积分域取整个半导体（远场准中性，电荷 ~0）。

    定义 2（空间电荷衰减判据）：
        p 侧 rho 从 -qN_a 衰减到 -qN_a/e 的位置 x_p，
        n 侧 rho 从 +qN_d 衰减到 +qN_d/e 的位置 x_n，
        W = x_n - x_p。该判据有任意性，作对照。

    返回 (W_charge, W_criterion)
    """
    rho_q = p - n + np.where(x >= 0.0, Nd, 0.0) - np.where(x < 0.0, Na, 0.0)
    hw = np.diff(x)

    # ---- 定义 1：电荷法 ----
    # p 侧总电荷（界面通量近似：rho 在控制体积中心）
    # 注意 rho_q 为数密度 (m^-3)，Q_p = ∫rho_q dx 单位为 m^-2；
    # W = |Q_p|/N_a + |Q_n|/N_d（数密度积分 / 掺杂浓度 = 长度）。
    Q_p = np.sum(0.5 * (rho_q[:-1] + rho_q[1:]) * hw * (x[1:] < 0.0))
    Q_n = np.sum(0.5 * (rho_q[:-1] + rho_q[1:]) * hw * (x[:-1] >= 0.0))
    W_charge = abs(Q_p) / Na + abs(Q_n) / Nd

    # ---- 定义 2：rho 衰减判据 ----
    x_p, x_n = x[0], x[-1]
    idx_j = int(np.argmin(np.abs(x)))
    # p 侧
    seg = rho_q[:idx_j + 1]
    target = -Na / np.e
    idx = np.where(seg >= target)[0]
    if len(idx) > 0:
        x_p = x[idx[0]]
    # n 侧
    seg = rho_q[idx_j:]
    target = Nd / np.e
    idx = np.where(seg <= target)[0]
    if len(idx) > 0:
        x_n = x[idx_j + idx[0]]
    W_crit = x_n - x_p

    return float(W_charge), float(W_crit)


def band_diagram(phi, V_a, Na, Nd):
    """
    平衡能带图数据（相对 p 侧深部，单位 eV）。

    返回 E_c, E_v, E_i, E_F（平衡态单费米能级）或
    偏压态 E_Fn, E_Fp（准费米势，相对同一参考）。
    """
    Ec = -phi                      # 电子伏特（q 已约去）
    Ev = Ec - EG
    Ei = Ec - 0.56
    return Ec, Ev, Ei
