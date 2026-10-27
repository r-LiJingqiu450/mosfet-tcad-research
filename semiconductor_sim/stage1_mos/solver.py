# -*- coding: utf-8 -*-
"""stage1_mos.solver —— MOS 电容非线性泊松方程求解器（Newton-Raphson + 追赶法）。

控制方程（x 垂直于 SiO2/Si 界面，指向衬底；phi 相对衬底体内 phi_bulk = 0）：
    d^2 phi/dx^2 = -(q/eps_si) (p - n - N_a)
    p = p0 exp(-phi/V_T),  n = n0 exp(+phi/V_T)
边界：phi(0) = phi_s（Dirichlet），phi(L) = 0（Dirichlet）。

数值格式：非均匀网格二阶中心差分（均匀网格退化为标准形式），三对角矩阵
严格对角占优（M-矩阵）；收敛控制三层：Armijo 阻尼回溯线搜索、步长限幅
clamp、参数延续（调用方以上一步表面势的解为初值）。

依赖：common.constants（常数）、common.thomas（追赶法）。
本文件不含绘图与主流程（见 make_figures.py / run.py）。"""

import warnings

import numpy as np

from common.constants import EPS_SI, K_B, NI, Q_E, TEMP, VT
from common.thomas import thomas_solve

def _residual(phi, p0, n0, Na, aa, bb0, cc):
    """离散泊松方程残差 F(phi)（用于阻尼线搜索）。

    F[j] 对应内部点 i=j+1：
        F[j] = aa[j]*phi_j + bb0[j]*phi_{j+1} + cc[j]*phi_{j+2} + 源项
    aa 乘左邻居，cc 乘右邻居（均匀网格时两者相等，非均匀时不可交换）。
    """
    phi_c = np.clip(phi, -3.0, 3.0)
    p = p0 * np.exp(-phi_c / VT)
    n = n0 * np.exp(+phi_c / VT)
    lap = aa * phi[:-2] + bb0 * phi[1:-1] + cc * phi[2:]
    return lap + (Q_E / EPS_SI) * (p[1:-1] - n[1:-1] - Na)


def solve_poisson(phi_s, x, Na, phi_init=None,
                  tol=1e-11, max_iter=200, clamp=0.3, warn=True,
                  damped=True, max_backtrack=12):
    """
    给定表面势 phi_s，求解一维泊松方程得到电势分布 phi(x)。

    控制方程（x 垂直于界面，指向衬底）：
        d^2 phi/dx^2 = -(q/eps) (p - n - Na)
        p = p0 exp(-phi/VT),  n = n0 exp(+phi/VT)

    边界：phi(0)=phi_s (Dirichlet), phi(-1)=0 (Dirichlet)

    数值格式
    --------
    非均匀二阶中心差分（均匀网格时退化为标准形式）：
        phi''_i = 2/(hL+hR) [ (phi_{i+1}-phi_i)/hR - (phi_i-phi_{i-1})/hL ]
    => 三对角系数
        a_i = 2/(hL(hL+hR)),  b_i = -2/(hL·hR),  c_i = 2/(hR(hL+hR))
    均匀网格时 a_i = c_i = 1/h^2, b_i = -2/h^2，与标准差分完全一致。

    对角占优性（M-矩阵性质，保证 Newton 每步良态）：
        |b_i| = 2/(hL·hR) + (q^2/(eps·kT))(p+n)  >  a_i + c_i = 2/(hL·hR)
    主对角元严格大于非对角元之和 -> 严格对角占优，三对角矩阵非奇异。

    收敛控制（三层）
    ----------------
    1. 阻尼回溯线搜索（Armijo）：若全残差 2-范数不下降，步长减半；
       解决非均匀网格/大初值扰动下的指数项发散问题；
    2. 步长限幅 clamp：单次 |Delta phi| 截断，防 exp 溢出；
    3. 参数延续（continuation）：调用方以上一步表面势的解作为初值。

    参数
    ----
    phi_s    : 表面势 (V)，相对于衬底体内
    x        : 空间网格 (m)
    Na       : 受主掺杂浓度 (m^-3)
    phi_init : 初值猜测
    tol      : 收敛判据（相邻迭代 |Delta phi| 最大值）
    max_iter : 最大迭代次数，未收敛时发 RuntimeWarning 并返回
    clamp    : 单步限幅 (V)
    damped   : 是否启用阻尼线搜索

    返回
    ----
    phi, niter, converged(bool)
    """
    N = len(x)
    hs = np.diff(x)
    hL = hs[:-1]
    hR = hs[1:]
    sumLR = hL + hR

    aa = 2.0 / (hL * sumLR)
    bb0 = -2.0 / (hL * hR)
    cc = 2.0 / (hR * sumLR)

    p0 = 0.5 * (Na + np.sqrt(Na ** 2 + 4.0 * NI ** 2))
    n0 = NI ** 2 / p0

    if phi_init is None:
        phi = phi_s * (1.0 - x / x[-1])
    else:
        phi = np.interp(x, np.linspace(0, x[-1], len(phi_init)), phi_init)
    phi = phi.astype(float)
    phi[0] = phi_s
    phi[-1] = 0.0

    converged = False
    niter = max_iter
    F = _residual(phi, p0, n0, Na, aa, bb0, cc)
    normF = np.linalg.norm(F)

    for niter in range(1, max_iter + 1):
        phi_c = np.clip(phi, -3.0, 3.0)
        p = p0 * np.exp(-phi_c / VT)
        n = n0 * np.exp(+phi_c / VT)

        b_diag = bb0 - (Q_E ** 2 / (EPS_SI * K_B * TEMP)) * (p[1:-1] + n[1:-1])
        a_low = aa.copy(); a_low[0] = 0.0
        c_up = cc.copy(); c_up[-1] = 0.0

        delta = thomas_solve(a_low, b_diag, c_up, -F)
        delta = np.clip(delta, -clamp, clamp)

        # ---- 阻尼回溯线搜索：保证残差 2-范数单调下降 ----
        alpha = 1.0
        if damped:
            for _ in range(max(1, int(max_backtrack))):
                phi_try = phi.copy()
                phi_try[1:-1] += alpha * delta
                F_try = _residual(phi_try, p0, n0, Na, aa, bb0, cc)
                normF_try = np.linalg.norm(F_try)
                if normF_try <= (1.0 - 2e-4 * alpha) * normF:
                    break
                alpha *= 0.5
            # 残差地板豁免：离散残差是量级 ~ q*Na/eps_si*N^(1/2) 的大项相减，
            # 当其相对量级低于 1e-8 时已受 double 舍入主导，Armijo 条件
            # 无法被分辨，属数值地板而非线搜索失败
            src_scale = (Q_E / EPS_SI) * Na * np.sqrt(max(len(phi) - 2, 1))
            if normF_try > (1.0 - 2e-4 * alpha) * normF and normF > 1e-8 * src_scale:
                # 回溯 12 次仍未满足 Armijo 条件且残差远离地板：告警
                if warn:
                    warnings.warn(
                        f"solve_poisson: Armijo 线搜索折半 "
                        f"{max(1,int(max_backtrack))} 次仍未"
                        f"满足充分下降条件 (phi_s={phi_s:.4f} V, "
                        f"alpha={alpha:.2e})，该步接受最小步长解。",
                        RuntimeWarning, stacklevel=2)
            phi, F, normF = phi_try, F_try, normF_try
            dmove = alpha * np.max(np.abs(delta))
        else:
            phi[1:-1] += delta
            F = _residual(phi, p0, n0, Na, aa, bb0, cc)
            normF = np.linalg.norm(F)
            dmove = np.max(np.abs(delta))

        phi[0] = phi_s
        phi[-1] = 0.0

        if dmove < tol:
            converged = True
            break

    if not converged and warn:
        warnings.warn(
            f"solve_poisson: Newton 未在 {max_iter} 次迭代内收敛 "
            f"(phi_s={phi_s:.4f} V), 返回末次迭代解。",
            RuntimeWarning, stacklevel=2)

    return phi, niter, converged
