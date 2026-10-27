# -*- coding: utf-8 -*-
"""stage2_pn.solver —— 一维对称突变 PN 结直流求解器（全 Newton，准费米势变量）。

物理模型（300 K，Si，N_a = N_d = 1e16 cm^-3）：
    泊松    d/dx(eps dphi/dx) = -q (p - n + N_d - N_a)
    连续性  dJ_n/dx = qU,  dJ_p/dx = -qU
    统计    n = n_i exp((phi-phi_n)/V_T),  p = n_i exp((phi_p-phi)/V_T)
    SRH     U = (np - n_i^2) / [tau_p(n+n_i) + tau_n(p+n_i)]

数值方法：
    电流离散 Scharfetter-Gummel 指数插值（对任意网格间距二阶精度）；
    未知量取 (phi, phi_n, phi_p) 准费米势，块矩阵条件数从 ~1e36 降到 O(1)
    （列变换见 common.qfp）；块 Thomas 追赶法解块三对角系统（common.thomas）；
    Armijo 回溯 + 参数延续。Gummel 解耦迭代在 SRH 复合下实测停滞，
    仅保留作对照（gummel_solve）。

依赖：common.constants / common.thomas / common.bernoulli / common.qfp。
后处理见 postprocess.py，主流程见 run.py。"""

import warnings

import numpy as np

from common.constants import (D_N, D_P, EPS_SI, NI, Q_E, TAU_N, TAU_P, VT,
                              V_bi_analytic, W_analytic)
from common.bernoulli import bernoulli, bernoulli_deriv
from common.thomas import thomas_solve, block_thomas
from common.qfp import qfp_columns


def doping_profile(x, Na, Nd):
    """突变结掺杂分布 (m^-3)：x<0 为 p 区 (N_a)，x>0 为 n 区 (N_d)。"""
    return np.where(x < 0.0, Na, Nd)


# ================================================================================
# 4. 载流子统计与 SRH 复合
# ================================================================================

def carrier_from_phi(phi, phi_n, phi_p):
    """
    由电势 phi 与准费米势 phi_n / phi_p 求载流子浓度（玻尔兹曼统计）：
        n = n_i exp((phi - phi_n)/V_T),  p = n_i exp((phi_p - phi)/V_T)
    phi_n, phi_p 可为标量（接触钉扎值）或数组（全空间分布）。
    """
    n = NI * np.exp((np.clip(phi, -50.0, 50.0) - phi_n) / VT)
    p = NI * np.exp((phi_p - np.clip(phi, -50.0, 50.0)) / VT)
    return n, p


def srh_rate(n, p, tau_n=TAU_N, tau_p=TAU_P, ni=NI):
    """
    SRH 复合率 U (m^-3 s^-1) 及其对 n、p 的偏导（中间能级陷阱）。

        U = (np - n_i^2) / [tau_p (n + n_i) + tau_n (p + n_i)]
        dU/dn = p/D - (np - n_i^2) tau_p / D^2
        dU/dp = n/D - (np - n_i^2) tau_n / D^2

    可证明 dU/dn、dU/dp 恒 > 0（分子展开后负项相消），
    该性质保证连续性方程线性化后的三对角矩阵主对角占优（见 6 节）。
    """
    den = tau_p * (n + ni) + tau_n * (p + ni)
    np_excess = n * p - ni * ni
    U = np_excess / den
    dU_dn = p / den - np_excess * tau_p / den ** 2
    dU_dp = n / den - np_excess * tau_n / den ** 2
    return U, dU_dn, dU_dp


# ================================================================================
# 5. 边界条件
# ================================================================================

def contact_values(V_a, Na, Nd):
    """
    理想欧姆接触边界值。

    返回
    ----
    phi_left  : p 侧接触电势 (V) = -V_T ln(N_a/n_i) + V_a
    phi_right : n 侧接触电势 (V) = +V_T ln(N_d/n_i)
    phi_n_L, phi_p_L : p 侧接触准费米势 = V_a
    phi_n_R, phi_p_R : n 侧接触准费米势 = 0

    自洽性：接触处载流子
        p 侧: p(0) = N_a, n(0) = n_i^2/N_a（准中性，与 V_a 无关）
        n 侧: n(L) = N_d, p(L) = n_i^2/N_d（准中性，与 V_a 无关）
    跨结压降 = phi_right - phi_left = V_bi - V_a。
    """
    phi_left = -VT * np.log(Na / NI) + V_a
    phi_right = VT * np.log(Nd / NI)
    return dict(phi_left=phi_left, phi_right=phi_right,
                phi_n_L=V_a, phi_p_L=V_a, phi_n_R=0.0, phi_p_R=0.0)


def contact_carriers(bc, Na, Nd):
    """由接触边界求接触节点载流子浓度（m^-3），返回 (n_L, p_L, n_R, p_R)。"""
    n_L = NI * np.exp((bc["phi_left"] - bc["phi_n_L"]) / VT)   # = n_i^2/N_a
    p_L = NI * np.exp((bc["phi_p_L"] - bc["phi_left"]) / VT)   # = N_a
    n_R = NI * np.exp((bc["phi_right"] - bc["phi_n_R"]) / VT)  # = N_d
    p_R = NI * np.exp((bc["phi_p_R"] - bc["phi_right"]) / VT)  # = n_i^2/N_d
    return n_L, p_L, n_R, p_R


# ================================================================================
# 6. 连续性方程：SG 离散的三对角求解
# ================================================================================

def _continuity_coeff(phi, hw, D, which):
    """
    组装连续性方程的三对角系数（SG 格式）。

    电子方程 dJ_n/dx = qU：
        J_{i+1/2} = (qD/h) [ B(delta) n_{i+1} - B(-delta) n_i ]
    空穴方程 dJ_p/dx = -qU：
        J_{i+1/2} = -(qD/h) [ B(-delta) p_{i+1} - B(delta) p_i ]
    其中 delta = (phi_i - phi_{i+1})/V_T，h 为界面间距。

    参数
    ----
    phi : 电势 (V)
    hw  : 界面间距数组（长度 N-1），hw[i] = x[i+1]-x[i]
    D   : 扩散系数 (m^2/s)
    which : 'n' 或 'p'

    返回 (a_lo, b_diag, c_up, ...) 的公共部分（线性化复合项由调用方附加）。
    """
    delta = (phi[:-1] - phi[1:]) / VT          # 长度 N-1，界面通量
    Bp = bernoulli(delta)
    Bm = bernoulli(-delta)
    qD_over_h = Q_E * D / hw                   # 长度 N-1
    if which == "n":
        # J_{i+1/2} = qD/h [B(d) n_{i+1} - B(-d) n_i]
        return qD_over_h * Bp, qD_over_h * Bm, delta, Bp, Bm
    else:
        # J_{i+1/2} = -(qD/h) [B(-d) p_{i+1} - B(d) p_i]
        return -qD_over_h * Bm, -qD_over_h * Bp, delta, Bp, Bm


def solve_electron_continuity(phi, p_curr, n_curr, x, bc, cc,
                              max_iter=50, tol=1e-10):
    """
    固定 phi 与 p，解电子连续性方程 dJ_n/dx = qU（三对角线性系统）。

    SRH 复合项线性化（半隐式）：
        U_i ~ U_i^* + (dU/dn)_i^* (n_i - n_i^*)
    线性化系数取当前迭代值（p_curr, n_curr），该处理保持三对角结构，
    且 dU/dn > 0 保证主对角占优（-A 为 M-矩阵）。

    边界：n(0) = n_i^2/N_a, n(L) = N_d（接触准中性值，见 contact_values）。

    返回新电子浓度分布 n_new (m^-3)。
    """
    N = len(x)
    hw = np.diff(x)                    # 界面间距，长度 N-1
    h_i = 0.5 * (hw[:-1] + hw[1:])     # 内部控制体积宽度，长度 N-2

    n_L, _, n_R, _ = contact_carriers(bc, cc["Na"], cc["Nd"])

    cBp, cBm, _, _, _ = _continuity_coeff(phi, hw, D_N, "n")
    # cBp[i] = (qD/h) B(delta_i)  乘 n_{i+1}
    # cBm[i] = (qD/h) B(-delta_i) 乘 n_i

    U_star, dU_dn, _ = srh_rate(n_curr, p_curr)

    # 三对角组装（线性系统：系数只依赖 phi、p_curr、n_curr，一次解出）
    # 内部节点 j = 0..N-3 对应原节点 i = j+1：
    #   右界面 (i+1/2) -> hw 索引 i -> cBp[i] 乘 n_{i+1}
    #   左界面 (i-1/2) -> hw 索引 i-1 -> cBm[i-1] 乘 n_{i-1}
    a_lo = cBm[:-1].copy()                     # 长度 N-2
    c_up = cBp[1:].copy()                      # 长度 N-2
    b_di = -(cBm[1:] + cBp[:-1]) - Q_E * h_i * dU_dn[1:-1]
    rhs = Q_E * h_i * (U_star[1:-1] - dU_dn[1:-1] * n_curr[1:-1])

    # 边界项移到右端
    rhs[0] -= a_lo[0] * n_L                    # j=0 (i=1) 左邻居 n_0
    rhs[-1] -= c_up[-1] * n_R                  # j=N-3 (i=N-2) 右邻居 n_{N-1}
    a_lo[0] = 0.0
    c_up[-1] = 0.0

    n_new = np.empty(N)
    n_new[1:-1] = np.clip(thomas_solve(a_lo, b_di, c_up, rhs), 0.0, None)
    n_new[0] = n_L
    n_new[-1] = n_R
    return n_new


def solve_hole_continuity(phi, n_curr, p_curr, x, bc, cc,
                          max_iter=50, tol=1e-10):
    """
    固定 phi 与 n，解空穴连续性方程 dJ_p/dx = -qU（三对角线性系统）。

    与电子子步骤完全对称：
        J_{i+1/2} = -(qD_p/h) [ B(-delta) p_{i+1} - B(delta) p_i ]
    边界：p(0) = N_a, p(L) = n_i^2/N_d。

    返回新空穴浓度分布 p_new (m^-3)。
    """
    N = len(x)
    hw = np.diff(x)
    h_i = 0.5 * (hw[:-1] + hw[1:])

    _, p_L, _, p_R = contact_carriers(bc, cc["Na"], cc["Nd"])

    cBp, cBm, _, _, _ = _continuity_coeff(phi, hw, D_P, "p")
    # 空穴: J_{i+1/2} = -(qD/h)[B(-d) p_{i+1} - B(d) p_i]
    # cBp[i] = -(qD/h) B(-delta_i)  乘 p_{i+1}
    # cBm[i] = -(qD/h) B(delta_i)   乘 p_i

    U_star, _, dU_dp = srh_rate(n_curr, p_curr)

    # 三对角组装（与电子子步骤对称，一次解出）
    # 注意：空穴方程右边为 -qU（电子为 +qU），故 SRH 线性化项
    # 移入主对角线时符号为 + q h_i (dU/dp)，与电子版的 - q h_i (dU/dn)
    # 相反。该符号错误会被病态矩阵（行和~0）放大，必须严格对应。
    a_lo = cBm[:-1].copy()
    c_up = cBp[1:].copy()
    b_di = -(cBm[1:] + cBp[:-1]) + Q_E * h_i * dU_dp[1:-1]
    rhs = -Q_E * h_i * (U_star[1:-1] - dU_dp[1:-1] * p_curr[1:-1])

    rhs[0] -= a_lo[0] * p_L
    rhs[-1] -= c_up[-1] * p_R
    a_lo[0] = 0.0
    c_up[-1] = 0.0

    p_new = np.empty(N)
    p_new[1:-1] = np.clip(thomas_solve(a_lo, b_di, c_up, rhs), 0.0, None)
    p_new[0] = p_L
    p_new[-1] = p_R
    return p_new


# ================================================================================
# 7. 泊松方程：阻尼 Newton 迭代（复用第一阶段方案）
# ================================================================================

def _poisson_residual(phi, n, p, Na_arr, Nd_arr, aa, bb0, cc, h_i):
    """
    离散泊松方程残差 F(phi)（有限体积形式，内部节点）。

        F_i = eps [(phi_{i+1}-phi_i)/hR - (phi_i-phi_{i-1})/hL]
              - q h_i (p_i - n_i + N_d,i - N_a,i)

    即"穿出控制体积的电位移通量 - 体积电荷"，单位 C/m^2。
    aa 乘左邻居，cc 乘右邻居（非均匀网格时不可交换）。
    """
    lap = aa * phi[:-2] + bb0 * phi[1:-1] + cc * phi[2:]
    return lap - Q_E * h_i * (
        p[1:-1] - n[1:-1] + Nd_arr[1:-1] - Na_arr[1:-1])


def solve_poisson_pn(phi, n, p, Na_arr, Nd_arr, x, bc,
                     tol=1e-11, max_iter=100, clamp=0.3,
                     damped=True, max_backtrack=12, warn=True):
    """
    阻尼 Newton 迭代求解泊松方程（n、p 固定为 Gummel 外层当前值）。

    离散格式（非均匀网格二阶中心差分，均匀网格退化为标准形式）：
        phi''_i = 2/(hL+hR) [ (phi_{i+1}-phi_i)/hR - (phi_i-phi_{i-1})/hL ]
    乘以 eps 与控制体积宽度 h_i = (hL+hR)/2 后并入源项。

    雅可比三对角（对内部节点，Gummel 解耦下 n、p 固定）：
        a_lo = -eps/hL
        b_di =  eps/hL + eps/hR
        c_up = -eps/hR

    矩阵性质：|b_di| = |a_lo| + |c_up|，对角占优取等；该矩阵为
    不可约弱对角占优的 M-矩阵（非奇异），且 n、p 固定时残差对 phi
    严格线性 —— Newton 一步收敛至机器精度（实测 2 步内）。

    收敛控制（三层，与第一阶段 MOS 仿真相同）：
        1) Armijo 阻尼回溯线搜索（残差 2-范数充分下降）；
        2) 步长限幅 clamp（单次 |Delta phi| 截断，防 exp 溢出）；
        3) 参数延续（调用方以上一偏压解为初值）。

    返回 (phi_new, niter, converged)
    """
    N = len(x)
    hw = np.diff(x)
    hL = hw[:-1]
    hR = hw[1:]
    h_i = 0.5 * (hL + hR)

    aa = -EPS_SI / hL
    bb0 = EPS_SI / hL + EPS_SI / hR
    cc = -EPS_SI / hR

    phi = phi.astype(float).copy()
    phi[0] = bc["phi_left"]
    phi[-1] = bc["phi_right"]

    F = _poisson_residual(phi, n, p, Na_arr, Nd_arr, aa, bb0, cc, h_i)
    normF = np.linalg.norm(F)
    converged = False
    niter = max_iter

    for niter in range(1, max_iter + 1):
        # ---- Jacobian（内部节点）----
        # 注意：Gummel 解耦下 n、p 为已知量（由连续性方程单独求解），
        # 泊松方程对 phi 是线性的，源项不含 phi 依赖：
        #     dF_i/dphi_i = eps(1/hL + 1/hR)
        # （这与第一阶段 MOS 仿真不同：那里 n、p 是 phi 的函数，
        #   雅可比含 q h_i (p+n)/V_T 项。）
        b_diag = bb0.copy()
        a_low = aa.copy(); a_low[0] = 0.0
        c_up = cc.copy(); c_up[-1] = 0.0

        delta = thomas_solve(a_low, b_diag, c_up, -F)
        delta = np.clip(delta, -clamp, clamp)

        # ---- 阻尼回溯线搜索（Armijo）----
        alpha = 1.0
        normF_try = normF
        if damped:
            for _ in range(max(1, int(max_backtrack))):
                phi_try = phi.copy()
                phi_try[1:-1] += alpha * delta
                F_try = _poisson_residual(phi_try, n, p, Na_arr, Nd_arr,
                                          aa, bb0, cc, h_i)
                normF_try = np.linalg.norm(F_try)
                if normF_try <= (1.0 - 2e-4 * alpha) * normF:
                    break
                alpha *= 0.5
            # 残差地板豁免（同第一阶段）：离散残差为大项相减，
            # 低于 1e-8 相对量级时受 double 舍入主导
            src_scale = Q_E * max(np.max(Na_arr), np.max(Nd_arr)) \
                * np.mean(h_i) * np.sqrt(max(len(phi) - 2, 1))
            if (normF_try > (1.0 - 2e-4 * alpha) * normF
                    and normF > 1e-8 * src_scale and warn):
                warnings.warn(
                    f"solve_poisson_pn: Armijo 线搜索折半 "
                    f"{max(1, int(max_backtrack))} 次仍未满足充分下降条件 "
                    f"(alpha={alpha:.2e})，该步接受最小步长解。",
                    RuntimeWarning, stacklevel=2)
            phi, F, normF = phi_try, F_try, normF_try
            dmove = alpha * np.max(np.abs(delta))
        else:
            phi[1:-1] += delta
            F = _poisson_residual(phi, n, p, Na_arr, Nd_arr,
                                  aa, bb0, cc, h_i)
            normF = np.linalg.norm(F)
            dmove = np.max(np.abs(delta))

        phi[0] = bc["phi_left"]
        phi[-1] = bc["phi_right"]

        if dmove < tol:
            converged = True
            break

    if not converged and warn:
        warnings.warn(
            f"solve_poisson_pn: Newton 未在 {max_iter} 次迭代内收敛，"
            f"返回末次迭代解。", RuntimeWarning, stacklevel=2)

    return phi, niter, converged


# ================================================================================


def _newton_system(phi, n, p, x, Na_arr, Nd_arr, hw, h_i):
    """
    组装全 Newton 系统的残差 F 与块三对角雅可比。

    未知量按节点分块 u_i = [phi_i, n_i, p_i]（内部节点 i = 1..N-2）。
    方程组：
        F_phi = 电位移通量差 - q h_i (p - n + N_d - N_a)
        F_n   = J_n,i+1/2 - J_n,i-1/2 - q h_i U
        F_p   = J_p,i+1/2 - J_p,i-1/2 + q h_i U

    雅可比块的推导要点：
        J_n,i+1/2 = (q D_n/h)[B(delta) n_{i+1} - B(-delta) n_i],
            delta = (phi_i - phi_{i+1})/V_T
        dB/ddelta = B(delta)(1 - B(-delta))/delta   （见 bernoulli_deriv）
        U 的线性化保留 dU/dn、dU/dp 一阶项（半隐式，保证对角占优）。

    返回 F (N-2, 3)，A_diag/A_low/A_up (N-2, 3, 3)。
    """
    N = len(x)
    M = N - 2

    # ---- 界面量（长度 N-1）----
    delta = (phi[:-1] - phi[1:]) / VT
    Bp = bernoulli(delta)
    Bm = bernoulli(-delta)
    Bpd = bernoulli_deriv(delta)
    Bmd = bernoulli_deriv(-delta)

    # 电子/空穴电流（界面 i+1/2）
    Jn = (Q_E * D_N / hw) * (Bp * n[1:] - Bm * n[:-1])
    Jp = -(Q_E * D_P / hw) * (Bm * p[1:] - Bp * p[:-1])

    U_star, dU_dn, dU_dp = srh_rate(n, p)

    # ---- 残差（内部节点）----
    lap = (-EPS_SI / hw[:-1]) * phi[:-2] + \
        (EPS_SI / hw[:-1] + EPS_SI / hw[1:]) * phi[1:-1] + \
        (-EPS_SI / hw[1:]) * phi[2:]
    F_phi = lap - Q_E * h_i * (p[1:-1] - n[1:-1]
                               + Nd_arr[1:-1] - Na_arr[1:-1])
    F_n = Jn[1:] - Jn[:-1] - Q_E * h_i * U_star[1:-1]
    F_p = Jp[1:] - Jp[:-1] + Q_E * h_i * U_star[1:-1]
    F = np.stack([F_phi, F_n, F_p], axis=1)

    # ---- 雅可比块 ----
    qDn_h = Q_E * D_N / hw
    qDp_h = Q_E * D_P / hw
    qh = Q_E * h_i

    A_diag = np.zeros((M, 3, 3))
    A_low = np.zeros((M, 3, 3))
    A_up = np.zeros((M, 3, 3))

    # 预计算各内部节点左右界面的导数耦合量
    # 节点 i=j+1：左界面 hw 索引 j，右界面 hw 索引 j+1
    # 电子：dJ_n/dphi_i（左节点）= (qD_n/h)/V_T [B'(d) n_{i+1} + B'(-d) n_i]
    #       （和式！FD 校验见 test_jacobian.py）
    dphi_dn_L = (qDn_h[:-1] / VT) * (Bpd[:-1] * n[1:-1] + Bmd[:-1] * n[:-2])
    dphi_dn_R = (qDn_h[1:] / VT) * (Bpd[1:] * n[2:] + Bmd[1:] * n[1:-1])
    # 空穴：J_p 的 SG 带总负号，但 dJ_p/dphi 的结构与电子完全相同
    #   dJ_p/dphi_k（左节点）= (qD_p/h)/V_T [B'(-d) p_{k+1} + B'(d) p_k]
    # 注意必须是和式：准中性区 p 为常数、d->0 时，差式
    # B'(-d)p_{k+1} - B'(d)p_k 会因大数相消精确为 0，而真实导数为
    # ~qD_p p/(h V_T) ~ 1e11（FD 校验实锤，见 test_jacobian.py）。
    dphi_dp_L = (qDp_h[:-1] / VT) * (Bmd[:-1] * p[1:-1] + Bpd[:-1] * p[:-2])
    dphi_dp_R = (qDp_h[1:] / VT) * (Bmd[1:] * p[2:] + Bpd[1:] * p[1:-1])

    # 主对角块
    A_diag[:, 0, 0] = EPS_SI / hw[:-1] + EPS_SI / hw[1:]
    A_diag[:, 0, 1] = qh
    A_diag[:, 0, 2] = -qh
    A_diag[:, 1, 0] = dphi_dn_R + dphi_dn_L
    A_diag[:, 1, 1] = -(qDn_h[1:] * Bm[1:] + qDn_h[:-1] * Bp[:-1]) \
        - qh * dU_dn[1:-1]
    A_diag[:, 1, 2] = -qh * dU_dp[1:-1]
    A_diag[:, 2, 0] = dphi_dp_R + dphi_dp_L
    A_diag[:, 2, 1] = qh * dU_dn[1:-1]
    A_diag[:, 2, 2] = qDp_h[1:] * Bp[1:] + qDp_h[:-1] * Bm[:-1] \
        + qh * dU_dp[1:-1]

    # 次对角块（节点 i 对 i-1）
    A_low[:, 0, 0] = -EPS_SI / hw[:-1]
    A_low[:, 1, 0] = -dphi_dn_L
    A_low[:, 1, 1] = qDn_h[:-1] * Bm[:-1]
    A_low[:, 2, 0] = -dphi_dp_L
    A_low[:, 2, 2] = -qDp_h[:-1] * Bp[:-1]

    # 超对角块（节点 i 对 i+1）
    A_up[:, 0, 0] = -EPS_SI / hw[1:]
    A_up[:, 1, 0] = -dphi_dn_R
    A_up[:, 1, 1] = qDn_h[1:] * Bp[1:]
    A_up[:, 2, 0] = -dphi_dp_R
    A_up[:, 2, 2] = -qDp_h[1:] * Bm[1:]

    # ---- 行缩放由调用方施加 ----
    # 此处返回原始（未缩放）残差与雅可比。缩放方式取决于未知量选型：
    #   - (phi, n, p) 变量：旧版按 qD/h 固定尺度（条件数 ~1e36，已弃用）；
    #   - (phi, phi_n, phi_p) 准费米势变量（现主求解器）：逐节点
    #     n_i/V_T、p_i/V_T 尺度，由 newton_solve.assemble 在列变换后施加。
    return F, A_diag, A_low, A_up


def _qfp_assemble(phi, phi_n, phi_p, x, Na_arr, Nd_arr, hw, h_i,
                  n_L, p_L, n_R, p_R):
    """
    组装准费米势变量下的（行缩放后）残差 F 与块三对角雅可比。

    步骤：
        1) 由 (phi, phi_n, phi_p) 经玻尔兹曼统计还原 (n, p)（边界钉扎）；
        2) _newton_system 组装 (phi, n, p) 变量下的原始残差与雅可比
           （已经 test_jacobian.py 有限差分校验）；
        3) _qfp_columns 列变换到 (phi, phi_n, phi_p) 变量；
        4) 逐节点行缩放：
             F_phi 行 / (eps * inv_hsum)
             F_n   行 / (q D_n n_i inv_hsum / V_T)
             F_p   行 / (q D_p p_i inv_hsum / V_T)
           使每个 3x3 块内元素量级归一（块消元不损失精度）。

    返回 (F, A_diag, A_low, A_up)，均为缩放后；(M,3)/(M,3,3)。
    """
    N = len(x)
    # ---- 1) 载流子还原（指数限幅防溢出 + 接触钉扎）----
    nn = NI * np.exp(np.clip((phi - phi_n) / VT, -150.0, 150.0))
    pp = NI * np.exp(np.clip((phi_p - phi) / VT, -150.0, 150.0))
    nn[0], pp[0], nn[-1], pp[-1] = n_L, p_L, n_R, p_R

    # ---- 2) 原始残差与雅可比（(phi, n, p) 变量）----
    F, Ad, Al, Au = _newton_system(phi, nn, pp, x, Na_arr, Nd_arr, hw, h_i)

    # ---- 3) 列变换 -> (phi, phi_n, phi_p) ----
    Ad = qfp_columns(Ad, nn[1:-1] / VT, pp[1:-1] / VT)
    Al = qfp_columns(Al, nn[:-2] / VT, pp[:-2] / VT)
    Au = qfp_columns(Au, nn[2:] / VT, pp[2:] / VT)

    # ---- 4) 逐节点行缩放 ----
    inv_hsum = 1.0 / (1.0 / hw[:-1] + 1.0 / hw[1:])
    n_eff = np.maximum(nn[1:-1], 1e8)
    p_eff = np.maximum(pp[1:-1], 1e8)
    s_phi = inv_hsum / EPS_SI
    s_n = inv_hsum * VT / (Q_E * D_N * n_eff)
    s_p = inv_hsum * VT / (Q_E * D_P * p_eff)
    scale = np.stack([s_phi, s_n, s_p], axis=1)          # (M,3)
    return (F * scale, Ad * scale[:, :, None],
            Al * scale[:, :, None], Au * scale[:, :, None])


def newton_solve(x, Na, Nd, V_a, phi=None, n=None, p=None,
                 tol=1e-9, max_iter=30, clamp=0.5, verbose=False):
    """
    全 Newton 求解 PN 结稳态（主求解器）。

    未知量按节点分块 u_i = [phi_i, phi_n,i, phi_p,i]（准费米势变量），
    载流子由玻尔兹曼统计还原：
        n_i = n_i exp((phi_i - phi_n,i)/V_T),  p_i = n_i exp((phi_p,i - phi_i)/V_T)
    残差与雅可比仍在 (phi, n, p) 下组装（_newton_system，已经
    test_jacobian.py 有限差分校验），再经 _qfp_columns 列变换到
    准费米势变量。该变换把指数耦合解析吸收，块矩阵条件数从 ~1e36
    降到 O(1)，Newton 二次收敛（实测平衡态 3~5 步，偏压延续 2~3 步）。

    行缩放（变换后）：
        F_phi 行 / (eps * inv_hsum)          —— 电位移尺度
        F_n   行 / (q D_n n_i inv_hsum/V_T)  —— 电子电流尺度（逐节点）
        F_p   行 / (q D_p p_i inv_hsum/V_T)  —— 空穴电流尺度（逐节点）
    使每个 3x3 块内元素量级归一，块消元不损失精度。

    参数
    ----
    phi, n, p : 初值（None 时自动生成；给 (phi,n,p) 时反推准费米势，
                用于参数延续）
    tol       : 收敛判据（步长 |Delta u|_inf，单位 V）
    clamp     : 单 Newton 步限幅 (V)，防 exp 溢出

    返回 dict：phi, n, p, iters, converged
    """
    Na_arr = np.where(x >= 0.0, 0.0, Na)
    Nd_arr = np.where(x >= 0.0, Nd, 0.0)
    bc = contact_values(V_a, Na, Nd)
    cc = dict(Na=Na, Nd=Nd)
    n_L, p_L, n_R, p_R = contact_carriers(bc, Na, Nd)

    N = len(x)
    hw = np.diff(x)
    h_i = 0.5 * (hw[:-1] + hw[1:])

    # ---- 初值（未知量：phi, phi_n, phi_p）----
    if phi is None:
        # 平衡态初值：耗尽近似解析形状（见 initial_equilibrium_guess）
        phi = initial_equilibrium_guess(x, Na, Nd)
    phi = np.asarray(phi, dtype=float).copy()
    phi[0], phi[-1] = bc["phi_left"], bc["phi_right"]
    if n is None or p is None:
        # 非零偏压冷启动初值：准费米势与电势同步加"偏压斜坡"——
        # p 侧准中性区整体平移 V_a，经耗尽区（宽度 W）平滑过渡到
        # n 侧的 0。若改用全域线性插值，器件中部会出现 n=p~n_i 的
        # 假载流子（残差 ~7， scaled），把大残差注入低载流子节点，
        # 行缩放地板处行近奇异 -> Newton 方向爆炸（dphi_p ~1e4 V）。
        # 斜坡初值下 (phi-phi_n)/VT 保持平衡形状，载流子剖面物理，
        # 残差仅为电流发散量级。
        V_shift = bc["phi_n_L"] - bc["phi_n_R"]      # = V_a
        W_ramp = max(W_analytic(Na, Nd, 0.0), 4.0 * np.max(hw))
        t = np.clip((0.5 * W_ramp - x) / W_ramp, 0.0, 1.0)
        ramp = 0.5 * (1.0 - np.cos(np.pi * t))       # 1(p侧) -> 0(n侧)
        phi = phi + V_shift * ramp
        phi_n = np.full_like(x, bc["phi_n_R"]) + V_shift * ramp
        phi_p = np.full_like(x, bc["phi_p_R"]) + V_shift * ramp
    else:
        # 参数延续：由上一偏压解的 (phi, n, p) 反推准费米势
        n = np.asarray(n, dtype=float)
        p = np.asarray(p, dtype=float)
        phi_n = phi - VT * np.log(np.maximum(n, 1.0) / NI)
        phi_p = phi + VT * np.log(np.maximum(p, 1.0) / NI)
    phi_n = np.asarray(phi_n, dtype=float).copy()
    phi_p = np.asarray(phi_p, dtype=float).copy()
    phi_n[0], phi_n[-1] = bc["phi_n_L"], bc["phi_n_R"]
    phi_p[0], phi_p[-1] = bc["phi_p_L"], bc["phi_p_R"]

    # 电势安全限幅（防 exp 溢出；延续求解下几乎不触发）
    # 电势/准费米势安全限幅（防 exp 溢出；延续求解下几乎不触发）。
    # 注意：
    #   1) 限幅范围必须包住边界值本身——曾用硬编码 p_lo=-1.2，V_a<=-0.72
    #      时 phi_left < p_lo，边界条件被裁坏（线搜索在错误问题上评估）；
    #   2) 只裁内部节点——边界节点已钉为 Dirichlet 值，裁它等于改边界条件。
    q_lo, q_hi = min(V_a, 0.0) - 0.6, max(V_a, 0.0) + 0.6
    p_lo = min(bc["phi_left"], bc["phi_right"]) - 0.8
    p_hi = max(bc["phi_left"], bc["phi_right"]) + 0.8
    phi[1:-1] = np.clip(phi[1:-1], p_lo, p_hi)
    phi_n[1:-1] = np.clip(phi_n[1:-1], q_lo, q_hi)
    phi_p[1:-1] = np.clip(phi_p[1:-1], q_lo, q_hi)

    def carriers(phi_, phin_, phip_):
        """由电势/准费米势还原载流子浓度（指数限幅 + 边界钉扎）。"""
        nn = NI * np.exp(np.clip((phi_ - phin_) / VT, -150.0, 150.0))
        pp = NI * np.exp(np.clip((phip_ - phi_) / VT, -150.0, 150.0))
        nn[0], pp[0], nn[-1], pp[-1] = n_L, p_L, n_R, p_R
        return nn, pp

    def assemble(phi_, phin_, phip_):
        """组装（行缩放后）残差 F 与块三对角雅可比（QFP 变量）。"""
        return _qfp_assemble(phi_, phin_, phip_, x, Na_arr, Nd_arr,
                             hw, h_i, n_L, p_L, n_R, p_R)

    F, A_diag, A_low, A_up = assemble(phi, phi_n, phi_p)
    normF0 = np.linalg.norm(F)
    normF = normF0

    converged = False
    niter = 0
    for niter in range(1, max_iter + 1):
        rhs = -F.copy()
        # 边界处理：Dirichlet 节点（0 与 N-1）的值在整个迭代中固定，
        # 即 du_bnd = 0；且 F 组装时已代入真实边界值（phi[0]=phi_left 等）。
        # 消去边界块后右端仅为 -F，不能再扣减 A_low[0] @ bnd_q——
        # 该扣减只在"F 以零边界值组装、再把边界值当未知量提升"的
        # formulation 下成立。此处若扣减，会把 ~0.18 的假残差注入
        # 首末内部节点，经 Laplace 型耦合放大成 O(0.3 V) 的假步长
        # （实测：max|dphi| 从 ~1e-4 V 假升至 0.29 V，Armijo 失效）。
        A_low[0] = 0.0
        A_up[-1] = 0.0

        du = block_thomas(A_diag, A_low, A_up, rhs)
        du_max = float(np.abs(du).max())
        sc = min(1.0, clamp / max(du_max, 1e-30)) if du_max > clamp else 1.0

        # Armijo 回溯（在缩放残差范数上）
        alpha = 1.0
        normF_try = normF
        for _ in range(24):
            phi_t = phi.copy()
            phin_t = phi_n.copy()
            phip_t = phi_p.copy()
            phi_t[1:-1] += alpha * sc * du[:, 0]
            phin_t[1:-1] += alpha * sc * du[:, 1]
            phip_t[1:-1] += alpha * sc * du[:, 2]
            phi_t[0], phi_t[-1] = bc["phi_left"], bc["phi_right"]
            phin_t[0], phin_t[-1] = bc["phi_n_L"], bc["phi_n_R"]
            phip_t[0], phip_t[-1] = bc["phi_p_L"], bc["phi_p_R"]
            # 同样只裁内部（边界已钉，见初始化处注释）
            phi_t[1:-1] = np.clip(phi_t[1:-1], p_lo, p_hi)
            phin_t[1:-1] = np.clip(phin_t[1:-1], q_lo, q_hi)
            phip_t[1:-1] = np.clip(phip_t[1:-1], q_lo, q_hi)
            F_t, Ad_t, Al_t, Au_t = assemble(phi_t, phin_t, phip_t)
            normF_try = np.linalg.norm(F_t)
            if normF_try <= (1.0 - 2e-4 * alpha) * normF:
                break
            alpha *= 0.5

        phi, phi_n, phi_p = phi_t, phin_t, phip_t
        F, A_diag, A_low, A_up = F_t, Ad_t, Al_t, Au_t
        dmove = alpha * sc * du_max
        normF_prev = normF
        normF = normF_try

        if verbose:
            print(f"    Newton {niter:3d}: ||F||={normF:.3e} "
                  f"(rel {normF/max(normF0,1e-300):.2e}), "
                  f"dphi={dmove:.3e}, alpha={alpha:.2e}")
        # 收敛判据（二者满足其一）：
        #   a) 步长进 tol 且残差相对下降 8 个数量级；
        #   b) 步长已机器精度级、残差停在组装噪声底（F 含 O(1)~O(100)
        #      项相消，||F|| 无法低于 ~1e-9；此时电荷密度误差已在
        #      ~1e-16 相对量级，物理上完全收敛，继续迭代只会空转）。
        #   注意 b 必须加残差绝对底限 1e-6：死锁态（线搜索 24 次全败、
        #   alpha~6e-8）同样满足"步长极小 + 残差不降"，若无底限会被
        #   误判为收敛（实测反偏支曾因此冻结在垃圾解上）。
        if dmove < tol and (normF < 1e-8 * max(normF0, 1e-300)
                            or (normF >= normF_prev * (1.0 - 1e-2)
                                and normF < 1e-6)):
            converged = True
            break

    if not converged:
        warnings.warn(
            f"newton_solve: Newton 未在 {max_iter} 步内收敛 "
            f"(V_a={V_a:.3f} V, ||F||={normF:.3e})。",
            RuntimeWarning, stacklevel=2)

    n, p = carriers(phi, phi_n, phi_p)
    return dict(phi=phi, n=n, p=p, iters=niter, converged=converged)


def solve_bias(x, Na, Nd, V_a, state=None, **kw):
    """
    单个偏压点求解（参数延续入口）。

    state = (phi, n, p) 为上一步解，作为初值传入；None 时冷启动。

    冷启动的鲁棒性兜底：直接冷启动（斜坡初值）在部分偏压点会死锁——
    斜坡与真解的准费米过渡区不匹配时，首步把贴边节点载流子压垮至
    行缩放地板以下，电子/空穴行近奇异，Newton 方向爆炸（dphi_p ~1e4
    V），Armijo 回溯 24 次全线失败。此时自动改用参数延续：先解平衡态，
    再以 <=0.1 V 步长分步推进到 V_a（实测全偏压范围收敛）。
    """
    if state is None:
        res = newton_solve(x, Na, Nd, V_a, **kw)
        if res["converged"] or abs(V_a) < 1e-12:
            return res
        # ---- 冷启动失败：参数延续兜底 ----
        warnings.warn(
            f"solve_bias: V_a={V_a:.3f} V 冷启动未收敛，改用参数延续。",
            RuntimeWarning, stacklevel=2)
        core = {k: v for k, v in kw.items()
                if k in ("tol", "max_iter", "clamp")}
        res_eq = newton_solve(x, Na, Nd, 0.0, **core)
        if not res_eq["converged"]:
            return res
        tol = core.get("tol", 1e-9)
        clamp = core.get("clamp", 0.5)
        max_iter = core.get("max_iter", 30)
        st = (res_eq["phi"], res_eq["n"], res_eq["p"])
        total = res_eq["iters"]
        n_sub = max(1, int(np.ceil(abs(V_a) / 0.1)))
        ok = True
        for k in range(1, n_sub + 1):
            V_prev = V_a * (k - 1) / n_sub
            V_next = V_a * k / n_sub
            r = newton_solve(x, Na, Nd, V_next, phi=st[0], n=st[1],
                             p=st[2], tol=tol, clamp=clamp, max_iter=max_iter)
            total += r["iters"]
            # 单步失败：区间减半重试（最多 4 轮）
            for _ in range(4):
                if r["converged"]:
                    break
                V_mid = 0.5 * (V_prev + V_next)
                rm = newton_solve(x, Na, Nd, V_mid, phi=st[0], n=st[1],
                                  p=st[2], tol=tol, clamp=clamp,
                                  max_iter=max_iter)
                total += rm["iters"]
                if not rm["converged"]:
                    r = rm
                    break
                st = (rm["phi"], rm["n"], rm["p"])
                r = newton_solve(x, Na, Nd, V_next, phi=st[0], n=st[1],
                                 p=st[2], tol=tol, clamp=clamp,
                                 max_iter=max_iter)
                total += r["iters"]
                V_next = V_mid
            if not r["converged"]:
                ok = False
                break
            st = (r["phi"], r["n"], r["p"])
        return dict(phi=st[0], n=st[1], p=st[2], iters=total, converged=ok)
    return newton_solve(x, Na, Nd, V_a, phi=state[0], n=state[1],
                        p=state[2], **kw)

def initial_equilibrium_guess(x, Na, Nd):
    """
    平衡态初值：耗尽近似解析形状（保证 Gummel 迭代快速收敛）。

        V_bi = V_T ln(N_a N_d/n_i^2),  W = sqrt(2 eps V_bi/q (1/N_a+1/N_d))
        p 侧耗尽区 (-W/2 <= x <= 0): phi = phi_L + (q N_a/2eps)(x + W/2)^2
        n 侧耗尽区 (0 <= x <= W/2) : phi = phi_L + V_bi - (q N_d/2eps)(W/2 - x)^2
        远场 (|x| > W/2)            ：常数 phi_L / phi_R（准中性区无压降）

    注意抛物线公式只在耗尽区内有效：远场若外推抛物线将发散
    （二次项系数 qN/2eps ~ 7.7e12 V/m^2，数微米外即达数十伏），
    必须切换为常数。
    """
    V_bi = V_bi_analytic(Na, Nd)
    W = W_analytic(Na, Nd, 0.0)
    phi_L = -VT * np.log(Na / NI)
    phi_R = VT * np.log(Nd / NI)
    wh = 0.5 * W

    phi = np.empty_like(x)
    far_p = x <= -wh
    far_n = x >= wh
    mid = ~(far_p | far_n)
    phi[far_p] = phi_L
    phi[far_n] = phi_R
    xm = x[mid]
    phi[mid] = np.where(
        xm < 0.0,
        phi_L + (Q_E * Na / (2.0 * EPS_SI)) * (xm + wh) ** 2,
        phi_L + V_bi - (Q_E * Nd / (2.0 * EPS_SI)) * (wh - xm) ** 2)
    return phi


def gummel_solve(x, Na, Nd, V_a, phi=None, n=None, p=None,
                 tol=1e-9, max_outer=200, poisson_kw=None, verbose=False,
                 omega=0.5):
    """
    给定偏压 V_a，用 Gummel 解耦迭代求解 (phi, n, p)。

    .. deprecated::
        主求解器已换为全 Newton（newton_solve，准费米势变量）。
        实测：本器件 SRH 复合下解耦迭代停滞——以平衡态为初值，
        原始迭代（omega=1）外层 dphi 迅速卡在 ~9e-2 V 不再下降，
        载流子浓度剧烈振荡（max(n) 在 N_d ~ 49 N_d 间跳动）；
        欠松弛 omega=0.5 把 dphi 压到 ~7e-3 V 但仍不收敛
        （200 步外层迭代未达 tol=1e-9）。原因是连续性方程的解
        对 phi 呈指数依赖，而泊松方程又把 phi 拉回，解耦迭代的
        收敛半径极小。本函数保留作方法对照与教学演示，
        pipeline 中不再调用。

    迭代结构（外循环）：
        1) 电子连续性（三对角，p 与线性化系数固定）；
        2) 空穴连续性（三对角，n 固定）；
        3) 泊松方程（线性，n、p 固定，Newton 一步收敛）；
        4) 收敛判据：|Delta phi|_inf < tol (V)。

    欠松弛（damping）
    ----------------
    连续性方程的解对 phi 呈指数依赖（n = n_i exp(phi/V_T)），而泊松
    方程又把 phi 拉回，二者的强耦合使原始 Gummel 迭代收敛半径极小。
    实测（V_a=0 平衡态，200 步外层迭代）：omega=1 时 dphi 卡在
    ~9e-2 V，max(n) 在 N_d ~ 49 N_d 间剧烈跳动（停滞 + 振荡）；
    omega=0.5 时 dphi 降至 ~7e-3 V 但仍远高于 tol=1e-9。
    对载流子更新引入松弛因子 omega：
        n <- (1-omega) n_old + omega n_new
    omega=0.5 只能减缓而不能消除停滞；全 Newton（newton_solve）
    在同一初值下 6 步收敛，二者对照见 make_figures.py 图 8。

    参数
    ----
    phi, n, p : 初值（None 时自动生成：平衡近似 + 准费米势线性插值）
    tol       : 外层收敛判据 (V)
    omega     : 欠松弛因子 (0, 1]

    返回 dict：phi, n, p, iters, converged
    """
    Na_arr = np.full_like(x, Na)
    Nd_arr = np.full_like(x, Nd)
    Na_arr[x >= 0.0] = 0.0
    Nd_arr[x < 0.0] = 0.0
    bc = contact_values(V_a, Na, Nd)
    cc = dict(Na=Na, Nd=Nd)

    # ---- 初值 ----
    if phi is None:
        # 耗尽近似初值 + 按偏压平移（p 侧接触整体移动 V_a）
        phi = initial_equilibrium_guess(x, Na, Nd)
        phi = phi + (bc["phi_left"] - (-VT * np.log(Na / NI)))
    phi = phi.astype(float).copy()
    phi[0], phi[-1] = bc["phi_left"], bc["phi_right"]

    if n is None or p is None:
        # 准费米势线性插值初值
        s = (x - x[0]) / (x[-1] - x[0])
        phi_n = bc["phi_n_L"] * (1.0 - s) + bc["phi_n_R"] * s
        phi_p = bc["phi_p_L"] * (1.0 - s) + bc["phi_p_R"] * s
        n, p = carrier_from_phi(phi, phi_n, phi_p)
    n = n.astype(float).copy()
    p = p.astype(float).copy()
    n[0], p[0], n[-1], p[-1] = contact_carriers(bc, Na, Nd)

    pk = dict(tol=1e-11, max_iter=100, clamp=0.3, damped=True)
    if poisson_kw:
        pk.update(poisson_kw)

    converged = False
    niter = 0
    total_newton = 0
    for niter in range(1, max_outer + 1):
        # 1) 电子连续性 + 欠松弛
        n_s = solve_electron_continuity(phi, p, n, x, bc, cc)
        n_new = (1.0 - omega) * n + omega * n_s
        # 2) 空穴连续性（用最新 n）+ 欠松弛
        p_s = solve_hole_continuity(phi, n_new, p, x, bc, cc)
        p_new = (1.0 - omega) * p + omega * p_s
        # 3) 泊松 + 欠松弛（压 2-cycle：松弛后谱半径 1-omega+omega*g'，
        #    对 g'=-1 的二周期模态，omega=0.5 时精确归零）
        phi_s, n_newton, _ = solve_poisson_pn(
            phi, n_new, p_new, Na_arr, Nd_arr, x, bc, **pk)
        phi_new = (1.0 - omega) * phi + omega * phi_s
        total_newton += n_newton

        dphi = np.max(np.abs(phi_new - phi))
        dn = np.max(np.abs(n_new - n)) / max(np.max(n_new), 1.0)
        dp = np.max(np.abs(p_new - p)) / max(np.max(p_new), 1.0)

        phi, n, p = phi_new, n_new, p_new

        if verbose and niter % 10 == 0:
            print(f"    Gummel iter {niter:3d}: dphi={dphi:.3e} V, "
                  f"dn={dn:.3e}, dp={dp:.3e}")
        if dphi < tol and niter > 1:
            converged = True
            break

    if not converged:
        warnings.warn(
            f"gummel_solve: Gummel 迭代未在 {max_outer} 次外层迭代内收敛 "
            f"(V_a={V_a:.3f} V, dphi={dphi:.3e})。",
            RuntimeWarning, stacklevel=2)

    return dict(phi=phi, n=n, p=p, iters=niter,
                newton_iters=total_newton, converged=converged,
                dphi=dphi)
