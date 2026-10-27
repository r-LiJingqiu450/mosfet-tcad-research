# -*- coding: utf-8 -*-
"""
stage2_pn/tests/test_stage2.py —— PN 结全 Newton 系统测试。

由原件 pn_junction/test_qfp.py 与 test_jacobian.py 合并迁移：
函数级导入替换 pj. 前缀，_qfp_columns/_block_thomas 已升入 common/。

覆盖：
    1. 块 Thomas vs 稠密解（common.thomas，另见 tests/test_common.py）
    2. 原始 (phi,n,p) 雅可比 FD 校验（粗网格）
    3. QFP 列变换链式法则（跨 14 个数量级浓度）
    4. 行缩放透明性（相似变换不改变 Newton 步）
    5. 条件数：QFP 变换 + 行缩放后 O(1e3)
    6. 冷启动首步 Newton 降低缩放残差
    7. 端到端：平衡态 V_bi / np=n_i^2
    8. 端到端：正偏注入电流 = 总复合（离散恒等式）
    9. 端到端：反偏产生电流与耗尽层展宽
    10. 变换后系统在收敛正偏解上的 FD 校验
"""

import numpy as np
import pytest

from semiconductor_sim.common import constants as C
from semiconductor_sim.common.grid import make_junction_grid
from semiconductor_sim.common.qfp import qfp_columns
from semiconductor_sim.common.thomas import block_thomas
from semiconductor_sim.stage2_pn.solver import (
    carrier_from_phi, contact_carriers, contact_values, _newton_system,
    _qfp_assemble, newton_solve, srh_rate)
from semiconductor_sim.stage2_pn.postprocess import current_density, depletion_width
from semiconductor_sim.stage2_pn.run import solve_equilibrium

# 粗网格：FD 循环与端到端求解都快（细网格 FD ~24 s / 求解 ~1 s）
GRID = dict(L=2.0e-6, d0=20e-9, transition=0.15e-6, h_max=40e-9)


# ============================================================================
# 初态构造（迁移自原 test_qfp.build_state / test_jacobian.build_state）
# ============================================================================

def build_state(Na, Nd, V_a=0.0, grid=None):
    """平衡态 QFP 初态：(phi, phi_n, phi_p)，n/p 由玻尔兹曼统计还原。"""
    g = dict(GRID if grid is None else grid)
    x = make_junction_grid(**g)
    bc = contact_values(V_a, Na, Nd)
    from semiconductor_sim.stage2_pn.solver import initial_equilibrium_guess
    phi = initial_equilibrium_guess(x, Na, Nd)
    phi = phi + (bc["phi_left"] - (-C.VT * np.log(Na / C.NI)))
    phi[0], phi[-1] = bc["phi_left"], bc["phi_right"]
    s = (x - x[0]) / (x[-1] - x[0])
    phi_n = bc["phi_n_L"] * (1.0 - s) + bc["phi_n_R"] * s
    phi_p = bc["phi_p_L"] * (1.0 - s) + bc["phi_p_R"] * s
    return x, phi, phi_n, phi_p, bc


def build_state_np(Na, Nd, V_a=0.0, grid=None):
    """原始 (phi, n, p) 初态（迁移自原 test_jacobian.build_state）。"""
    g = dict(GRID if grid is None else grid)
    x = make_junction_grid(**g)
    bc = contact_values(V_a, Na, Nd)
    from semiconductor_sim.stage2_pn.solver import initial_equilibrium_guess
    phi = initial_equilibrium_guess(x, Na, Nd)
    phi = phi + (bc["phi_left"] - (-C.VT * np.log(Na / C.NI)))
    phi[0], phi[-1] = bc["phi_left"], bc["phi_right"]
    s = (x - x[0]) / (x[-1] - x[0])
    phi_n = bc["phi_n_L"] * (1.0 - s) + bc["phi_n_R"] * s
    phi_p = bc["phi_p_L"] * (1.0 - s) + bc["phi_p_R"] * s
    n, p = carrier_from_phi(phi, phi_n, phi_p)
    n[0], p[0], n[-1], p[-1] = contact_carriers(bc, Na, Nd)
    return x, phi, n, p


def raw_transformed(phi, phi_n, phi_p, x, bc, Na, Nd):
    """(phi,phi_n,phi_p) -> (n,p) -> 原始残差 + 列变换后的原始雅可比。"""
    Na_arr = np.where(x >= 0.0, 0.0, Na)
    Nd_arr = np.where(x >= 0.0, Nd, 0.0)
    hw = np.diff(x)
    h_i = 0.5 * (hw[:-1] + hw[1:])
    nn = C.NI * np.exp(np.clip((phi - phi_n) / C.VT, -150, 150))
    pp = C.NI * np.exp(np.clip((phi_p - phi) / C.VT, -150, 150))
    n_L, p_L, n_R, p_R = contact_carriers(bc, Na, Nd)
    nn[0], pp[0], nn[-1], pp[-1] = n_L, p_L, n_R, p_R
    F, Ad, Al, Au = _newton_system(phi, nn, pp, x, Na_arr, Nd_arr, hw, h_i)
    Ad = qfp_columns(Ad, nn[1:-1] / C.VT, pp[1:-1] / C.VT)
    Al = qfp_columns(Al, nn[:-2] / C.VT, pp[:-2] / C.VT)
    Au = qfp_columns(Au, nn[2:] / C.VT, pp[2:] / C.VT)
    return F, Ad, Al, Au, nn, pp


def raw_system_np(phi, n, p, x, Na, Nd):
    """原始 (phi,n,p) 系统残差与雅可比。"""
    Na_arr = np.where(x >= 0.0, 0.0, Na)
    Nd_arr = np.where(x >= 0.0, Nd, 0.0)
    hw = np.diff(x)
    h_i = 0.5 * (hw[:-1] + hw[1:])
    return _newton_system(phi, n, p, x, Na_arr, Nd_arr, hw, h_i)


def solve_step(F, Ad, Al, Au):
    """块三对角求解一步（与 newton_solve 相同的边界代数）。

    Dirichlet 节点值固定（du_bnd = 0）且 F 组装时已代入真实边界值，
    故右端仅为 -F——不能再扣减 A_low[0] @ bnd_q（会注入 ~0.18 假残差，
    经 Laplace 型耦合放大成 O(0.3 V) 假步长，见 newton_solve 注释）。
    """
    rhs = -F.copy()
    Al = Al.copy()
    Al[0] = 0.0
    Au = Au.copy()
    Au[-1] = 0.0
    return block_thomas(Ad, Al, Au, rhs)


def dense(Ad, Al, Au):
    """块三对角 -> 稠密矩阵（条件数计算用）。"""
    M = len(Ad)
    J = np.zeros((3 * M, 3 * M))
    for j in range(M):
        J[3 * j:3 * j + 3, 3 * j:3 * j + 3] = Ad[j]
        if j > 0:
            J[3 * j:3 * j + 3, 3 * (j - 1):3 * (j - 1) + 3] = Al[j]
        if j < M - 1:
            J[3 * j:3 * j + 3, 3 * (j + 1):3 * (j + 1) + 3] = Au[j]
    return J


# ============================================================================
# 1) 原始系统雅可比 FD 校验（迁移自 test_jacobian）
# ============================================================================

def test_analytic_jacobian_vs_fd():
    """原始 (phi,n,p) 系统雅可比的 FD 校验（粗网格）。

    阈值 1e-6：实测最大列归一化误差 ~6e-8，留一个多数量级余量；
    足以捕获符号错误、遗漏项等 gross 错误。
    """
    Na = C.NA_CM3 * 1e6
    Nd = C.ND_CM3 * 1e6
    x, phi, n, p = build_state_np(Na, Nd)
    N = len(x)
    M = N - 2
    F0, Ad, Al, Au = raw_system_np(phi, n, p, x, Na, Nd)
    J = np.zeros((3 * M, 3 * M))
    for j in range(M):
        J[3 * j:3 * j + 3, 3 * j:3 * j + 3] = Ad[j]
        if j > 0:
            J[3 * j:3 * j + 3, 3 * (j - 1):3 * (j - 1) + 3] = Al[j]
        if j < M - 1:
            J[3 * j:3 * j + 3, 3 * (j + 1):3 * (j + 1) + 3] = Au[j]
    names = ["phi", "n", "p"]
    fails = []
    max_rel = 0.0
    for var in range(3):
        base = [phi, n, p][var]
        for i in range(1, N - 1):
            col = 3 * (i - 1) + var
            if var == 0:
                h = 1e-6
            else:
                h = 1e-6 * max(abs(base[i]), 1.0)
            arr_p = base.copy()
            arr_p[i] += h
            arr_m = base.copy()
            arr_m[i] -= h
            args_p = [phi, n, p]
            args_m = [phi, n, p]
            args_p[var] = arr_p
            args_m[var] = arr_m
            Fp, _, _, _ = raw_system_np(args_p[0], args_p[1], args_p[2],
                                        x, Na, Nd)
            Fm, _, _, _ = raw_system_np(args_m[0], args_m[1], args_m[2],
                                        x, Na, Nd)
            Jfd = (Fp - Fm) / (2.0 * h)
            Jana = J[:, col].reshape(M, 3)
            col_norm = max(np.abs(Jana).max(), 1e-300)
            rel = np.abs(Jfd - Jana) / col_norm
            max_rel = max(max_rel, rel.max())
            if rel.max() > 1e-4:
                bad = np.argwhere(rel > 1e-4)
                for (r, c) in bad:
                    fails.append((names[var], i, r, Jana[r, c], Jfd[r, c]))
    assert max_rel < 1e-6, f"雅可比 FD 校验失败: {max_rel:.3e}"
    assert len(fails) == 0


# ============================================================================
# 2) QFP 列变换链式法则（代数单测）
# ============================================================================

def test_qfp_columns_chain_rule():
    """qfp_columns(A, cn, cp) 应等于 A @ du/dv（随机块 + 跨 14 个数量级
    载流子浓度，验证实现与数学定义一致）。"""
    rng = np.random.default_rng(42)
    A = rng.normal(size=(50, 3, 3)) * rng.choice(
        [1e-20, 1.0, 1e10], size=(50, 1, 1))
    cn = 10.0 ** rng.uniform(9, 23, size=50)
    cp = 10.0 ** rng.uniform(9, 23, size=50)
    New = qfp_columns(A, cn, cp)
    Mmat = np.zeros((50, 3, 3))
    Mmat[:, 0, 0] = 1.0
    Mmat[:, 1, 0] = cn
    Mmat[:, 1, 1] = -cn
    Mmat[:, 2, 0] = -cp
    Mmat[:, 2, 2] = cp
    ref = np.einsum("kij,kjl->kil", A, Mmat)
    rel = np.abs(New - ref).max() / np.abs(ref).max()
    assert rel < 1e-12


# ============================================================================
# 3) 行缩放透明性
# ============================================================================

def test_row_scaling_transparent():
    """缩放系统与未缩放系统解出的 Newton 步应相同（行缩放是逐行对角
    相似变换，不改变解）。"""
    Na = C.NA_CM3 * 1e6
    Nd = C.ND_CM3 * 1e6
    x, phi, phi_n, phi_p, bc = build_state(Na, Nd)
    hw = np.diff(x)
    h_i = 0.5 * (hw[:-1] + hw[1:])
    Na_arr = np.where(x >= 0.0, 0.0, Na)
    Nd_arr = np.where(x >= 0.0, Nd, 0.0)
    Fs, Ads, Als, Aus = _qfp_assemble(
        phi, phi_n, phi_p, x, Na_arr, Nd_arr, hw, h_i,
        *contact_carriers(bc, Na, Nd))
    F0, Ad, Al, Au, _, _ = raw_transformed(phi, phi_n, phi_p, x, bc, Na, Nd)
    du_s = solve_step(Fs, Ads, Als, Aus)
    du_u = solve_step(F0, Ad, Al, Au)
    scale = float(np.abs(du_u).max())
    diff = float(np.abs(du_s - du_u).max())
    assert diff < 1e-10 * max(scale, 1e-30)


# ============================================================================
# 4) 条件数
# ============================================================================

def test_condition_number_reduced():
    """(phi,n,p) 直拼条件数 ~1e36；QFP 变换 + 行缩放后应降到 O(1e3)。"""
    Na = C.NA_CM3 * 1e6
    Nd = C.ND_CM3 * 1e6
    x, phi, phi_n, phi_p, bc = build_state(Na, Nd)
    hw = np.diff(x)
    h_i = 0.5 * (hw[:-1] + hw[1:])
    Na_arr = np.where(x >= 0.0, 0.0, Na)
    Nd_arr = np.where(x >= 0.0, Nd, 0.0)
    _, Ad, Al, Au, _, _ = raw_transformed(phi, phi_n, phi_p, x, bc, Na, Nd)
    Fs, Ads, Als, Aus = _qfp_assemble(
        phi, phi_n, phi_p, x, Na_arr, Nd_arr, hw, h_i,
        *contact_carriers(bc, Na, Nd))
    c_un = np.linalg.cond(dense(Ad, Al, Au))
    c_sc = np.linalg.cond(dense(Ads, Als, Aus))
    assert c_sc < 1e8, f"缩放系统条件数仍过大: {c_sc:.3e}"
    assert c_un / c_sc > 1e8, "行缩放未带来数量级改善"


# ============================================================================
# 5) Newton 步方向正确性（冷启动首步）
# ============================================================================

def test_newton_step_reduces_residual():
    """从平衡态初值出发，一个完整 Newton 步应显著降低缩放残差。

    残差必须在缩放范数下度量（newton_solve 的 Armijo 与收敛判据都用
    缩放残差；原始范数混排 C/m^2 与 A/m^2 两种量纲）。
    """
    Na = C.NA_CM3 * 1e6
    Nd = C.ND_CM3 * 1e6
    x, phi, phi_n, phi_p, bc = build_state(Na, Nd)
    hw = np.diff(x)
    h_i = 0.5 * (hw[:-1] + hw[1:])
    Na_arr = np.where(x >= 0.0, 0.0, Na)
    Nd_arr = np.where(x >= 0.0, Nd, 0.0)
    n_L, p_L, n_R, p_R = contact_carriers(bc, Na, Nd)

    def scaled_norm(phi_, pn_, pp_):
        F, _, _, _ = _qfp_assemble(
            phi_, pn_, pp_, x, Na_arr, Nd_arr, hw, h_i, n_L, p_L, n_R, p_R)
        return float(np.linalg.norm(F))

    F, Ad, Al, Au = _qfp_assemble(
        phi, phi_n, phi_p, x, Na_arr, Nd_arr, hw, h_i, n_L, p_L, n_R, p_R)
    du = solve_step(F, Ad, Al, Au)
    phi_t = phi.copy()
    pn_t = phi_n.copy()
    pp_t = phi_p.copy()
    phi_t[1:-1] += du[:, 0]
    pn_t[1:-1] += du[:, 1]
    pp_t[1:-1] += du[:, 2]
    phi_t[0], phi_t[-1] = bc["phi_left"], bc["phi_right"]
    pn_t[0], pn_t[-1] = bc["phi_n_L"], bc["phi_n_R"]
    pp_t[0], pp_t[-1] = bc["phi_p_L"], bc["phi_p_R"]
    n0, n1 = scaled_norm(phi, phi_n, phi_p), scaled_norm(phi_t, pn_t, pp_t)
    assert n1 < 0.5 * n0, "完整 Newton 步未降低缩放残差"


# ============================================================================
# 6) 端到端：平衡态物理解
# ============================================================================

def _equilibrium():
    """粗网格平衡态解（测试 6/7/8 共用，避免重复求解）。"""
    Na = C.NA_CM3 * 1e6
    Nd = C.ND_CM3 * 1e6
    x = make_junction_grid(**GRID)
    return solve_equilibrium(x, Na, Nd, verbose=False), x


def test_equilibrium_matches_analytic():
    """冷启动全 Newton 解平衡态：收敛快、内电势与解析一致、np=n_i^2。"""
    res, x = _equilibrium()
    Vbi_a = C.V_bi_analytic(C.NA_CM3 * 1e6, C.ND_CM3 * 1e6)
    Vbi_n = res["phi"][-1] - res["phi"][0]
    np_err = float(np.max(np.abs(res["n"] * res["p"] - C.NI ** 2)) / C.NI ** 2)
    assert res["converged"]
    assert res["iters"] <= 15, f"收敛过慢: {res['iters']} 步"
    assert abs(Vbi_n / Vbi_a - 1.0) < 1e-9
    assert np_err < 1e-6


# ============================================================================
# 7) 端到端：正偏注入电流 = 总复合
# ============================================================================

def test_forward_current_equals_recombination():
    """V_a=+0.5 V：从平衡态单步延续，验证离散恒等式
        Jn 右接触 - Jn 左接触 = q * sum(h_i U_i)
    （电子从 n 接触注入、在器件内复合的总量；总电流恒定是另一条恒等式，
    不能用来验证复合量）。
    """
    Na = C.NA_CM3 * 1e6
    Nd = C.ND_CM3 * 1e6
    eq, x = _equilibrium()
    hw = np.diff(x)
    h_i = 0.5 * (hw[:-1] + hw[1:])
    r = newton_solve(x, Na, Nd, 0.5, phi=eq["phi"], n=eq["n"], p=eq["p"],
                     tol=1e-9)
    bc = contact_values(0.5, Na, Nd)
    phi = r["phi"].copy()
    phi[0], phi[-1] = bc["phi_left"], bc["phi_right"]
    J, Jt, Jn, Jp = current_density(phi, r["n"], r["p"], x)
    U = srh_rate(r["n"], r["p"])[0]
    qU = C.Q_E * np.sum(U[1:-1] * h_i)
    ratio = (Jn[-1] - Jn[0]) / qU
    assert r["converged"]
    assert J > 0
    assert abs(ratio - 1.0) < 1e-2, f"注入电流与总复合不符: {ratio:.4f}"


# ============================================================================
# 8) 端到端：反偏产生电流
# ============================================================================

def test_reverse_bias_generation_current():
    """V_a=-1.0 V（0.1 V 步进延续）：J<0、|J|=q*int|U|（产生）、
    耗尽层宽于平衡态（W 随反偏增长）。"""
    Na = C.NA_CM3 * 1e6
    Nd = C.ND_CM3 * 1e6
    eq, x = _equilibrium()
    hw = np.diff(x)
    h_i = 0.5 * (hw[:-1] + hw[1:])
    st = (eq["phi"], eq["n"], eq["p"])
    r = None
    for k in range(1, 11):
        Va = -0.1 * k
        r = newton_solve(x, Na, Nd, Va, phi=st[0], n=st[1], p=st[2],
                         tol=1e-9)
        assert r["converged"], f"V_a={Va:.1f} 未收敛"
        st = (r["phi"], r["n"], r["p"])
    bc = contact_values(-1.0, Na, Nd)
    phi = r["phi"].copy()
    phi[0], phi[-1] = bc["phi_left"], bc["phi_right"]
    J, Jt, Jn, Jp = current_density(phi, r["n"], r["p"], x)
    U = srh_rate(r["n"], r["p"])[0]
    qU = C.Q_E * np.sum(U[1:-1] * h_i)
    W = depletion_width(phi, r["n"], r["p"], Na, Nd, x)[0]
    W0 = depletion_width(eq["phi"], eq["n"], eq["p"], Na, Nd, x)[0]
    J_gen = C.Q_E * C.NI * W / (2.0 * C.TAU_N)   # 产生电流解析估计
    ratio = (Jn[-1] - Jn[0]) / qU
    assert J < 0
    assert abs(ratio - 1.0) < 2e-2, f"产生电流恒等式不符: {ratio:.4f}"
    assert W > 1.2 * W0, "反偏耗尽层未显著展宽"
    assert 0.3 < abs(J) / J_gen < 1.0, "产生电流与 q*n_i*W/2tau 估计不符"


# ============================================================================
# 9) 变换后系统的 FD 校验（收敛正偏解，恒等相消已被注入梯度破坏）
# ============================================================================

def test_transformed_jacobian_vs_fd():
    """在 V_a=+0.5 收敛解上对 (phi, phi_n, phi_p) 做 FD，校验
    qfp_columns 后的雅可比。

    为何不用平衡态：平台区玻尔兹曼-界面恒等 B(d)=B(-d)e^d 使 SG 电流
    对 phi 的真实导数恒为 0，FD 只能测到 bernoulli 直接式的舍入噪声。
    正偏解的少子注入梯度破坏恒等，FD 恢复意义；残余误差 ~2e-3 为
    准中性区 FD 固有噪声底。
    """
    Na = C.NA_CM3 * 1e6
    Nd = C.ND_CM3 * 1e6
    eq, x = _equilibrium()
    r = newton_solve(x, Na, Nd, 0.5, phi=eq["phi"], n=eq["n"], p=eq["p"],
                     tol=1e-9)
    assert r["converged"]
    bc = contact_values(0.5, Na, Nd)
    phi = r["phi"].copy()
    phi[0], phi[-1] = bc["phi_left"], bc["phi_right"]
    phi_n = phi - C.VT * np.log(np.maximum(r["n"], 1.0) / C.NI)
    phi_p = phi + C.VT * np.log(np.maximum(r["p"], 1.0) / C.NI)
    phi_n[0], phi_n[-1] = bc["phi_n_L"], bc["phi_n_R"]
    phi_p[0], phi_p[-1] = bc["phi_p_L"], bc["phi_p_R"]

    N = len(x)
    M = N - 2
    _, Ad, Al, Au, _, _ = raw_transformed(phi, phi_n, phi_p, x, bc, Na, Nd)
    J = dense(Ad, Al, Au)
    h = 1e-6
    worst = 0.0
    n_clean = 0
    n_col = 0
    for var in range(3):
        base = [phi, phi_n, phi_p][var]
        for i in range(1, N - 1):
            j = i - 1
            col = 3 * j + var
            ap = base.copy()
            ap[i] += h
            am = base.copy()
            am[i] -= h
            up = [phi, phi_n, phi_p]
            um = [phi, phi_n, phi_p]
            up[var] = ap
            um[var] = am
            Fp = raw_transformed(up[0], up[1], up[2], x, bc, Na, Nd)[0]
            Fm = raw_transformed(um[0], um[1], um[2], x, bc, Na, Nd)[0]
            Jfd = (Fp - Fm) / (2 * h)
            Jana = J[:, col].reshape(M, 3)
            col_norm = max(np.abs(Jana).max(), 1e-300)
            rel = float((np.abs(Jfd - Jana) / col_norm).max())
            worst = max(worst, rel)
            n_clean += rel < 1e-4
            n_col += 1
    assert worst < 5e-3, f"变换后雅可比 FD 校验失败: {worst:.3e}"
    assert n_clean >= 0.8 * n_col, "干净列占比过低（疑恒等相消区）"
