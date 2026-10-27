#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
stage2_pn.make_figures —— 生成 PN 结论文（第二阶段）全部插图。

输出到 results/stage2/（系列统一结果目录，见 docs/spec/reproduce_spec.md）：

    fig0_structure.png      器件结构与掺杂分布示意
    fig1_equilibrium.png    平衡态能带图与电势/电荷分布
    fig2_bias_profiles.png  正偏 (+0.5 V) 与反偏 (-3 V) 的载流子/电场分布
    fig3_jv_curve.png       J-V 特性（半对数）与 Shockley 解析对照
    fig4_quasi_fermi.png    正偏准费米势分裂与少子注入
    fig5_recombination.png  SRH 复合/产生率分布与电流守恒恒等式
    fig6_depletion.png      耗尽层宽度-偏压关系与耗尽近似解析对照
    fig7_newton.png         Newton 收敛历史与雅可比条件数对比
    fig8_gummel.png         Gummel 解耦迭代停滞对照（论文 3.3/4.8 节）

结构（系列图表规范 docs/spec/figure_spec.md）：
    _prepare()   单次计算全部插图所需数据（网格->平衡态->特征点->扫描->历史）
    figN_*()     纯绘图函数，输入 _prepare() 的数据 dict
    make_figures() 编排入口：data 缺省时现场 _prepare()

依赖：numpy, matplotlib, common.*, stage2_pn.solver/postprocess。
运行：python3 -m stage2_pn.make_figures（约 1~2 分钟，主要是 66 点偏压扫描）
"""

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

from common.constants import (D_N, D_P, EG, EPS_SI, NI, NA_CM3, ND_CM3,
                              Q_E, TAU_N, TAU_P, VT,
                              Js_analytic, V_bi_analytic, W_analytic)
from common.grid import make_junction_grid
from common.io import ensure_results_dir
from common.qfp import qfp_columns
from common.thomas import block_thomas
from .solver import (carrier_from_phi, contact_carriers, contact_values,
                     initial_equilibrium_guess, newton_solve, _newton_system,
                     _qfp_assemble, solve_electron_continuity,
                     solve_hole_continuity, solve_poisson_pn, srh_rate)
from .postprocess import current_density
from .run import sweep_bias

# ------------------------------------------------------------------
# 全局绘图风格（与第一阶段一致：STIX 数学字体、英文标注、网格）
# ------------------------------------------------------------------
plt.rcParams.update({
    "figure.dpi": 130, "font.size": 11, "axes.grid": True,
    "grid.alpha": 0.3, "axes.linewidth": 1.0,
    "mathtext.fontset": "stix", "font.family": "STIXGeneral",
    "axes.titlesize": 11.5, "axes.labelsize": 11,
    "legend.fontsize": 9, "xtick.labelsize": 10, "ytick.labelsize": 10,
    "figure.constrained_layout.use": True,
})

C_N = "#1f77b4"     # 电子
C_P = "#d62728"     # 空穴
C_PHI = "#000000"   # 电势
C_E = "#2ca02c"     # 电场
C_ANA = "#7f7f7f"   # 解析曲线


def _outdir(outdir):
    """outdir 缺省解析为 results/stage2/（系列统一结果目录）；显式传入时确保目录存在。"""
    if outdir is None:
        return ensure_results_dir("stage2")
    os.makedirs(outdir, exist_ok=True)
    return outdir


# ------------------------------------------------------------------
# 器件与网格
# ------------------------------------------------------------------
L = 5.0e-6
Na = NA_CM3 * 1e6
Nd = ND_CM3 * 1e6


def solve_steps(x_, Na_, Nd_, V_target, state, v_step=0.1):
    """从 state 出发以 <= v_step 的步长延续到 V_target。"""
    V_cur = 0.0
    n_sub = max(1, int(np.ceil(abs(V_target) / v_step)))
    for k in range(1, n_sub + 1):
        V_next = V_target * k / n_sub
        r = newton_solve(x_, Na_, Nd_, V_next, phi=state[0], n=state[1],
                         p=state[2])
        assert r["converged"], f"V_a={V_next} 未收敛"
        state = (r["phi"], r["n"], r["p"])
        V_cur = V_next
    return r


def newton_history(x_, Na_, Nd_, V_a, state=None, max_iter=30):
    """返回 (iter 列表, 缩放残差范数历史)。逻辑与 solver.newton_solve 一致。"""
    Na_arr = np.where(x_ >= 0.0, 0.0, Na_)
    Nd_arr = np.where(x_ >= 0.0, Nd_, 0.0)
    bc = contact_values(V_a, Na_, Nd_)
    n_L, p_L, n_R, p_R = contact_carriers(bc, Na_, Nd_)
    hw_ = np.diff(x_)
    h_i_ = 0.5 * (hw_[:-1] + hw_[1:])

    if state is None:
        phi = initial_equilibrium_guess(x_, Na_, Nd_)
    else:
        phi = np.asarray(state[0], dtype=float).copy()
    phi[0], phi[-1] = bc["phi_left"], bc["phi_right"]
    if state is None:
        V_shift = V_a
        W_ramp = max(W_analytic(Na_, Nd_, 0.0), 4.0 * np.max(hw_))
        t = np.clip((0.5 * W_ramp - x_) / W_ramp, 0.0, 1.0)
        ramp = 0.5 * (1.0 - np.cos(np.pi * t))
        phi = phi + V_shift * ramp
        phi_n = np.full_like(x_, bc["phi_n_R"]) + V_shift * ramp
        phi_p = np.full_like(x_, bc["phi_p_R"]) + V_shift * ramp
    else:
        n_, p_ = np.asarray(state[1]), np.asarray(state[2])
        phi_n = phi - VT * np.log(np.maximum(n_, 1.0) / NI)
        phi_p = phi + VT * np.log(np.maximum(p_, 1.0) / NI)
    phi_n[0], phi_n[-1] = bc["phi_n_L"], bc["phi_n_R"]
    phi_p[0], phi_p[-1] = bc["phi_p_L"], bc["phi_p_R"]

    def assemble(ph, pn, pp):
        return _qfp_assemble(ph, pn, pp, x_, Na_arr, Nd_arr, hw_, h_i_,
                             n_L, p_L, n_R, p_R)

    F, Ad, Al, Au = assemble(phi, phi_n, phi_p)
    hist = [float(np.linalg.norm(F))]
    normF = hist[0]
    for it in range(1, max_iter + 1):
        rhs = -F.copy()
        Al[0] = 0.0
        Au[-1] = 0.0
        du = block_thomas(Ad, Al, Au, rhs)
        du_max = float(np.abs(du).max())
        sc = min(1.0, 0.5 / max(du_max, 1e-30)) if du_max > 0.5 else 1.0
        alpha, normF_try = 1.0, normF
        for _ in range(24):
            phi_t = phi.copy()
            pn_t = phi_n.copy()
            pp_t = phi_p.copy()
            phi_t[1:-1] += alpha * sc * du[:, 0]
            pn_t[1:-1] += alpha * sc * du[:, 1]
            pp_t[1:-1] += alpha * sc * du[:, 2]
            phi_t[0], phi_t[-1] = bc["phi_left"], bc["phi_right"]
            pn_t[0], pn_t[-1] = bc["phi_n_L"], bc["phi_n_R"]
            pp_t[0], pp_t[-1] = bc["phi_p_L"], bc["phi_p_R"]
            F_t, Ad_t, Al_t, Au_t = assemble(phi_t, pn_t, pp_t)
            normF_try = float(np.linalg.norm(F_t))
            if normF_try <= (1.0 - 2e-4 * alpha) * normF:
                break
            alpha *= 0.5
        phi, phi_n, phi_p = phi_t, pn_t, pp_t
        F, Ad, Al, Au = F_t, Ad_t, Al_t, Au_t
        dmove = alpha * sc * du_max
        normF_prev = normF
        normF = normF_try
        hist.append(normF)
        # 与 newton_solve 完全一致的收敛判据
        if dmove < 1e-9 and (normF < 1e-8 * hist[0] or
                             (normF >= normF_prev * (1.0 - 1e-2)
                              and normF < 1e-6)):
            break
    return np.arange(len(hist)), np.array(hist)


def dense(Ad, Al, Au):
    M = len(Ad)
    J = np.zeros((3 * M, 3 * M))
    for j in range(M):
        J[3 * j:3 * j + 3, 3 * j:3 * j + 3] = Ad[j]
        if j > 0:
            J[3 * j:3 * j + 3, 3 * (j - 1):3 * (j - 1) + 3] = Al[j]
        if j < M - 1:
            J[3 * j:3 * j + 3, 3 * (j + 1):3 * (j + 1) + 3] = Au[j]
    return J


def gummel_hist(x_, V_a, omega, max_outer=60):
    """复刻 solver.gummel_solve 外层循环，返回历史。"""
    bc = contact_values(V_a, Na, Nd)
    phi = initial_equilibrium_guess(x_, Na, Nd)
    phi = phi + (bc["phi_left"] - (-VT * np.log(Na / NI)))
    phi = phi.astype(float).copy()
    phi[0], phi[-1] = bc["phi_left"], bc["phi_right"]
    s = (x_ - x_[0]) / (x_[-1] - x_[0])
    phi_n = bc["phi_n_L"] * (1 - s) + bc["phi_n_R"] * s
    phi_p = bc["phi_p_L"] * (1 - s) + bc["phi_p_R"] * s
    n, p = carrier_from_phi(phi, phi_n, phi_p)
    n[0], p[0], n[-1], p[-1] = contact_carriers(bc, Na, Nd)
    pk = dict(tol=1e-11, max_iter=100, clamp=0.3, damped=True)
    dphi, nmax = [], []
    for it in range(1, max_outer + 1):
        n_s = solve_electron_continuity(phi, p, n, x_, bc,
                                        dict(Na=Na, Nd=Nd))
        n_new = (1 - omega) * n + omega * n_s
        p_s = solve_hole_continuity(phi, n_new, p, x_, bc,
                                    dict(Na=Na, Nd=Nd))
        p_new = (1 - omega) * p + omega * p_s
        phi_s, _, _ = solve_poisson_pn(phi, n_new, p_new,
                                       np.where(x_ >= 0, 0, Na),
                                       np.where(x_ >= 0, Nd, 0), x_, bc, **pk)
        phi_new = (1 - omega) * phi + omega * phi_s
        dphi.append(float(np.max(np.abs(phi_new - phi))))
        nmax.append(float(np.max(n_new)) / Nd)
        phi, n, p = phi_new, n_new, p_new
    return np.arange(1, max_outer + 1), np.array(dphi), np.array(nmax)


# ================================================================================
# 数据准备（单次计算，供全部插图与论文数据汇总使用）
# ================================================================================

def _prepare(verbose=True):
    """计算全部插图所需数据，返回 dict。计算顺序与原件脚本一致。"""
    d = {}
    d["Na"], d["Nd"], d["L"] = Na, Nd, L

    # ---- 网格 ----
    x = make_junction_grid(L)
    xc = 0.5 * (x[:-1] + x[1:])          # 界面中心
    hw = np.diff(x)
    h_i = 0.5 * (hw[:-1] + hw[1:])       # 控制体厚度（内部节点）
    d.update(x=x, xc=xc, hw=hw, h_i=h_i)
    if verbose:
        print(f"[网格] N = {len(x)}, 结区 {np.min(hw)*1e9:.2f} nm, "
              f"远场 {np.max(hw)*1e9:.2f} nm")

    # ---- 平衡态 ----
    if verbose:
        print("[平衡态] ...")
    res_eq = newton_solve(x, Na, Nd, 0.0, verbose=False)
    phi_eq, n_eq, p_eq = res_eq["phi"], res_eq["n"], res_eq["p"]
    from .run import verify_equilibrium
    ver = verify_equilibrium(x, Na, Nd, res_eq)
    V_bi = ver["V_bi_num"]
    d.update(res_eq=res_eq, phi_eq=phi_eq, n_eq=n_eq, p_eq=p_eq,
             ver=ver, V_bi=V_bi)
    if verbose:
        print(f"  V_bi 数值 {V_bi:.6f} V | 解析 {ver['V_bi_exact']:.6f} V | "
              f"max|np-n_i^2|/n_i^2 = {ver['np_err']:.2e} | "
              f"W = {ver['W_charge']*1e9:.1f} nm | iters {ver['iters']}")

    # 平衡态电场与电荷密度
    E_eq = -np.diff(phi_eq) / hw                     # 界面电场 (V/m)
    rho_eq = Q_E * (p_eq - n_eq + np.where(x >= 0, Nd, 0)
                    - np.where(x >= 0, 0, Na))       # 电荷密度 (C/m^3)
    d.update(E_eq=E_eq, rho_eq=rho_eq)

    # ---- 偏压特征点（参数延续：0.1 V 分步推进，与 sweep_bias 生产路径一致；
    # 单步 0.5 V 大跳变会进入慢收敛爬行——实测需 31 步才越过收敛判据，
    # 这正是本文采用参数延续策略的原因，见论文 3.4 节）----
    if verbose:
        print("[偏压点] ...")
    bias_points = {}
    for V in (+0.5, -3.0):
        # 每个分支从平衡态独立出发（与 sweep_bias 生产路径一致：
        # 从正偏解穿越 0 V 继续走反偏会因残留过量载流子卡死）
        r = solve_steps(x, Na, Nd, V, (phi_eq, n_eq, p_eq))
        bias_points[V] = r
        if verbose:
            print(f"  V_a = {V:+.2f} V : iters {r['iters']}, "
                  f"J = {current_density(r['phi'], r['n'], r['p'], x)[0]:+.4e} A/m^2")
    d["bias_points"] = bias_points
    r5 = bias_points[+0.5]
    d["r5"] = r5

    # ---- 偏压扫描（J-V 与 W-V）----
    if verbose:
        print("[偏压扫描] ...")
    v_fwd = np.round(np.arange(0.05, 0.85, 0.05), 4)
    v_rev = np.round(-np.arange(0.1, 5.05, 0.1), 4)
    records = sweep_bias(x, Na, Nd, v_fwd, v_rev, verbose=False)
    Va_arr = np.array([r["V_a"] for r in records])
    J_arr = np.array([r["J"] for r in records])
    W_arr = np.array([r["W_charge"] for r in records])
    iters_arr = np.array([r["iters"] for r in records])
    n_bad = sum(1 for r in records if not r["converged"])
    d.update(records=records, Va_arr=Va_arr, J_arr=J_arr, W_arr=W_arr,
             iters_arr=iters_arr, n_bad=n_bad)
    if verbose:
        print(f"  {len(records)} 点, 未收敛 {n_bad}, "
              f"平均 Newton {iters_arr.mean():.1f} 步/点, "
              f"J 范围 [{J_arr.min():.3e}, {J_arr.max():.3e}] A/m^2")

    Js_exact = Js_analytic(Na, Nd)
    d["Js_exact"] = Js_exact
    # 三段非重叠区间拟合（论文表 3；与 run.run() 输出一致）
    seg_fits = []
    for lo, hi in [(0.10, 0.30), (0.35, 0.65), (0.70, 0.80)]:
        m = (Va_arr >= lo - 1e-9) & (Va_arr <= hi + 1e-9)
        c = np.polyfit(Va_arr[m], np.log(np.abs(J_arr[m])), 1)
        seg_fits.append((lo, hi, c[0], 1.0 / (c[0] * VT)))
        if verbose:
            print(f"  分段拟合 {lo:.2f}-{hi:.2f} V: 斜率 {c[0]:.2f} /V, "
                  f"n_id = {seg_fits[-1][3]:.3f}")
    d["seg_fits"] = seg_fits

    # ---- 恒等式验证量（图 5 与论文汇总用）：Jn(L)-Jn(0) = q*int(U)dx ----
    phi5, n5, p5 = r5["phi"], r5["n"], r5["p"]
    U5, _, _ = srh_rate(n5, p5)
    intU5 = float(np.sum(h_i * U5[1:-1]))          # 内部节点控制体积分
    J5, Jt5, Jn5, Jp5 = current_density(phi5, n5, p5, x)
    id5 = Jn5[-1] - Jn5[0]
    d.update(intU5=intU5, id5=id5)

    # ---- 图 7：Newton 收敛历史 + 雅可比条件数 ----
    # 7a：自建历史记录循环（复刻 newton_solve 主循环，仅记录 ||F||）
    it_eq, h_eq = newton_history(x, Na, Nd, 0.0)
    # 偏压延续单步历史（0.40 -> 0.45，0.05 V 步长）
    r40 = solve_steps(x, Na, Nd, 0.40, (phi_eq, n_eq, p_eq))
    it_b, h_b = newton_history(x, Na, Nd, 0.45,
                               state=(r40["phi"], r40["n"], r40["p"]))
    d.update(it_eq=it_eq, h_eq=h_eq, it_b=it_b, h_b=h_b)
    if verbose:
        print(f"  [Newton] 平衡态 {len(it_eq)-1} 步, ||F|| {h_eq[0]:.2e} -> "
              f"{h_eq[-1]:.2e}; 延续步 {len(it_b)-1} 步")

    # 7b：条件数对比（粗网格，稠密矩阵）
    xc2 = make_junction_grid(2.0e-6, d0=20e-9, transition=0.15e-6,
                             h_max=40e-9)
    Na2 = np.where(xc2 >= 0.0, 0.0, Na)
    Nd2 = np.where(xc2 >= 0.0, Nd, 0.0)
    hw2 = np.diff(xc2)
    h_i2 = 0.5 * (hw2[:-1] + hw2[1:])
    phi_g = initial_equilibrium_guess(xc2, Na, Nd)
    n_g, p_g = carrier_from_phi(phi_g, np.zeros_like(phi_g),
                                np.zeros_like(phi_g))
    F_g, Ad_g, Al_g, Au_g = _newton_system(phi_g, n_g, p_g, xc2, Na2, Nd2,
                                           hw2, h_i2)
    Ad_t = qfp_columns(Ad_g, n_g[1:-1] / VT, p_g[1:-1] / VT)
    Al_t = qfp_columns(Al_g, n_g[:-2] / VT, p_g[:-2] / VT)
    Au_t = qfp_columns(Au_g, n_g[2:] / VT, p_g[2:] / VT)

    J_raw = dense(Ad_g, Al_g, Au_g)
    J_qfp = dense(Ad_t, Al_t, Au_t)
    c_raw = np.linalg.cond(J_raw)
    c_qfp = np.linalg.cond(J_qfp)

    # 行缩放（与 _qfp_assemble 第 4 步一致）：求解器实际使用的矩阵
    inv_hsum = 1.0 / (1.0 / hw2[:-1] + 1.0 / hw2[1:])
    n_eff = np.maximum(n_g[1:-1], 1e8)
    p_eff = np.maximum(p_g[1:-1], 1e8)
    scale = np.stack([inv_hsum / EPS_SI,
                      inv_hsum * VT / (Q_E * D_N * n_eff),
                      inv_hsum * VT / (Q_E * D_P * p_eff)], axis=1)
    J_sc = dense(Ad_t * scale[:, :, None], Al_t * scale[:, :, None],
                 Au_t * scale[:, :, None])
    c_sc = np.linalg.cond(J_sc)
    d.update(xc2=xc2, c_raw=c_raw, c_qfp=c_qfp, c_sc=c_sc)
    if verbose:
        print(f"  [条件数] raw {c_raw:.3e} -> QFP {c_qfp:.3e} -> "
              f"QFP+scaling {c_sc:.3e}")

    # ---- 图 8：Gummel 解耦迭代对照 ----
    if verbose:
        print("[Gummel 对照] ...")
    it_g1, dp_g1, nm_g1 = gummel_hist(x, 0.0, 1.0, max_outer=60)
    it_g2, dp_g2, nm_g2 = gummel_hist(x, 0.0, 0.5, max_outer=60)
    d.update(it_g1=it_g1, dp_g1=dp_g1, nm_g1=nm_g1,
             it_g2=it_g2, dp_g2=dp_g2, nm_g2=nm_g2)
    if verbose:
        print(f"  omega=1.0: dphi 末值 {dp_g1[-1]:.2e} V, max(n)/N_d 范围 "
              f"[{nm_g1.min():.1f}, {nm_g1.max():.1f}]")
        print(f"  omega=0.5: dphi 末值 {dp_g2[-1]:.2e} V, max(n)/N_d 范围 "
              f"[{nm_g2.min():.1f}, {nm_g2.max():.1f}]")
    return d


# ================================================================================
# 插图（论文编号 S2-F1 ~ S2-F9）
# ================================================================================

def fig0_structure(d, outdir=None):
    """图 0（S2-F1）：器件结构示意。"""
    outdir = _outdir(outdir)
    ver = d["ver"]
    W0 = ver["W_charge"]
    fig, ax = plt.subplots(figsize=(7.2, 3.0))
    ax.set_xlim(-2.6, 2.6)
    ax.set_ylim(-1.15, 1.55)
    ax.axis("off")

    # 半导体体区
    ax.add_patch(Rectangle((-2.5, -0.55), 5.0, 1.1, facecolor="#f2f2f2",
                           edgecolor="k", lw=1.2, zorder=1))
    # 耗尽区（平衡态宽度）
    ax.add_patch(Rectangle((-W0 / 2 * 1e6, -0.55), W0 * 1e6, 1.1,
                           facecolor="#cfe3f7", edgecolor="none", zorder=2))
    ax.plot([0, 0], [-0.55, 0.55], "k--", lw=1.0, zorder=3)
    ax.annotate("", xy=(-W0 / 2 * 1e6, 0.78), xytext=(W0 / 2 * 1e6, 0.78),
                arrowprops=dict(arrowstyle="<->", color="k", lw=1.0))
    ax.text(0, 0.86, f"metallurgical junction\n$W_0$ = {W0*1e9:.0f} nm",
            ha="center", va="bottom", fontsize=9.5)

    ax.text(-1.9, 0.28, "p-type", fontsize=11, ha="center", weight="bold")
    ax.text(1.9, 0.28, "n-type", fontsize=11, ha="center", weight="bold")
    ax.text(-1.9, 0.02, "$N_A=10^{16}$ cm$^{-3}$", fontsize=10, ha="center")
    ax.text(1.9, 0.02, "$N_D=10^{16}$ cm$^{-3}$", fontsize=10, ha="center")
    ax.text(-0.9, -0.32, "depletion (negative)", fontsize=8.5, ha="center",
            color="#205080")
    ax.text(0.9, -0.32, "depletion (positive)", fontsize=8.5, ha="center",
            color="#205080")

    # 欧姆接触
    for xs, lab in ((-2.5, "p$^+$ contact\n(anode)"), (2.5, "n$^+$ contact\n(cathode)")):
        ax.add_patch(Rectangle((xs - 0.14 if xs < 0 else xs, -0.55), 0.14, 1.1,
                               facecolor="#d9d9d9", edgecolor="k", lw=1.0, zorder=4))
        ax.text(xs + (-0.28 if xs < 0 else 0.28), -0.78, lab,
                fontsize=9, ha="center", va="top")

    # 偏压箭头（置于最上方，避免与结标注重叠——旧版曾因 y=1.05 与
    # "metallurgical junction" 文本重叠，白底框盖住 "gical"）
    ax.add_patch(FancyArrowPatch((-3.0, 1.32), (-2.62, 1.32),
                                 arrowstyle="-|>", mutation_scale=12, color="k"))
    ax.add_patch(FancyArrowPatch((2.62, 1.32), (3.0, 1.32),
                                 arrowstyle="-|>", mutation_scale=12, color="k"))
    ax.text(0, 1.32, "$V_a$", fontsize=12, ha="center", va="center",
            bbox=dict(fc="white", ec="k", boxstyle="round,pad=0.25"))
    ax.text(0, -1.02, "$L$ = 5 $\\mu$m", fontsize=9, ha="center", va="top")
    ax.annotate("", xy=(-2.5, -0.98), xytext=(2.5, -0.98),
                arrowprops=dict(arrowstyle="<->", color="k", lw=0.8))
    fig.savefig(os.path.join(outdir, "fig0_structure.png"), bbox_inches="tight")
    plt.close(fig)


def fig1_equilibrium(d, outdir=None):
    """图 1（S2-F2）：平衡态能带图 + 电势/电荷。"""
    outdir = _outdir(outdir)
    x, xm = d["x"], d["x"] * 1e6
    phi_eq = d["phi_eq"]
    rho_eq = d["rho_eq"]
    ver, V_bi = d["ver"], d["V_bi"]
    W0 = ver["W_charge"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.0))

    # 参考：平衡态单费米能级 E_F = 0，E_i = -phi，E_c = E_i + E_g/2
    Ei = -phi_eq
    Ec = Ei + EG / 2
    Ev = Ei - EG / 2
    sel = np.abs(xm) <= 1.6
    ax1.plot(xm[sel], Ec[sel], color=C_N, lw=1.6, label="$E_c$")
    ax1.plot(xm[sel], Ei[sel], color="k", lw=1.2, ls="--", label="$E_i$")
    ax1.plot(xm[sel], Ev[sel], color=C_P, lw=1.6, label="$E_v$")
    ax1.axhline(0, color="#9467bd", lw=1.2, label="$E_F$")
    ax1.annotate("", xy=(-1.45, Ec[np.argmin(np.abs(xm + 1.45))]),
                 xytext=(-1.45, 0.0),
                 arrowprops=dict(arrowstyle="<->", color="k", lw=0.9))
    ax1.text(-1.38, 0.45, "$E_g/2+\\psi_F$\n0.917 eV", fontsize=8.5, va="center")
    ax1.annotate("", xy=(0.0, 0.917), xytext=(0.0, 0.203),
                 arrowprops=dict(arrowstyle="<->", color="k", lw=0.9))
    ax1.text(0.07, 0.56, "$qV_{bi}$\n0.714 eV", fontsize=8.5, va="center")
    ax1.set_xlabel("position $x$ ($\\mu$m)")
    ax1.set_ylabel("energy (eV)")
    ax1.set_title("(a) equilibrium band diagram ($V_a=0$)")
    ax1.set_ylim(-0.75, 1.05)
    ax1.legend(loc="center left", framealpha=0.9)

    # 右：电势与电荷密度
    ax2.plot(xm, phi_eq, color=C_PHI, lw=1.6, label="numerical $\\varphi(x)$")
    xa = np.linspace(-2.5, 2.5, 400)
    phi_ana = np.where(
        xa < -W0 / 2, -V_bi / 2,
        np.where(xa > W0 / 2, +V_bi / 2,
                 V_bi * xa / W0))   # 耗尽近似：三角势
    ax2.plot(xa, phi_ana, color=C_ANA, lw=1.2, ls="--",
             label="depletion approximation")
    ax2.axhline(0, color="k", lw=0.6)
    ax2.set_xlabel("position $x$ ($\\mu$m)")
    ax2.set_ylabel("electrostatic potential $\\varphi$ (V)", color=C_PHI)
    ax2.set_title("(b) potential profile at equilibrium")
    ax2.set_xlim(-2.5, 2.5)
    ax2.legend(loc="upper left")

    ax2b = ax2.twinx()
    ax2b.plot(xm, rho_eq * 1e-3, color=C_E, lw=1.2, alpha=0.85)
    ax2b.set_ylabel("charge density (10$^3$ C/m$^3$)", color=C_E)
    ax2b.tick_params(axis="y", labelcolor=C_E)
    ax2b.grid(False)
    ax2b.set_ylim(-1.5 * Q_E * Na * 1e-3 * 1.02, 1.5 * Q_E * Na * 1e-3 * 1.02)
    fig.savefig(os.path.join(outdir, "fig1_equilibrium.png"), bbox_inches="tight")
    plt.close(fig)


def fig2_bias_profiles(d, outdir=None):
    """图 2（S2-F3）：正偏 / 反偏的载流子浓度与电场。"""
    outdir = _outdir(outdir)
    x, xc, hw = d["x"], d["xc"], d["hw"]
    bias_points = d["bias_points"]
    xm = x * 1e6

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.2))
    for row, V in enumerate((+0.5, -3.0)):
        r = bias_points[V]
        n, p = r["n"] * 1e-6, r["p"] * 1e-6        # cm^-3
        E = -np.diff(r["phi"]) / hw * 1e-2         # V/cm（取界面）
        axL, axR = axes[row]

        axL.semilogy(xm, n, color=C_N, lw=1.5, label="$n$")
        axL.semilogy(xm, p, color=C_P, lw=1.5, label="$p$")
        axL.axhline(Na * 1e-6, color=C_P, lw=0.9, ls=":", label="$N_A=N_D$")
        axL.axhline(NI * 1e-6, color="k", lw=0.9, ls=":", label="$n_i$")
        axL.set_xlim(-2.5, 2.5)
        axL.set_ylim(1e2, 1e19)
        axL.set_ylabel("carrier density (cm$^{-3}$)")
        axL.set_title(f"({'ab'[row]}) carriers, $V_a={V:+.1f}$ V")
        axL.legend(loc="lower center", ncol=2)

        axR.plot(xc * 1e6, E, color=C_E, lw=1.5)
        axR.axhline(0, color="k", lw=0.6)
        axR.set_xlim(-2.5, 2.5)
        axR.set_xlabel("position $x$ ($\\mu$m)")
        axR.set_ylabel("electric field (V/cm)")
        axR.set_title(f"({'cd'[row]}) electric field, $V_a={V:+.1f}$ V")
        Emax = np.abs(E).max()
        axR.set_ylim(-1.25 * Emax, 1.25 * Emax)
    fig.savefig(os.path.join(outdir, "fig2_bias_profiles.png"), bbox_inches="tight")
    plt.close(fig)


def fig3_jv_curve(d, outdir=None):
    """图 3（S2-F4）：J-V 特性。"""
    outdir = _outdir(outdir)
    Va_arr, J_arr = d["Va_arr"], d["J_arr"]
    Js_exact, seg_fits = d["Js_exact"], d["seg_fits"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.2))

    fwd = Va_arr > 0
    rev = Va_arr < 0
    ax1.semilogy(Va_arr[fwd], J_arr[fwd], "o", ms=4.5, color=C_N,
                 label="full-Newton (this work)")
    Vv = np.linspace(0.0, 0.85, 200)
    ax1.semilogy(Vv, Js_exact * (np.exp(Vv / VT) - 1.0), "--",
                 color=C_ANA, lw=1.3, label="Shockley $J_s[e^{V/V_T}-1]$")
    ax1.semilogy(-Vv[Vv > 0.05], -Js_exact * (np.exp(-Vv[Vv > 0.05] / VT) - 1.0),
                 "--", color=C_ANA, lw=1.3)
    ax1.semilogy(Va_arr[rev], -J_arr[rev], "s", ms=4, color=C_P,
                 label="reverse (generation)")
    ax1.axhline(Js_exact, color="#9467bd", lw=1.0, ls=":")
    ax1.text(-5.1, Js_exact * 1.5, "diffusion saturation $J_s$", fontsize=8.5,
             color="#9467bd")
    ax1.set_xlabel("applied bias $V_a$ (V)")
    ax1.set_ylabel("$|J|$ (A/m$^2$)")
    ax1.set_title("(a) semilog $J$–$V$")
    ax1.legend(loc="lower right")
    ax1.set_xlim(-5.2, 0.9)

    # 右：线性坐标正偏 + 分段理想因子说明
    ax2.plot(Va_arr[fwd], J_arr[fwd] * 1e-4, "o-", ms=4, color=C_N,
             label="numerical")
    ax2.plot(Vv, Js_exact * (np.exp(Vv / VT) - 1.0) * 1e-4, "--",
             color=C_ANA, lw=1.3, label="Shockley (ideal)")
    ax2.set_xlabel("applied bias $V_a$ (V)")
    ax2.set_ylabel("$J$ (A/cm$^2$)")
    ax2.set_title("(b) forward $J$–$V$ (linear scale)")
    ax2.legend(loc="upper left")
    ax2.text(0.04, 2.0e5,
             f"$J_s$ = {Js_exact:.3e} A/m$^2$\n"
             f"segmented $n_{{id}}$:\n"
             + "\n".join(f"  {lo:.2f}–{hi:.2f} V: {nid:.2f}"
                         for lo, hi, s, nid in seg_fits),
             fontsize=9, va="top",
             bbox=dict(fc="white", ec="#888", boxstyle="round,pad=0.3"))
    fig.savefig(os.path.join(outdir, "fig3_jv_curve.png"), bbox_inches="tight")
    plt.close(fig)


def fig4_quasi_fermi(d, outdir=None):
    """图 4（S2-F5）：正偏偏费米势分裂（能带图 + 准费米势）。"""
    outdir = _outdir(outdir)
    x = d["x"]
    xm = x * 1e6
    r5 = d["r5"]
    V5 = 0.5
    phi5, n5, p5 = r5["phi"], r5["n"], r5["p"]
    phi_n5 = phi5 - VT * np.log(np.maximum(n5, 1.0) / NI)
    phi_p5 = phi5 + VT * np.log(np.maximum(p5, 1.0) / NI)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.2))
    Ei5 = -phi5
    ax1.plot(xm, Ei5 + EG / 2, color=C_N, lw=1.5, label="$E_c$")
    ax1.plot(xm, Ei5, color="k", lw=1.1, ls="--", label="$E_i$")
    ax1.plot(xm, Ei5 - EG / 2, color=C_P, lw=1.5, label="$E_v$")
    ax1.plot(xm, -phi_n5, color=C_N, lw=1.8, ls="-.", label="$E_{Fn}$")
    ax1.plot(xm, -phi_p5, color="#8c564b", lw=1.8, ls=":", label="$E_{Fp}$")
    ax1.axhline(0, color="#9467bd", lw=1.0, ls="-", alpha=0.7)
    ax1.text(2.35, 0.03, "$E_F$ (equil.)", fontsize=8.5, color="#9467bd",
             ha="right")
    ax1.set_xlim(-2.5, 2.5)
    ax1.set_ylim(-0.9, 1.0)
    ax1.set_xlabel("position $x$ ($\\mu$m)")
    ax1.set_ylabel("energy (eV)")
    ax1.set_title(f"(a) quasi-Fermi level splitting, $V_a=+{V5}$ V")
    ax1.legend(loc="center left", fontsize=8.5, framealpha=0.9)

    # 右：少子注入浓度
    n_inj = n5 * 1e-6
    p_inj = p5 * 1e-6
    ax2.semilogy(xm, n_inj, color=C_N, lw=1.5, label="$n$")
    ax2.semilogy(xm, p_inj, color=C_P, lw=1.5, label="$p$")
    ax2.axhline(Na * 1e-6, color=C_P, lw=0.9, ls=":", label="$N_A=N_D$")
    ax2.axhline(NI * 1e-6, color="k", lw=0.9, ls=":", label="$n_i$")
    # 注入少子的指数衰减（p 侧电子、n 侧空穴）与解析扩散长度对照
    n_edge = NI**2 / Na * np.exp(V5 / VT) * 1e-6     # cm^-3
    p_edge = n_edge
    Ln = np.sqrt(D_N * TAU_N)
    Lp = np.sqrt(D_P * TAU_P)
    mask_p = xm < -0.2
    mask_n = xm > 0.2
    ax2.semilogy(xm[mask_p], n_edge * np.exp(-(np.abs(xm[mask_p]) - 0.2) /
                                              (Ln * 1e6)),
                 "--", color=C_ANA, lw=1.2,
                 label="$n_{p0}e^{V/V_T}e^{-x/L_n}$")
    ax2.semilogy(xm[mask_n], p_edge * np.exp(-(xm[mask_n] - 0.2) /
                                              (Lp * 1e6)),
                 "--", color=C_ANA, lw=1.2,
                 label="$p_{n0}e^{V/V_T}e^{-x/L_p}$")
    ax2.set_xlim(-2.5, 2.5)
    ax2.set_ylim(1e2, 1e19)
    ax2.set_xlabel("position $x$ ($\\mu$m)")
    ax2.set_ylabel("carrier density (cm$^{-3}$)")
    ax2.set_title(f"(b) minority-carrier injection, $V_a=+{V5}$ V")
    ax2.legend(loc="lower center", ncol=2, fontsize=8.5)
    fig.savefig(os.path.join(outdir, "fig4_quasi_fermi.png"), bbox_inches="tight")
    plt.close(fig)


def fig5_recombination(d, outdir=None):
    """图 5（S2-F6）：SRH 复合/产生率与电流守恒恒等式。"""
    outdir = _outdir(outdir)
    x, xc, h_i = d["x"], d["xc"], d["h_i"]
    r5 = d["r5"]
    xm = x * 1e6
    phi5, n5, p5 = r5["phi"], r5["n"], r5["p"]

    U5, _, _ = srh_rate(n5, p5)
    U3, _, _ = srh_rate(d["bias_points"][-3.0]["n"],
                        d["bias_points"][-3.0]["p"])
    # U 定义在节点上；绘图用节点值
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.2))
    ax1.semilogy(xm[50:-50], np.abs(U5[50:-50]), color=C_N, lw=1.5,
                 label="$|U|$, $V_a=+0.5$ V (recombination)")
    ax1.semilogy(xm[50:-50], np.abs(U3[50:-50]), color=C_P, lw=1.5,
                 label="$|U|$, $V_a=-3$ V (generation)")
    ax1.axvline(0, color="k", lw=0.6, ls=":")
    ax1.set_xlim(-2.5, 2.5)
    ax1.set_xlabel("position $x$ ($\\mu$m)")
    ax1.set_ylabel("$|U|$ (m$^{-3}$s$^{-1}$)")
    ax1.set_title("(a) SRH net recombination rate")
    ax1.legend(loc="upper center")

    # 恒等式验证：Jn[-1] - Jn[0] = q * ∫U dx
    intU5 = d["intU5"]
    id5 = d["id5"]
    J5, Jt5, Jn5, Jp5 = current_density(phi5, n5, p5, x)
    ax2.plot(xc * 1e6, Jt5 * 1e-4, color="k", lw=1.8, label="$J_t=J_n+J_p$")
    ax2.plot(xc * 1e6, Jn5 * 1e-4, color=C_N, lw=1.3, label="$J_n$")
    ax2.plot(xc * 1e6, Jp5 * 1e-4, color=C_P, lw=1.3, label="$J_p$")
    ax2.set_xlim(-2.5, 2.5)
    ax2.set_xlabel("position $x$ ($\\mu$m)")
    ax2.set_ylabel("current density (A/cm$^2$)")
    ax2.set_title(f"(b) current conservation at $V_a=+0.5$ V")
    ax2.legend(loc="center right", fontsize=8.5)
    ax2.text(-2.35, ax2.get_ylim()[0] + 0.06 * (ax2.get_ylim()[1] - ax2.get_ylim()[0]),
             f"$J_n(L)-J_n(0)$ = {id5*1e-4:.4e} A/cm$^2$\n"
             f"$q\\int U\\,dx$ = {Q_E*intU5*1e-4:.4e} A/cm$^2$",
             fontsize=8.5, va="bottom",
             bbox=dict(fc="white", ec="#888", boxstyle="round,pad=0.3"))
    print(f"  [恒等式] Jn(L)-Jn(0) = {id5:.6e}, q*int(U) = {Q_E*intU5:.6e}, "
          f"ratio = {id5/(Q_E*intU5):.4f}")
    fig.savefig(os.path.join(outdir, "fig5_recombination.png"), bbox_inches="tight")
    plt.close(fig)


def fig6_depletion(d, outdir=None):
    """图 6（S2-F7）：耗尽层宽度 vs 偏压。"""
    outdir = _outdir(outdir)
    Va_arr, W_arr = d["Va_arr"], d["W_arr"]
    ver = d["ver"]

    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    W_ana = np.array([W_analytic(Na, Nd, float(v)) for v in Va_arr])
    m_ana = Va_arr < V_bi_analytic(Na, Nd) - 0.02
    ax.plot(Va_arr, W_arr * 1e9, "o", ms=4, color=C_N, label="numerical $W$")
    ax.plot(Va_arr[m_ana], W_ana[m_ana] * 1e9, "--", color=C_ANA, lw=1.4,
            label="depletion approximation")
    ax.axhline(ver["W_charge"] * 1e9, color=C_P, lw=0.9, ls=":")
    ax.text(-4.9, ver["W_charge"] * 1e9 + 12, f"$W_0$ = {ver['W_charge']*1e9:.0f} nm",
            fontsize=9, color=C_P)
    ax.set_xlabel("applied bias $V_a$ (V)")
    ax.set_ylabel("depletion width $W$ (nm)")
    ax.set_title("depletion width vs. applied bias")
    ax.legend(loc="upper left")
    ax.set_xlim(-5.2, 0.9)
    fig.savefig(os.path.join(outdir, "fig6_depletion.png"), bbox_inches="tight")
    plt.close(fig)


def fig7_newton(d, outdir=None):
    """图 7（S2-F8）：Newton 收敛历史 + 雅可比条件数。"""
    outdir = _outdir(outdir)
    it_eq, h_eq = d["it_eq"], d["h_eq"]
    it_b, h_b = d["it_b"], d["h_b"]
    c_raw, c_qfp, c_sc = d["c_raw"], d["c_qfp"], d["c_sc"]
    xc2 = d["xc2"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.2))
    ax1.semilogy(it_eq, h_eq / h_eq[0], "o-", ms=5, color=C_N,
                 label="equilibrium, cold start")
    ax1.semilogy(it_b, h_b / h_b[0], "s-", ms=5, color=C_P,
                 label="continuation step 0.40$\\to$0.45 V")
    ax1.axhline(1e-8, color=C_ANA, ls="--", lw=1.0)
    ax1.text(0.6, 1.6e-8, "convergence floor $10^{-8}$", fontsize=8.5,
             color=C_ANA)
    ax1.set_xlabel("Newton iteration")
    ax1.set_ylabel(r"$\|F\| / \|F_0\|$ (scaled)")
    ax1.set_title("(a) Newton convergence history")
    ax1.set_ylim(1e-13, 3)
    ax1.legend(loc="lower left")

    bars = ax2.bar(["$(\\varphi,n,p)$\nraw",
                    "$(\\varphi,\\varphi_n,\\varphi_p)$\nQFP",
                    "QFP\n+ row scaling"],
                   [c_raw, c_qfp, c_sc], width=0.45,
                   color=[C_P, C_N, "#2ca02c"], edgecolor="k", lw=0.8)
    ax2.set_yscale("log")
    ax2.set_ylabel("condition number $\\kappa_2(J)$")
    ax2.set_title(f"(b) Jacobian condition number (coarse grid, $N$={len(xc2)})")
    for b, v in zip(bars, [c_raw, c_qfp, c_sc]):
        ax2.text(b.get_x() + b.get_width() / 2, v * 1.6, f"{v:.1e}",
                 ha="center", fontsize=9.5)
    ax2.set_ylim(1e0, 1e40)
    ax2.grid(axis="x")
    fig.savefig(os.path.join(outdir, "fig7_newton.png"), bbox_inches="tight")
    plt.close(fig)


def fig8_gummel(d, outdir=None):
    """图 8（S2-F9）：Gummel 解耦迭代对照（论文 3.3/4.8 节引用）。"""
    outdir = _outdir(outdir)
    it_g1, dp_g1, nm_g1 = d["it_g1"], d["dp_g1"], d["nm_g1"]
    it_g2, dp_g2, nm_g2 = d["it_g2"], d["dp_g2"], d["nm_g2"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.0))
    ax1.semilogy(it_g1, dp_g1, "o-", ms=4, color=C_P,
                 label="$\\omega$ = 1 (raw Gummel)")
    ax1.semilogy(it_g2, dp_g2, "s-", ms=4, color=C_N,
                 label="$\\omega$ = 0.5 (under-relaxed)")
    ax1.axhline(1e-9, color=C_ANA, ls="--", lw=1.0)
    ax1.text(1.2, 1.6e-9, "convergence tol $10^{-9}$ V", fontsize=8.5,
             color=C_ANA)
    ax1.set_xlabel("outer (Gummel) iteration")
    ax1.set_ylabel(r"$|\Delta\varphi|_\infty$ (V)")
    ax1.set_title("(a) Gummel outer-iteration stagnation ($V_a=0$)")
    ax1.set_ylim(1e-10, 1)
    ax1.legend(loc="lower left")

    ax2.semilogy(it_g1, nm_g1, "o-", ms=4, color=C_P,
                 label="$\\omega$ = 1")
    nm_g2_plot = np.clip(nm_g2, 0.5, 1e3)   # 爆冲值裁剪到坐标范围，另加标注
    ax2.semilogy(it_g2, nm_g2_plot, "s-", ms=4, color=C_N,
                 label="$\\omega$ = 0.5")
    i_blow = int(np.argmax(nm_g2 > 1e3))
    if nm_g2.max() > 1e3:
        ax2.annotate("carrier blow-up\n(clipped, max = %.0e)" % nm_g2.max(),
                     xy=(it_g2[i_blow], 1e3), xytext=(it_g2[i_blow] - 22, 3e2),
                     fontsize=8.5,
                     arrowprops=dict(arrowstyle="->", lw=0.9))
    ax2.axhline(1.0, color="k", lw=0.9, ls=":")
    ax2.text(1.2, 1.25, "$n = N_D$ (equilibrium)", fontsize=8.5)
    ax2.set_xlabel("outer (Gummel) iteration")
    ax2.set_ylabel(r"max$(n)/N_D$")
    ax2.set_title("(b) carrier oscillation during decoupled iteration")
    ax2.set_ylim(0.5, 3e3)
    ax2.legend(loc="upper right")
    fig.savefig(os.path.join(outdir, "fig8_gummel.png"), bbox_inches="tight")
    plt.close(fig)


# ================================================================================
# 论文正文所需的汇总数据
# ================================================================================

def print_summary(d):
    """打印论文数据汇总（与 run.run() 的 verbose 输出互补）。"""
    ver = d["ver"]
    Js_exact = d["Js_exact"]
    Va_arr, J_arr, W_arr = d["Va_arr"], d["J_arr"], d["W_arr"]
    records = d["records"]
    iters_arr = d["iters_arr"]
    id5, intU5 = d["id5"], d["intU5"]
    c_raw, c_qfp, c_sc = d["c_raw"], d["c_qfp"], d["c_sc"]

    print("\n" + "=" * 60)
    print("论文数据汇总")
    print("=" * 60)
    print(f"V_bi 数值 = {ver['V_bi_num']:.10f} V, 解析 = {ver['V_bi_exact']:.10f} V, "
          f"相对偏差 = {abs(ver['V_bi_num']-ver['V_bi_exact'])/ver['V_bi_exact']:.2e}")
    print(f"平衡态 max|np-n_i^2|/n_i^2 = {ver['np_err']:.3e}")
    print(f"平衡态 J = {ver['J_eq']:.3e} A/m^2 (尺度 {ver['J_scale']:.3e})")
    print(f"W 电荷法 = {ver['W_charge']*1e9:.2f} nm, "
          f"解析 = {ver['W_exact']*1e9:.2f} nm, "
          f"偏差 = {abs(ver['W_charge']-ver['W_exact'])/ver['W_exact']*100:.2f}%")
    print(f"平衡态 Newton 迭代 {ver['iters']} 步")
    print(f"J_s = {Js_exact:.4e} A/m^2 = {Js_exact*1e-4:.4e} A/cm^2")
    print("分段理想因子（不重叠区间）:")
    for lo, hi, slope, nid in d["seg_fits"]:
        print(f"  {lo:.2f}-{hi:.2f} V: 斜率 {slope:.2f} /V, n_id = {nid:.3f}")
    # 理想 Shockley 失效倍数
    for V in (0.70, 0.75, 0.80):
        i = int(np.argmin(np.abs(Va_arr - V)))
        ideal = Js_exact * np.exp(Va_arr[i] / VT)
        print(f"  V={Va_arr[i]:.2f}: J={J_arr[i]:.4e}, 理想={ideal:.4e}, "
              f"高估 {ideal/J_arr[i]:.0f} 倍")
    rec_rev = [r for r in records if r["V_a"] <= -3.0]
    Jr = np.mean([abs(r["J"]) for r in rec_rev])
    print(f"反偏 (-3~-5V) 平均 |J| = {Jr:.4e} A/m^2 = {Jr*1e-4:.4e} A/cm^2")
    print(f"  其中扩散分量 J_s 占比 {Js_exact/Jr:.2e}")
    print(f"  W(-3V) = {[r['W_charge']*1e9 for r in records if r['V_a']==-3.0][0]:.1f} nm, "
          f"W(-5V) = {W_arr[-1]*1e9:.1f} nm")
    J_gen = Q_E * NI * W_arr[-1] / (2 * TAU_N)
    print(f"  q n_i W/(2 tau) 估计 = {J_gen:.4e} A/m^2, 数值/估计 = {Jr/J_gen:.3f}")
    print(f"恒等式 V=+0.5V: Jn(L)-Jn(0) = {id5:.6e}, q*int U = {Q_E*intU5:.6e}, "
          f"ratio = {id5/(Q_E*intU5):.4f}")
    print(f"条件数: raw {c_raw:.3e} -> QFP {c_qfp:.3e} -> QFP+scaling {c_sc:.3e} "
          f"(raw/scaled 比值 {c_raw/c_sc:.1e})")
    print(f"扫描: {len(records)} 点, 未收敛 {d['n_bad']}, 平均 {iters_arr.mean():.1f} 步/点, "
          f"最多 {iters_arr.max()} 步")
    print(f"正偏 J(0.8V) = {J_arr[np.argmax(Va_arr)]:.4e} A/m^2")
    print(f"扩散长度 L_n = {np.sqrt(D_N*TAU_N)*1e6:.3f} um, "
          f"L_p = {np.sqrt(D_P*TAU_P)*1e6:.3f} um")
    print(f"德拜长度 L_D = {np.sqrt(EPS_SI*VT/(Q_E*Na))*1e9:.1f} nm")


# ================================================================================
# 编排入口
# ================================================================================

def make_figures(data=None, outdir=None, verbose=True):
    """生成第二阶段全部插图（fig0-fig8）到 results/stage2/。

    data 缺省时现场调用 _prepare()（约 1~2 分钟，主要是 66 点偏压扫描）。
    """
    outdir = _outdir(outdir)
    if data is None:
        data = _prepare(verbose=verbose)
    fig0_structure(data, outdir)
    fig1_equilibrium(data, outdir)
    fig2_bias_profiles(data, outdir)
    fig3_jv_curve(data, outdir)
    fig4_quasi_fermi(data, outdir)
    fig5_recombination(data, outdir)
    fig6_depletion(data, outdir)
    fig7_newton(data, outdir)
    fig8_gummel(data, outdir)
    if verbose:
        print_summary(data)
        print(f"完成，图已输出到 {os.path.abspath(outdir)}")
    return outdir


if __name__ == "__main__":
    make_figures()
