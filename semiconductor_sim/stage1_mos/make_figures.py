# -*- coding: utf-8 -*-
"""stage1_mos.make_figures —— 第一阶段全部插图（输出到 results/stage1/）。

图清单（论文编号 S1-F1 ~ S1-F9，文件名与论文引用一一对应）：
    fig0_structure.png            S1-F1  MOS 电容结构示意图
    fig1_band_diagram.png         S1-F2  能带图
    fig2_carrier_field.png        S1-F3  载流子/电场分布
    fig3_cv_curve.png             S1-F4  低频 C-V 特性
    fig4_qv_curve.png             S1-F5  Q-V 与表面势
    fig5_parametric_study.png     S1-F6  V_T 参数敏感性
    fig6_hf_vs_lf_cv.png          S1-F7  高频 vs 低频 C-V
    fig7_charge_decomposition.png S1-F8  电荷分解与表面势钉扎
    fig8_grid_convergence.png     S1-F9  网格收敛（由 run.analyze_grid_convergence 调出）
    fig9_grid_comparison.png            均匀 vs 非均匀（由 run.compare_grids 调出）

约定（系列图表规范 docs/spec/figure_spec.md）：所有函数接受 outdir 参数，
缺省写入 results/stage1/；文件名即上图清单，报告构建脚本按名引用。"""

import os

import numpy as np
import matplotlib.pyplot as plt

from common.constants import EC_MINUS_EI, EI_MINUS_EV, NI_CM3, VT
from common.io import ensure_results_dir


def _outdir(outdir):
    """outdir 缺省解析为 results/stage1/（系列统一结果目录）；显式传入时确保目录存在。"""
    if outdir is None:
        return ensure_results_dir("stage1")
    os.makedirs(outdir, exist_ok=True)
    return outdir
def _style():
    plt.rcParams.update({
        "figure.dpi": 130, "font.size": 11, "axes.grid": True, "grid.alpha": 0.3,
        "axes.linewidth": 1.0, "mathtext.fontset": "stix", "font.family": "STIXGeneral",
    })


def plot_band_diagram(res, outdir=None):
    """能带图：不同表面势下 Ec / Ei / Ev 相对费米能级的位置。"""
    outdir = _outdir(outdir)
    x = res["x"] * 1e6
    psi_F = res["psi_F"]
    cases = [(-0.30, "Accumulation"), (0.0, "Flat-band"), (0.35, "Depletion"),
             (res["phi_s_inv_actual"], f"Threshold (2psi_F={res['phi_s_inv']:.3f} V)"),
             (0.85, "Strong inversion")]
    colors = plt.cm.viridis(np.linspace(0.05, 0.85, len(cases)))

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    for (phi_s, label), col in zip(cases, colors):
        phi = res["phi_profiles"][round(float(phi_s), 3)]
        Ei = psi_F - phi
        Ec = Ei + EC_MINUS_EI
        Ev = Ei - EI_MINUS_EV
        ax.plot(x, Ec, color=col, lw=1.8, label=label)
        ax.plot(x, Ei, color=col, lw=1.2, ls="--", alpha=0.75)
        ax.plot(x, Ev, color=col, lw=1.8)

    ax.plot([], [], color="gray", lw=1.8, label="$E_C$")
    ax.plot([], [], color="gray", lw=1.2, ls="--", label="$E_i$")
    ax.plot([], [], color="gray", lw=1.8, label="$E_V$")
    ax.axhline(0, color="k", lw=1.0, ls=":", alpha=0.6)
    ax.text(x[-1] * 0.99, 0.02, "$E_F$", ha="right", va="bottom", fontsize=11)
    ax.set_xlabel("Depth $x$ ($\\mu$m)")
    ax.set_ylabel("Energy relative to $E_F$ (eV)")
    ax.set_title("MOS Capacitor Band Diagram (p-type Si, $N_A=10^{16}$ cm$^{-3}$)")
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(f"{outdir}/fig1_band_diagram.png", dpi=200)
    plt.close(fig)


def plot_carrier_field(res, outdir=None):
    """载流子浓度分布与电场分布。"""
    outdir = _outdir(outdir)
    x = res["x"] * 1e6
    Na_cm3 = res["Na_cm3"]
    cases = [(-0.30, "Accumulation"), (0.0, "Flat-band"),
             (res["phi_s_inv_actual"], "Threshold"), (0.85, "Strong inversion")]
    colors = plt.cm.plasma(np.linspace(0.1, 0.85, len(cases)))
    p0 = Na_cm3
    n0 = NI_CM3 ** 2 / p0

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8))
    ax = axes[0]
    for (phi_s, label), col in zip(cases, colors):
        phi = res["phi_profiles"][round(float(phi_s), 3)]
        p = p0 * np.exp(-phi / VT)
        n = n0 * np.exp(+phi / VT)
        ax.semilogy(x, np.maximum(p, 1e0), color=col, lw=1.8, label=f"p ({label})")
        ax.semilogy(x, np.maximum(n, 1e0), color=col, lw=1.8, ls="--", label=f"n ({label})")
    ax.axhline(Na_cm3, color="k", ls=":", lw=1.0, alpha=0.6)
    ax.text(x[-1] * 0.98, Na_cm3 * 1.7, "$N_A$", ha="right", fontsize=10)
    ax.set_xlabel("Depth $x$ ($\\mu$m)")
    ax.set_ylabel("Carrier concentration (cm$^{-3}$)")
    ax.set_title("Carrier Distribution")
    ax.set_ylim(1e0, max(1e20, Na_cm3 * 1e4))
    ax.legend(fontsize=7.5, ncol=2)

    ax = axes[1]
    for (phi_s, label), col in zip(cases, colors):
        phi = res["phi_profiles"][round(float(phi_s), 3)]
        E = -np.gradient(phi, res["x"]) * 1e-2
        ax.plot(x, E, color=col, lw=1.8, label=label)
    ax.axhline(0, color="k", lw=0.8, ls=":")
    ax.set_yscale("symlog", linthresh=1e4)
    ax.set_xlabel("Depth $x$ ($\\mu$m)")
    ax.set_ylabel("Electric field (V/cm)")
    ax.set_title("Electric Field Distribution (symlog)")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(f"{outdir}/fig2_carrier_field.png", dpi=200)
    plt.close(fig)


def plot_cv_curve(res, outdir=None):
    """低频 C-V 特性曲线。"""
    outdir = _outdir(outdir)
    Vg = res["Vg_arr"]
    C = res["C_lf_norm"]
    Vt = res["Vt_numeric"]

    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    ax.plot(Vg, C, color="#c0392b", lw=2.2, label="Numerical (quasi-static / low-freq)")
    ax.axvline(0.0, color="gray", ls="--", lw=1.0, alpha=0.7)
    ax.axvline(Vt, color="#2c3e50", ls="--", lw=1.2, alpha=0.9)
    ax.text(0.02, 0.30, "$V_{FB}$", transform=ax.transAxes, fontsize=11, color="gray")
    ax.text(0.80, 0.30, f"$V_T$ = {Vt:.3f} V", transform=ax.transAxes,
            fontsize=11, color="#2c3e50")
    # 平带电容标注（与表 2 第 3 行对应）
    i_fb = int(np.argmin(np.abs(Vg - 0.0)))
    ax.plot([Vg[i_fb]], [C[i_fb]], "o", ms=7, color="#c0392b", zorder=5)
    ax.annotate(f"$C_{{FB}}$ = {C[i_fb]:.4f} $C_{{ox}}$",
                xy=(Vg[i_fb], C[i_fb]), xytext=(Vg[i_fb] + 0.6, C[i_fb] + 0.28),
                fontsize=10, color="#c0392b",
                arrowprops=dict(arrowstyle="->", color="#c0392b", lw=1.2))
    ax.set_ylim(0.0, 1.15)
    ax.axvspan(Vg.min(), 0.0, ymax=0.88, alpha=0.10, color="#2980b9")
    ax.axvspan(0.0, Vt, ymax=0.88, alpha=0.10, color="#27ae60")
    ax.axvspan(Vt, Vg.max(), ymax=0.88, alpha=0.10, color="#8e44ad")
    ax.text(0.06, 0.92, "Accumulation", transform=ax.transAxes, fontsize=10, color="#2980b9")
    ax.text(0.47, 0.92, "Depletion", transform=ax.transAxes, fontsize=10, color="#27ae60")
    ax.text(0.82, 0.92, "Inversion", transform=ax.transAxes, fontsize=10, color="#8e44ad")
    ax.set_xlabel("Gate voltage $V_G$ (V)")
    ax.set_ylabel("Normalized capacitance $C/C_{ox}$")
    ax.set_title("MOS Capacitor C-V Characteristic ($N_A=10^{16}$ cm$^{-3}$, $t_{ox}$=20 nm)")
    ax.legend(loc="lower right", fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{outdir}/fig3_cv_curve.png", dpi=200)
    plt.close(fig)


def plot_qv_curve(res, outdir=None):
    """Q-V 特性与表面势关系。"""
    outdir = _outdir(outdir)
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8))
    ax = axes[0]
    ax.plot(res["Vg_arr"], res["Qs_arr"] * 1e2, color="#16a085", lw=2.0)
    ax.axhline(0, color="k", lw=0.8, ls=":")
    ax.axvline(res["Vt_numeric"], color="#2c3e50", ls="--", lw=1.2)
    ax.set_xlabel("Gate voltage $V_G$ (V)")
    ax.set_ylabel("Semiconductor charge $Q_s$ ($\\mu$C/cm$^2$)")
    ax.set_title("Q-V Characteristic")

    ax = axes[1]
    ax.plot(res["Vg_arr"], res["phi_s_arr"], color="#8e44ad", lw=2.0)
    ax.axhline(res["phi_s_inv"], color="#c0392b", ls="--", lw=1.2,
               label=f"$2\\psi_F$ = {res['phi_s_inv']:.3f} V")
    ax.set_xlabel("Gate voltage $V_G$ (V)")
    ax.set_ylabel("Surface potential $\\phi_s$ (V)")
    ax.set_title("Surface Potential vs Gate Voltage")
    ax.legend(fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{outdir}/fig4_qv_curve.png", dpi=200)
    plt.close(fig)


def plot_parametric_study(outdir=None):
    """V_T 随掺杂浓度与氧化层厚度的参数化研究。"""
    outdir = _outdir(outdir)
    # 延迟导入：run 模块在模块级导入本文件（绘图函数），此处反向引用
    # 必须推迟到调用期，避免循环导入（见 docs/spec/code_spec.md）
    from .run import run
    Na_list = np.array([1e15, 5e15, 1e16, 5e16, 1e17])
    tox_list = np.array([5.0, 10.0, 20.0, 40.0])
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8))

    ax = axes[0]
    colors = plt.cm.coolwarm(np.linspace(0.1, 0.9, len(tox_list)))
    for tox, col in zip(tox_list, colors):
        vts = [run(NA_CM3=Na, TOX_NM=tox, NGRID=601, verbose=False)["Vt_numeric"]
               for Na in Na_list]
        ax.semilogx(Na_list, vts, "o-", color=col, lw=1.8, ms=5,
                    label=f"$t_{{ox}}$ = {tox:.0f} nm")
    ax.set_xlabel("Substrate doping $N_A$ (cm$^{-3}$)")
    ax.set_ylabel("Threshold voltage $V_T$ (V)")
    ax.set_title("Threshold Voltage vs Substrate Doping")
    ax.legend(fontsize=9)

    ax = axes[1]
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(Na_list)))
    for Na, col in zip(Na_list, colors):
        vts = [run(NA_CM3=Na, TOX_NM=tox, NGRID=601, verbose=False)["Vt_numeric"]
               for tox in tox_list]
        ax.plot(tox_list, vts, "s-", color=col, lw=1.8, ms=5,
                label=f"$N_A$ = {Na:.0e} cm$^{{-3}}$")
    ax.set_xlabel("Oxide thickness $t_{ox}$ (nm)")
    ax.set_ylabel("Threshold voltage $V_T$ (V)")
    ax.set_title("Threshold Voltage vs Oxide Thickness")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(f"{outdir}/fig5_parametric_study.png", dpi=200)
    plt.close(fig)


def plot_hf_vs_lf_cv(res, outdir=None):
    """低频 vs 高频 C-V 对比（新增）。"""
    outdir = _outdir(outdir)
    Vg = res["Vg_arr"]
    Vt = res["Vt_numeric"]
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    ax.plot(Vg, res["C_lf_norm"], color="#c0392b", lw=2.2,
            label="Low frequency (quasi-static)")
    ax.plot(Vg, res["C_hf_norm"], color="#2980b9", lw=2.2, ls="--",
            label="High frequency (inversion charge frozen)")
    ax.axvline(Vt, color="#2c3e50", ls=":", lw=1.2, alpha=0.9)
    ax.text(Vt + 0.03, 0.60, f"$V_T$ = {Vt:.3f} V", fontsize=10, color="#2c3e50")
    ax.annotate("HF curve saturates:\ninversion charge cannot\nrespond",
                xy=(Vg.max() * 0.75, res["C_hf_norm"][-1]),
                xytext=(Vg.min() * 0.42, 0.78), fontsize=9, color="#2980b9",
                arrowprops=dict(arrowstyle="->", color="#2980b9", lw=1.2))
    ax.set_ylim(0.0, 1.15)
    ax.set_xlabel("Gate voltage $V_G$ (V)")
    ax.set_ylabel("Normalized capacitance $C/C_{ox}$")
    ax.set_title("Low-Frequency vs High-Frequency C-V Characteristics")
    ax.legend(fontsize=9.5, loc="lower right")
    fig.tight_layout()
    fig.savefig(f"{outdir}/fig6_hf_vs_lf_cv.png", dpi=200)
    plt.close(fig)


def plot_charge_decomposition(res, outdir=None):
    """反型电荷 / 耗尽电荷分解与表面势钉扎（新增）。"""
    outdir = _outdir(outdir)
    Vg = res["Vg_arr"]
    phi_s = res["phi_s_arr"]
    Vt = res["Vt_numeric"]

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8))
    ax = axes[0]
    m = phi_s > 0
    ax.plot(Vg[m], res["Qd_arr"][m] * 1e2, color="#2980b9", lw=2.2, label="$Q_d$ (depletion)")
    ax.plot(Vg[m], -res["Qn_arr"][m] * 1e2, color="#c0392b", lw=2.2, label="$-Q_n$ (inversion)")
    ax.plot(Vg[m], res["Qs_arr"][m] * 1e2, "k--", lw=1.4, label="$Q_s$ (total)")
    ax.axvline(Vt, color="#2c3e50", ls=":", lw=1.2)
    ax.text(Vt + 0.02, ax.get_ylim()[1] * 0.80, "$V_T$", fontsize=10, color="#2c3e50")
    ax.set_xlabel("Gate voltage $V_G$ (V)")
    ax.set_ylabel("Charge ($\\mu$C/cm$^2$)")
    ax.set_title("Charge Decomposition beyond Threshold")
    ax.legend(fontsize=9)

    ax = axes[1]
    ax.plot(Vg, phi_s, color="#8e44ad", lw=2.2)
    ax.axhline(res["phi_s_inv"], color="#c0392b", ls="--", lw=1.4,
               label=f"$2\\psi_F$ = {res['phi_s_inv']:.3f} V (pinning level)")
    ax.fill_between(Vg, res["phi_s_inv"], phi_s,
                    where=phi_s > res["phi_s_inv"], alpha=0.12, color="#c0392b")
    ax.annotate("pinned: extra $\\Delta V_G$\nscreened by $\\Delta Q_n$",
                xy=(Vg.max() * 0.88, res["phi_s_inv"] + 0.02),
                xytext=(Vg.min() * 0.45, 0.80), fontsize=9,
                arrowprops=dict(arrowstyle="->", color="#7f8c8d", lw=1.2))
    ax.set_xlabel("Gate voltage $V_G$ (V)")
    ax.set_ylabel("Surface potential $\\phi_s$ (V)")
    ax.set_title("Surface-Potential Pinning in Strong Inversion")
    ax.legend(fontsize=9, loc="lower right")
    fig.tight_layout()
    fig.savefig(f"{outdir}/fig7_charge_decomposition.png", dpi=200)
    plt.close(fig)


def plot_structure(res, outdir=None):
    """图 1：MOS 电容结构示意图（教材风格截面图）。"""
    outdir = _outdir(outdir)
    _style()
    fig, ax = plt.subplots(figsize=(9.2, 4.4))
    tox_nm = res["tox_nm"]
    W_nm = res["W_max_analytic"] * 1e9
    L_nm = 400.0

    # 各层（示意厚度，不按比例）
    y_gate = (1.00, 1.22)
    y_ox = (0.72, 1.00)
    y_si = (-0.20, 0.72)
    ax.add_patch(plt.Rectangle((0, y_gate[0]), 10, y_gate[1] - y_gate[0],
                               fc="#95a5a6", ec="k", lw=1.2))
    ax.add_patch(plt.Rectangle((0, y_ox[0]), 10, y_ox[1] - y_ox[0],
                               fc="#d6eaf8", ec="k", lw=1.2, hatch="///"))
    ax.add_patch(plt.Rectangle((0, y_si[0]), 10, y_si[1] - y_si[0],
                               fc="#fdebd0", ec="k", lw=1.2))

    # 耗尽区（示意）
    ax.add_patch(plt.Rectangle((0, y_si[1] - 0.34), 10, 0.34,
                               fc="#fad7a0", ec="#b9770e", lw=1.0, ls="--"))
    # 反型层（示意）
    ax.add_patch(plt.Rectangle((0, y_si[1] - 0.05), 10, 0.05,
                               fc="#5dade2", ec="#1a5276", lw=0.8))

    # 标注
    ax.annotate("Metal gate (G)", xy=(10.15, 1.11), fontsize=11, va="center")
    ax.annotate(f"SiO$_2$  $t_{{ox}}$ = {tox_nm:.0f} nm", xy=(10.15, 0.86),
                fontsize=11, va="center")
    ax.annotate("Inversion layer $Q_n$ (electrons)", xy=(10.15, 0.69),
                fontsize=10, va="center", color="#1a5276")
    ax.annotate("Depletion region  $W$", xy=(10.15, 0.50), fontsize=11,
                va="center", color="#873600")
    ax.annotate(f"p-type Si  $N_A$ = {res['Na_cm3']:.0e} cm$^{{-3}}$",
                xy=(10.15, 0.18), fontsize=11, va="center")

    # 尺寸标注线
    ax.annotate("", xy=(5.0, 1.00), xytext=(5.0, 0.72),
                arrowprops=dict(arrowstyle="<->", color="#333", lw=1.1))
    ax.text(5.25, 0.86, "$t_{ox}$", fontsize=11)
    ax.annotate("", xy=(5.0, 0.72), xytext=(5.0, 0.72 - 0.34),
                arrowprops=dict(arrowstyle="<->", color="#333", lw=1.1))
    ax.text(5.25, 0.55, "$W$", fontsize=11)
    ax.annotate("", xy=(0.0, -0.20), xytext=(10.0, -0.20),
                arrowprops=dict(arrowstyle="<->", color="#333", lw=1.0))
    ax.text(4.4, -0.32, "1-$\\mu$m simulation domain $L$", fontsize=9,
            ha="center", color="#333")

    # 偏压与电荷方向
    ax.annotate("$V_G$", xy=(5, 1.38), fontsize=13, ha="center", color="#c0392b")
    ax.plot([5, 5], [1.30, 1.22], color="#c0392b", lw=1.6)
    ax.plot([4.6, 5.4], [1.30, 1.30], color="#c0392b", lw=1.6)
    ax.text(2.0, 1.45, f"$\\phi_s$ = 0.85 V (example)", fontsize=10,
            color="#8e44ad")

    ax.text(2.0, -0.13, "field $\\mathcal{E}(x)$", fontsize=9.5,
            color="#666", style="italic")
    ax.annotate("", xy=(3.2, -0.05), xytext=(2.0, -0.13),
                arrowprops=dict(arrowstyle="->", color="#666", lw=1.2))

    ax.set_xlim(-0.4, 17.5)
    ax.set_ylim(-0.55, 1.62)
    ax.set_aspect("auto")
    ax.axis("off")
    ax.set_title("MOS Capacitor Structure (schematic, not to scale)",
                 fontsize=12, pad=10)
    fig.tight_layout()
    fig.savefig(f"{outdir}/fig0_structure.png", dpi=200)
    plt.close(fig)


def plot_grid_convergence(results, p_VT, outdir=None):
    """图 8：网格收敛性（log-log 离散误差 vs 步长 h，二阶参考线）。

    results / p_VT 来自 run.analyze_grid_convergence 的计算结果
    （本函数只绘图，不重算）。
    """
    outdir = _outdir(outdir)
    _style()
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    hs = np.array([r["h"] for r in results]) * 1e9
    ev = np.array([max(r["err_v"], 1e-12) for r in results])
    ew = np.array([max(r["err_w"], 1e-12) for r in results])
    ax.loglog(hs, ev, "o-", color="#c0392b", lw=1.9, ms=6, label=r"$|V_T - V_T^{ref}|$")
    ax.loglog(hs, ew, "s-", color="#2980b9", lw=1.9, ms=6, label=r"$|W_{max} - W_{max}^{ref}|$")
    href = np.array([hs[1], hs[1] / 2])
    scale = ev[1] if ev[1] > 0 else 1e-6
    ax.loglog(href, scale * (href / href[0]) ** 2, "--", color="#555", lw=1.4,
              label=r"$\propto h^2$ reference")
    ax.set_xlabel("Mesh spacing $h$ (nm)")
    ax.set_ylabel("Discretization error")
    ax.set_title(f"Grid Convergence Study (observed order $p$ ≈ {p_VT:.2f})")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(f"{outdir}/fig8_grid_convergence.png", dpi=200)
    plt.close(fig)

def plot_grid_comparison(ru, rn, outdir=None):
    """图 9：均匀 801 点 vs 非均匀网格（表面势剖面与间距分布）。

    ru / rn 为 run.compare_grids 返回的两个 run() 结果 dict。
    """
    outdir = _outdir(outdir)
    _style()
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.6))
    phi_s = 0.85
    phi_u = ru["phi_profiles"][round(phi_s, 3)]
    phi_n = rn["phi_profiles"][round(phi_s, 3)]
    ax = axes[0]
    m = ru["x"] < 0.25e-6
    ax.plot(ru["x"][m] * 1e9, phi_u[m], "o-", ms=2.5, lw=1.0,
            color="#c0392b", label="Uniform, $N$=801 ($h$=1.25 nm)")
    mn = rn["x"] < 0.25e-6
    ax.plot(rn["x"][mn] * 1e9, phi_n[mn], "-", lw=1.6,
            color="#2980b9",
            label=f"Non-uniform, $N$={len(rn['x'])} (surface $h_0$=0.31 nm)")
    ax.set_xlabel("Depth $x$ (nm)")
    ax.set_ylabel("Surface potential $\\phi$ (V)")
    ax.set_title("Potential near the Surface ($\\phi_s$ = 0.85 V)")
    ax.legend(fontsize=9)

    ax = axes[1]
    ax.semilogy(ru["x"][:-1] * 1e9, np.diff(ru["x"]) * 1e9, "o-", ms=2,
                color="#c0392b", label="Uniform 801 pts")
    ax.semilogy(rn["x"][:-1] * 1e9, np.diff(rn["x"]) * 1e9, "-",
                color="#2980b9", label=f"Non-uniform {len(rn['x'])} pts")
    ax.set_xlabel("Position $x$ (nm)")
    ax.set_ylabel("Local spacing $\\Delta x$ (nm)")
    ax.set_title("Mesh Spacing Distribution")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(f"{outdir}/fig9_grid_comparison.png", dpi=200)
    plt.close(fig)


# ================================================================================
# 编排入口
# ================================================================================

def make_figures(res=None, outdir=None):
    """生成第一阶段 res 相关插图 fig0-fig7 到 results/stage1/。

    fig8 / fig9（网格收敛、网格对比）由 run.analyze_grid_convergence /
    run.compare_grids 在计算分析时调出绘图函数生成。
    res 缺省时现场调用 run.run()。
    """
    outdir = _outdir(outdir)
    if res is None:
        from .run import run
        res = run()
    _style()
    plot_structure(res, outdir)
    plot_band_diagram(res, outdir)
    plot_carrier_field(res, outdir)
    plot_cv_curve(res, outdir)
    plot_qv_curve(res, outdir)
    plot_hf_vs_lf_cv(res, outdir)
    plot_charge_decomposition(res, outdir)
    plot_parametric_study(outdir)
    print(f"第一阶段插图已输出到 {outdir}/（fig0-fig7；fig8/fig9 见网格分析）")
    return outdir


if __name__ == "__main__":
    make_figures()
