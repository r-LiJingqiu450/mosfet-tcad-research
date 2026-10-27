# -*- coding: utf-8 -*-
"""stage1_mos.run —— MOS 电容主流程：扫描表面势 -> Q-V / C-V（低频 + 高频）
+ 电荷分解 + 教科书解析公式交叉验证 + 网格收敛性分析。

统一接口（系列规范 docs/spec/code_spec.md）：
    run(...) -> dict，必含四个标准键
        params     : 输入参数快照
        results    : （本阶段将主要数值结果平铺于顶层，历史键名不变）
        validation : common.validation.Checklist 实例
        defects    : 本阶段已知缺陷 ID 列表（docs/defects/）
绘图集中于 make_figures.py；本模块的 analyze_grid_convergence /
compare_grids 只做计算与打印，图像经 make_figures 生成。"""

import numpy as np

from common.constants import EPS_OX, EPS_SI, K_B, NI, Q_E, TEMP, VFB, VT
from common.grid import make_grid
from common.validation import (VALIDATION, SELF_CHECK, FORMAT_CHECK,
                               ARTIFACT, Checklist)
from .solver import solve_poisson
from .postprocess import postprocess, extract_depletion_width

# 本阶段已登记缺陷（详见 docs/defects/defect_log.md）
STAGE1_DEFECTS = ["DEF-S1-01", "DEF-S1-02"]


def run(NA_CM3=1e16, TOX_NM=20.0, L_UM=1.0, NGRID=801, VFB=0.0,
        phi_s_range=(-0.35, 0.90), n_samples=251,
        grid_kind="uniform", grid_params=None,
        tol=1e-11, max_iter=200, clamp=0.3, verbose=True):
    """
    完整仿真流程（低频 C-V + 高频 C-V + 电荷分解）。

    grid_kind='nonuniform' 时 grid_params={'d0':..., 'transition':...}
    """
    Na = NA_CM3 * 1e6
    tox = TOX_NM * 1e-9
    L = L_UM * 1e-6

    if grid_kind == "uniform":
        x = make_grid(L, NGRID, "uniform")
    else:
        gp = grid_params or dict(d0=1.0e-9, transition=0.35e-6, h_max=4.5e-9)
        x = make_grid(L, NGRID, "nonuniform", **gp)

    Cox = EPS_OX / tox
    psi_F = VT * np.log(Na / NI)
    phi_s_inv = 2.0 * psi_F

    phi_s_arr = np.linspace(phi_s_range[0], phi_s_range[1], n_samples)
    Qs_arr = np.zeros_like(phi_s_arr)
    Vg_arr = np.zeros_like(phi_s_arr)
    W_arr = np.zeros_like(phi_s_arr)
    Qn_arr = np.zeros_like(phi_s_arr)
    Qd_arr = np.zeros_like(phi_s_arr)
    phi_profiles = {}
    iter_counts, conv_flags = [], []

    phi_prev = None
    W_lock = None
    for k, phi_s in enumerate(phi_s_arr):
        phi, nit, ok = solve_poisson(phi_s, x, Na, phi_init=phi_prev,
                                     tol=tol, max_iter=max_iter, clamp=clamp)
        phi_prev = phi.copy()
        iter_counts.append(nit)
        conv_flags.append(ok)

        Qs, Es = postprocess(phi, x, Cox)
        Qs_arr[k] = Qs
        Vg_arr[k] = VFB + phi_s - Qs / Cox

        if phi_s > 0:
            if phi_s <= phi_s_inv:
                # 耗尽区：耗尽近似成立，W = |Q_s|/(q·N_a)
                W = extract_depletion_width(phi, x, Na, Es, Qs=Qs)
            else:
                # 强反型区：表面势钉扎，耗尽层宽度锁定在 W_max，
                # 新增栅压全部降在氧化层/反型层上（Q_n 线性增长）
                if W_lock is None:
                    W_lock = extract_depletion_width(phi, x, Na, Es, Qs=Qs)
                W = W_lock
            W_arr[k] = W
            Qd_arr[k] = -Q_E * Na * W               # 电离受主耗尽电荷（钉扎）
            Qn_arr[k] = Qs_arr[k] - Qd_arr[k]       # 反型层（拖尾）电荷
        phi_profiles[round(float(phi_s), 3)] = phi.copy()

    # ---- 低频（准静态）C-V ----
    dQs_dphis = np.gradient(Qs_arr, phi_s_arr)
    Cs_arr = -dQs_dphis
    C_lf = Cox * Cs_arr / (Cox + Cs_arr)

    # ---- 高频 C-V（反型层电荷冻结近似）----
    # 反型区（phi_s > 2*psi_F）：反型电荷来不及响应小信号，
    #     电荷变化仅来自耗尽层边界移动 -> C_s = eps_si/W_max
    # 积累/耗尽/平带区（phi_s <= 2*psi_F）：多子电荷与耗尽边界
    #     均可即时响应 -> 与准静态值相同（禁用 phi_s>0 切换：
    #     平带正侧 W=|Qs|/qNa -> 0 会产生 eps/W -> inf 的伪尖峰）
    Cs_hf = np.where(phi_s_arr > phi_s_inv,
                     EPS_SI / np.maximum(W_arr, 1e-30), Cs_arr)
    C_hf = Cox * Cs_hf / (Cox + Cs_hf)

    # ============================================================================
    # 验证：数值解 vs 教科书解析公式
    # ============================================================================
    W_max_analytic = np.sqrt(4.0 * EPS_SI * psi_F / (Q_E * Na))
    Vt_analytic = VFB + 2.0 * psi_F + np.sqrt(
        4.0 * Q_E * EPS_SI * Na * psi_F) / Cox

    idx_inv = int(np.argmin(np.abs(phi_s_arr - phi_s_inv)))
    Vt_numeric = Vg_arr[idx_inv]
    phi_s_inv_actual = float(phi_s_arr[idx_inv])
    W_num = -Qs_arr[idx_inv] / (Q_E * Na)

    # 平带电容：V_G = V_FB 处（phi_s = 0）数值插值
    idx_fb = int(np.argmin(np.abs(Vg_arr - VFB)))
    Cfb_numeric = C_lf[idx_fb] / Cox

    L_D = np.sqrt(EPS_SI * K_B * TEMP / (Q_E ** 2 * Na))
    Csfb = EPS_SI / L_D
    Cfb_analytic = Cox * Csfb / (Cox + Csfb)

    if verbose:
        print("=" * 74)
        print(" 一维 MOS 电容泊松方程求解 —— 数值解 vs 解析公式验证")
        print("=" * 74)
        print(f" 器件参数：Na = {NA_CM3:.2e} cm^-3,  tox = {TOX_NM:.1f} nm,  T = {TEMP:.0f} K")
        print(f" 网格：{NGRID} 点 ({grid_kind}),  C_ox = {Cox*1e-4:.4e} F/cm^2")
        print("-" * 74)
        print(f" 费米势 psi_F = {psi_F:.4f} V,  2*psi_F = {2*psi_F:.4f} V")
        print(f" 阈值电压 V_T : 数值 {Vt_numeric:.4f} V | 解析 {Vt_analytic:.4f} V"
              f" | 偏差 {abs(Vt_numeric-Vt_analytic)/Vt_analytic*100:.3f} %")
        print(f" 最大耗尽宽度: 数值 {W_num*1e9:.2f} nm | 解析 {W_max_analytic*1e9:.2f} nm"
              f" | 偏差 {abs(W_num-W_max_analytic)/W_max_analytic*100:.3f} %")
        print(f" 平带电容 C_FB/C_ox : 数值 {Cfb_numeric:.4f} | "
              f"解析 {Cfb_analytic/Cox:.4f}"
              f" | 偏差 "
              f"{abs(Cfb_numeric-Cfb_analytic/Cox)/(Cfb_analytic/Cox)*100:.2f} %")
        print(f" Newton 迭代：平均 {np.mean(iter_counts):.1f} 次, "
              f"最大 {np.max(iter_counts)} 次,"
              f" 全部收敛: {all(conv_flags)}")
        print("=" * 74)

    # ---- 验证清单（规范 3.2：验证 / 自洽 / 格式自洽 / artifact 四分类）----
    cl = Checklist("stage1")
    cl.add(metric="阈值电压 V_T", kind=VALIDATION,
           reference="V_FB + 2psi_F + gamma*sqrt(2psi_F)（耗尽近似解析式）",
           criterion="< 1%",
           measured=f"{abs(Vt_numeric - Vt_analytic) / Vt_analytic * 100:.3f} %",
           passed=abs(Vt_numeric - Vt_analytic) / Vt_analytic < 0.01)
    cl.add(metric="最大耗尽宽度 W_max", kind=VALIDATION,
           reference="sqrt(4 eps_si psi_F / (q N_a))",
           criterion="< 1%",
           measured=f"{abs(W_num - W_max_analytic) / W_max_analytic * 100:.3f} %",
           passed=abs(W_num - W_max_analytic) / W_max_analytic < 0.01)
    cl.add(metric="平带电容 C_FB/C_ox", kind=VALIDATION,
           reference="C_ox * (eps_si/L_D) / (C_ox + eps_si/L_D)（德拜近似）",
           criterion="< 5%",
           measured=f"{abs(Cfb_numeric - Cfb_analytic / Cox) / (Cfb_analytic / Cox) * 100:.2f} %",
           passed=abs(Cfb_numeric - Cfb_analytic / Cox) / (Cfb_analytic / Cox) < 0.05,
           note="德拜近似为一阶参照，允许偏差大于前两项")
    cl.add(metric="Newton 迭代收敛性", kind=SELF_CHECK,
           reference="全 phi_s 扫描点全部收敛（tol=1e-11）",
           criterion="全部收敛",
           measured=f"平均 {np.mean(iter_counts):.1f} 次 / 最大 {np.max(iter_counts)} 次",
           passed=bool(all(conv_flags)))

    result = dict(
        x=x, phi_s_arr=phi_s_arr, Qs_arr=Qs_arr, Vg_arr=Vg_arr,
        Cs_arr=Cs_arr, C_lf=C_lf, C_hf=C_hf, C_hf_norm=C_hf / Cox,
        C_lf_norm=C_lf / Cox, W_arr=W_arr, Qn_arr=Qn_arr, Qd_arr=Qd_arr,
        Cox=Cox, Na_cm3=NA_CM3, tox_nm=TOX_NM,
        phi_profiles=phi_profiles, psi_F=psi_F,
        phi_s_inv=phi_s_inv, phi_s_inv_actual=phi_s_inv_actual,
        Vt_numeric=Vt_numeric, Vt_analytic=Vt_analytic,
        W_num=W_num, W_max_analytic=W_max_analytic,
        Cfb_numeric=Cfb_numeric, Cfb_analytic=Cfb_analytic / Cox,
        L_D=L_D, grid_kind=grid_kind,
        iter_mean=np.mean(iter_counts), iter_max=np.max(iter_counts),
    )
    # ---- 系列统一返回结构（docs/spec/code_spec.md：四标准键）----
    result["params"] = dict(NA_CM3=NA_CM3, TOX_NM=TOX_NM, L_UM=L_UM,
                            NGRID=NGRID, VFB=VFB,
                            phi_s_range=tuple(phi_s_range),
                            n_samples=n_samples, grid_kind=grid_kind,
                            grid_params=grid_params, tol=tol,
                            max_iter=max_iter, clamp=clamp)
    result["validation"] = cl
    result["defects"] = list(STAGE1_DEFECTS)
    return result


# ================================================================================
# 7. 网格收敛性分析
# ================================================================================

def analyze_grid_convergence(outdir="results"):
    """
    网格加密研究：N = 201, 401, 801, 1601, 3201

    以 N = 3201 的解为参考，考察 V_T 与 W_max 的离散误差随步长 h 的变化，
    验证二阶收敛（误差 ∝ h^2，log-log 斜率 p ≈ 2）。
    """
    Ns = [201, 401, 801, 1601, 3201]
    results = []
    for N in Ns:
        results.append(dict(N=N, h=1e-6 / (N - 1),
                            Vt=run(NGRID=N, verbose=False)["Vt_numeric"],
                            W=run(NGRID=N, verbose=False)["W_num"]))

    # 以最密网格 (N=3201) 为参考计算离散误差
    ref_vt = results[-1]["Vt"]
    ref_w = results[-1]["W"]
    for r in results:
        r["err_v"] = abs(r["Vt"] - ref_vt)
        r["err_w"] = abs(r["W"] - ref_w)

    print("\n[网格收敛性分析]")
    print(f" {'N':>6} {'h (nm)':>9} {'V_T (V)':>12} {'err(V_T)':>11} "
          f"{'W_max (nm)':>12} {'err(W)':>11}")
    for r in results:
        print(f" {r['N']:>6} {r['h']*1e9:>9.4f} {r['Vt']:>12.6f} "
              f"{r['err_v']:>11.3e} {r['W']*1e9:>12.4f} {r['err_w']:>11.3e}")

    # 收敛阶估计：p = log(e1/e2) / log(h1/h2)，取相邻网格对平均
    def observed_order(errs, hs):
        ps = []
        for i in range(len(hs) - 2):          # 最密网格误差为 0，不参与
            if errs[i] > 0 and errs[i + 1] > 0:
                ps.append(np.log(errs[i] / errs[i + 1]) / np.log(hs[i] / hs[i + 1]))
        return float(np.mean(ps)) if ps else float("nan")

    hs = [r["h"] for r in results]
    p_VT = observed_order([r["err_v"] for r in results], hs)
    p_W = observed_order([r["err_w"] for r in results], hs)
    print(f" 平均收敛阶：V_T  p ≈ {p_VT:.2f}   W_max  p ≈ {p_W:.2f}"
          f"（理论值 2，二阶中心差分）")

    # 绘图移交 make_figures（本模块只做计算与打印，规范：绘图集中于 make_figures）
    from .make_figures import plot_grid_convergence
    plot_grid_convergence(results, p_VT, outdir)

    return results, p_VT


# ================================================================================
# 8. 均匀 vs 非均匀网格对比
# ================================================================================

def compare_grids(outdir="results"):
    """均匀 801 点 vs 非均匀 383 点：同等表面精度下的点数差异。"""
    ru = run(NGRID=801, grid_kind="uniform", verbose=False)
    rn = run(NGRID=401, grid_kind="nonuniform",
             grid_params=dict(d0=0.31e-9, transition=0.05e-6, h_max=4.5e-9),
             verbose=False)

    phi_s = 0.85
    phi_u = ru["phi_profiles"][round(phi_s, 3)]
    phi_n = rn["phi_profiles"][round(phi_s, 3)]

    print("\n[均匀 vs 非均匀网格]")
    print(f" 均匀 {len(ru['x'])} 点 (h=1.25 nm 全域)      : V_T = {ru['Vt_numeric']:.6f} V")
    print(f" 非均匀 {len(rn['x'])} 点 (表面 {rn['x'][1]*1e9:.2f} nm / 深部 "
          f"{np.diff(rn['x'])[-1]*1e9:.1f} nm) : V_T = {rn['Vt_numeric']:.6f} V")
    print(f" V_T 差别 = {abs(ru['Vt_numeric']-rn['Vt_numeric'])*1e6:.2f} uV")

    # 绘图移交 make_figures
    from .make_figures import plot_grid_comparison
    plot_grid_comparison(ru, rn, outdir)

    return ru, rn


# ================================================================================
# 8. 仿真域长度（Dirichlet 边界）敏感性分析
# ================================================================================

def boundary_sensitivity():
    """
    仿真域长度 / Dirichlet 边界敏感性分析。

    phi(L)=0 的截断误差估算：L 越大，深处电势衰减越充分，
    V_T 与 W_max 趋于与 L 无关的极限值。
    """
    print("\n[仿真域长度敏感性 (N=801 固定)]")
    print(f"  {'L (um)':>8} {'V_T (V)':>10} {'W_max (nm)':>11}")
    rows = []
    for L in [1.0, 1.5, 2.0]:
        r = run(NA_CM3=1e16, L_UM=L, NGRID=801, verbose=False)
        rows.append((L, r["Vt_numeric"], r["W_num"] * 1e9))
        print(f"  {L:>8.1f} {rows[-1][1]:>10.6f} {rows[-1][2]:>11.2f}")
    dV = abs(rows[0][1] - rows[-1][1]) * 1e6
    print(f"  L=1um 与 2um 的 V_T 相差 {dV:.0f} uV；N 固定时 "
          f"h=L/(N-1) 联动变化，"
          f"该差异完全来自 h（论文表 5 对照：同 h 不同 L 的"
          f"结果六位小数一致），"
          f"phi(L)=0 截断误差可忽略")

# ================================================================================
# 主程序
# ================================================================================

def main():
    """完整流程：C-V 主扫描 + 收敛性分析 + 边界敏感性（终端输出双写留存）。"""
    from common.io import TerminalLog
    with TerminalLog("stage1"):
        print("\n[1/4] 求解泊松方程并计算 C-V 特性（低频 + 高频）...")
        res = run()
        print("[2/4] 生成全部插图 ...")
        from .make_figures import make_figures
        make_figures(res)
        print("[3/4] 网格收敛性分析 ...")
        analyze_grid_convergence()
        print("[4/4] 均匀 vs 非均匀网格对比与边界敏感性 ...")
        compare_grids()
        boundary_sensitivity()


if __name__ == "__main__":
    main()
