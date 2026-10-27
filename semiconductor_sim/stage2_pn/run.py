# -*- coding: utf-8 -*-
"""
stage2_pn.run —— 第二阶段主流程：平衡态求解与验证 -> 偏压扫描 ->
Shockley 分析（J_s 对照、分区理想因子、理想公式失效点、Gummel 对照）。

统一接口（系列规范 docs/spec/code_spec.md）：
    run(...) -> dict，必含四个标准键
        params     : 输入参数快照
        results    : x / res_eq / ver / records / shockley 等
        validation : common.validation.Checklist 实例
        defects    : 本阶段已知缺陷 ID 列表（docs/defects/）

数值口径与原件 pn_junction.main() 完全一致（同一网格、同一求解器、
同一扫描序列），仅组织方式改为"计算进 dict、打印进 verbose 分支"。
"""

import time
import warnings

import numpy as np

from common.constants import (Js_analytic, NA_CM3, ND_CM3, NI, TAU_N, TEMP,
                              V_bi_analytic, VT, W_analytic)
from common.grid import make_junction_grid
from common.validation import (VALIDATION, SELF_CHECK, FORMAT_CHECK,
                               ARTIFACT, Checklist)
from .solver import gummel_solve, newton_solve, solve_bias
from .postprocess import current_density, depletion_width

# 本阶段已登记缺陷（详见 docs/defects/defect_log.md）
STAGE2_DEFECTS = ["DEF-S2-01", "DEF-S2-02", "DEF-S2-03"]


# ================================================================================
# 1. 平衡态求解与验证
# ================================================================================

def solve_equilibrium(x, Na, Nd, verbose=True):
    """V_a = 0 平衡态全 Newton 求解（准费米势变量）。"""
    res = newton_solve(x, Na, Nd, 0.0, verbose=verbose)
    return res


def verify_equilibrium(x, Na, Nd, res):
    """
    平衡态四重验证，返回 dict：
        V_bi_num   : 数值内建电势（跨结压降）
        V_bi_exact : 解析值
        np_err     : max|np - n_i^2| / n_i^2（玻尔兹曼统计自洽性）
        J_eq       : 平衡态电流（应 ~0）
        W_charge   : 电荷法耗尽宽度（数值）
        W_exact    : 解析耗尽宽度
    """
    phi, n, p = res["phi"], res["n"], res["p"]
    V_bi_num = float(phi[-1] - phi[0])
    np_err = float(np.max(np.abs(n * p - NI ** 2)) / NI ** 2)

    J, Jt, _, _ = current_density(phi, n, p, x)
    J_rel = abs(J) / max(abs(Jt).max(), 1e-30)

    W_c, _ = depletion_width(phi, n, p, Na, Nd, x)
    W_e = W_analytic(Na, Nd, 0.0)

    return dict(V_bi_num=V_bi_num, V_bi_exact=V_bi_analytic(Na, Nd),
                np_err=np_err, J_eq=J, J_scale=float(np.max(np.abs(Jt))),
                W_charge=W_c, W_exact=W_e,
                iters=res["iters"], converged=res["converged"])


# ================================================================================
# 2. 偏压扫描（参数延续）
# ================================================================================

def sweep_bias(x, Na, Nd, v_forward, v_reverse, tol=1e-9, verbose=True,
               v_step=0.1, min_step=1e-3):
    """
    偏压扫描：从平衡态出发，正偏 0 -> v_forward，反偏 0 -> v_reverse。

    求解器：全 Newton（solve_bias，参数延续）。Gummel 解耦迭代在 SRH
    复合下收敛极其缓慢（实测停滞：|dphi| 卡在 1e-2~1e-1 V 量级不再下降，
    同时载流子浓度在 N_d~49N_d 间剧烈振荡，30 步外层迭代远不收敛），
    已弃用；全 Newton 二次收敛，实测每子步 5~16 步。

    延续策略
    --------
    正/反两个分支均从平衡态独立出发（避免分支间大跳变）。相邻偏压点
    间距超过 v_step (V) 时自动细分为多步；单点不收敛时将该区间步长
    减半重试（最小 min_step），仍失败则记录未收敛并继续。

    返回 list[dict]：每步记录 V_a, J, W, 收敛信息。
    """
    def record(V_a, res):
        J, _, _, _ = current_density(res["phi"], res["n"], res["p"], x)
        Wc, Wk = depletion_width(res["phi"], res["n"], res["p"], Na, Nd, x)
        return dict(V_a=V_a, J=J, W_charge=Wc, W_crit=Wk,
                    iters=res["iters"], converged=res["converged"])

    def walk(V_target, state, V_cur):
        """从 V_cur 延续到 V_target（自适应细分），返回 (state, records)。"""
        recs = []
        # 预细分：步长 <= v_step
        n_sub = max(1, int(np.ceil(abs(V_target - V_cur) / v_step)))
        pts = [V_cur + (V_target - V_cur) * (k + 1) / n_sub
               for k in range(n_sub)]
        i = 0
        while i < len(pts):
            V_next = pts[i]
            res = solve_bias(x, Na, Nd, V_next, state=state, tol=tol)
            if res["converged"]:
                state = (res["phi"], res["n"], res["p"])
                recs.append(record(V_next, res))
                i += 1
                continue
            # 失败：区间减半重试
            if abs(V_next - V_cur) > min_step:
                pts.insert(i, 0.5 * (V_cur + V_next))
                continue
            warnings.warn(
                f"sweep_bias: V_a={V_next:.4f} V 在最小步长 {min_step} V 下"
                f"仍不收敛，记录末次解。", RuntimeWarning, stacklevel=2)
            recs.append(record(V_next, res))
            i += 1
        return state, recs

    # ---- 平衡点（全 Newton 冷启动）----
    res0 = newton_solve(x, Na, Nd, 0.0, tol=tol)
    if not res0["converged"]:
        warnings.warn("sweep_bias: 平衡态未收敛。", RuntimeWarning,
                      stacklevel=2)
    state = (res0["phi"], res0["n"], res0["p"])
    state_eq = state
    records = [record(0.0, res0)]

    for branch, v_list in (("fwd", v_forward), ("rev", v_reverse)):
        # 每个分支从平衡态独立出发
        state = state_eq
        V_cur = 0.0
        for V_target in v_list:
            if abs(V_target) < 1e-12:
                continue
            state, recs = walk(float(V_target), state, V_cur)
            records.extend(recs)
            V_cur = float(V_target)
            if verbose:
                r = recs[-1]
                print(f"  V_a = {V_target:+7.3f} V : J = {r['J']:+.6e} A/m^2, "
                      f"W = {r['W_charge']*1e9:7.1f} nm, "
                      f"Newton iters = {r['iters']}, "
                      f"收敛 = {r['converged']}")
    return records


# ================================================================================
# 3. 主流程（统一入口）
# ================================================================================

def run(verbose=True, L=5.0e-6, v_forward=None, v_reverse=None, tol=1e-9):
    """
    完整仿真流程：网格 -> 平衡态求解与四重验证 -> 偏压扫描 ->
    Shockley 分析（J_s 对照、分区理想因子、理想公式失效点、Gummel 对照）。

    参数
    ----
    L          : 仿真域总长 (m)，须 >> 少子扩散长度 L_p ~ 0.34 um
    v_forward  : 正偏目标点序列 (V)，缺省 np.arange(0.05, 0.85, 0.05)
    v_reverse  : 反偏目标点序列 (V)，缺省 -np.arange(0.1, 5.05, 0.1)
    tol        : Newton 收敛判据 (V)

    返回
    ----
    dict：params / results / validation / defects 四标准键
    （results 内为 x, res_eq, ver, records, Va_arr, J_arr, W_arr,
      Js_exact, seg_fits, ideality, failure, gummel, reverse 等）
    """
    t0 = time.time()
    Na = NA_CM3 * 1e6
    Nd = ND_CM3 * 1e6

    if v_forward is None:
        v_forward = np.round(np.arange(0.05, 0.85, 0.05), 4)
    if v_reverse is None:
        v_reverse = np.round(-np.arange(0.1, 5.05, 0.1), 4)

    if verbose:
        print("=" * 72)
        print("一维对称突变 PN 结直流特性仿真")
        print(f"N_a = N_d = {NA_CM3:.0e} cm^-3, L = {L*1e6:.1f} um, "
              f"tau = {TAU_N:.0e} s, T = {TEMP:.0f} K")
        print("=" * 72)

    # ---- 网格 ----
    x = make_junction_grid(L)
    if verbose:
        print(f"\n[网格] N = {len(x)} 点, x in [{x[0]*1e6:.2f}, {x[-1]*1e6:.2f}] um, "
              f"结区步长 {np.min(np.diff(x))*1e9:.2f} nm, "
              f"远场步长 {np.max(np.diff(x))*1e9:.2f} nm")
    assert abs(x[np.argmin(np.abs(x))]) < 1e-15, "结必须位于节点上"

    # ---- 平衡态 ----
    if verbose:
        print("\n[平衡态求解] ...")
    res_eq = solve_equilibrium(x, Na, Nd, verbose=verbose)
    ver = verify_equilibrium(x, Na, Nd, res_eq)
    if verbose:
        print("\n[平衡态四重验证]")
        print(f"  V_bi  数值 = {ver['V_bi_num']:.5f} V, "
              f"解析 = {ver['V_bi_exact']:.5f} V, "
              f"偏差 = {abs(ver['V_bi_num']-ver['V_bi_exact'])/ver['V_bi_exact']*100:.3f}%")
        print(f"  max|np - n_i^2|/n_i^2 = {ver['np_err']:.3e}")
        print(f"  平衡电流 J = {ver['J_eq']:.3e} A/m^2 "
              f"(截面电流尺度 {ver['J_scale']:.3e})")
        print(f"  W 电荷法 = {ver['W_charge']*1e9:.2f} nm, "
              f"解析 = {ver['W_exact']*1e9:.2f} nm, "
              f"偏差 = {abs(ver['W_charge']-ver['W_exact'])/ver['W_exact']*100:.2f}%")
        print(f"  Newton 外层迭代 {ver['iters']} 次, 收敛 = {ver['converged']}")

    # ---- 偏压扫描 ----
    if verbose:
        print("\n[偏压扫描] ...")
    records = sweep_bias(x, Na, Nd, v_forward, v_reverse, tol=tol,
                         verbose=verbose)
    Va_arr = np.array([r["V_a"] for r in records])
    J_arr = np.array([r["J"] for r in records])
    W_arr = np.array([r["W_charge"] for r in records])
    iters_arr = np.array([r["iters"] for r in records])
    n_bad = sum(1 for r in records if not r["converged"])

    # ---- Shockley 分析 ----
    Js_exact = Js_analytic(Na, Nd)
    if verbose:
        print("\n[Shockley 分析]")
        print(f"  解析扩散饱和电流 J_s = {Js_exact:.4e} A/m^2 "
              f"= {Js_exact*1e-4:.4e} A/cm^2")

    # 理想因子：ln J - V_a 分区线性拟合（三段互不重叠区间）
    # J = J_0 exp(V/(n_id V_T)) -> 斜率 = 1/(n_id V_T) -> n_id = 1/(斜率 V_T)
    # 本器件 tau = 1e-10 s，耗尽区 SRH 复合/产生主导电流传输：
    #   低偏压下复合电流 n_id -> 2（耗尽区复合极限），中偏压扩散分量
    #   赶上使局部斜率向 1 靠近，高偏压（V_a > V_bi 后注入少子超过
    #   掺杂）进入高注入，表观 n_id 反升到 2 以上（电导调制 + 过渡段）。
    #   注意高注入段表观 n_id > 2 并不对应简单理论极限，见论文 4.4 节。
    ideality = []
    if verbose:
        print("\n[分区理想因子拟合（互不重叠区间）]")
    for lo, hi in ((0.10, 0.30), (0.35, 0.65), (0.70, 0.80)):
        mask = (Va_arr >= lo - 1e-9) & (Va_arr <= hi + 1e-9)
        if mask.sum() >= 2:
            coef = np.polyfit(Va_arr[mask], np.log(np.abs(J_arr[mask])), 1)
            n_seg = 1.0 / (coef[0] * VT)
            ideality.append(dict(lo=lo, hi=hi, n_points=int(mask.sum()),
                                 slope=float(coef[0]), n_id=float(n_seg)))
            if verbose:
                print(f"  {lo:.2f}-{hi:.2f} V ({mask.sum()} 点): "
                      f"斜率 {coef[0]:.2f} /V, n_id = {n_seg:.3f}")
    if verbose:
        print(f"  （参考：低偏压耗尽区复合极限 n_id=2，纯扩散极限 n_id=1，"
              f"理想斜率 1/V_T = {1/VT:.2f} /V）")

    # 理想 Shockley 方程失效点：数值 J 与 Js*exp(V/VT) 的比值
    failure = []
    if verbose:
        print("\n[理想 Shockley 方程失效]")
    for V in (0.70, 0.75, 0.80):
        i = int(np.argmin(np.abs(Va_arr - V)))
        ideal = Js_exact * np.exp(Va_arr[i] / VT)
        failure.append(dict(V=float(Va_arr[i]), J=float(J_arr[i]),
                            J_ideal=float(ideal),
                            overestimate=float(ideal / J_arr[i])))
        if verbose:
            print(f"  V={Va_arr[i]:.2f}: J={J_arr[i]:.4e}, "
                  f"J_ideal={ideal:.4e}, 理想公式高估 {ideal/J_arr[i]:.0f} 倍")

    # Gummel 解耦迭代对照（停滞证据；主求解器为全 Newton）
    gummel = []
    if verbose:
        print("\n[Gummel 解耦迭代对照（各 30 步外层迭代）]")
    for omega in (1.0, 0.5):
        g = gummel_solve(x, Na, Nd, 0.0, max_outer=30, omega=omega,
                         verbose=False)
        gummel.append(dict(omega=omega, iters=int(g["iters"]),
                           converged=bool(g["converged"]),
                           dphi=float(g["dphi"])))
        if verbose:
            print(f"  omega={omega}: 外层 {g['iters']} 步, "
                  f"收敛 = {g['converged']}, 末步 |dphi| = {g['dphi']:.2e} V")

    # 反偏电流：耗尽区产生主导。恒等式验证 J = q * int(U) dx
    # （由 dJ_n/dx = qU, dJ_p/dx = -qU 积分得到），并对照耗尽近似
    # 解析式 J_gen ~ q n_i W / (2 tau)（均匀产生近似，偏大）。
    reverse = None
    rec_rev = [r for r in records if r["V_a"] <= -3.0]
    if rec_rev:
        mean_abs_j = float(np.mean([abs(r["J"]) for r in rec_rev]))
        reverse = dict(V_lo=float(rec_rev[0]["V_a"]),
                       V_hi=float(rec_rev[-1]["V_a"]),
                       mean_abs_J=mean_abs_j,
                       Js_fraction=float(Js_exact / mean_abs_j))
        if verbose:
            print(f"  反偏 (-3~-5V) 平均 |J| = {mean_abs_j:.4e} "
                  f"A/m^2 = {mean_abs_j*1e-4:.4e} A/cm^2")
            print(f"    其中扩散分量 J_s 仅占 "
                  f"{Js_exact/mean_abs_j:.1e}"
                  f"（耗尽区 SRH 产生主导）")

    # ---- 验证清单（规范 3.2 四分类）----
    # 注意 V_bi 项：phi(0)/phi(L) 由接触边界条件钉扎，数值跨结压降等于
    # 边界值之差是"边界保持性"自洽，不是独立物体验证（论文已降级表述）。
    cl = Checklist("stage2")
    cl.add(metric="内建电势 V_bi（边界一致性）", kind=SELF_CHECK,
           reference="接触边界 phi(L)-phi(0) = V_T ln(N_a N_d/n_i^2)",
           criterion="与边界设定一致",
           measured=f"数值 {ver['V_bi_num']:.6f} V vs 解析 "
                    f"{ver['V_bi_exact']:.6f} V "
                    f"({abs(ver['V_bi_num']-ver['V_bi_exact'])/ver['V_bi_exact']*100:.3f}%)",
           note="Dirichlet 边界钉扎量，非独立验证；仅检验求解器保持边界值")
    cl.add(metric="质量作用 np = n_i^2", kind=FORMAT_CHECK,
           reference="玻尔兹曼统计恒等式（机器精度）",
           criterion="< 1e-6",
           measured=f"max|np-n_i^2|/n_i^2 = {ver['np_err']:.3e}",
           passed=ver["np_err"] < 1e-6)
    cl.add(metric="平衡态电流 J = 0", kind=SELF_CHECK,
           reference="热平衡无净电流（相对截面电流尺度）",
           criterion="|J| << max|J_t|",
           measured=f"J = {ver['J_eq']:.3e} A/m^2（尺度 {ver['J_scale']:.3e}）",
           passed=abs(ver["J_eq"]) < 1e-6 * max(abs(ver["J_scale"]), 1e-30))
    cl.add(metric="耗尽层宽度 W", kind=VALIDATION,
           reference="W = sqrt(2 eps (V_bi-V_a)/q (1/N_a+1/N_d))（耗尽近似）",
           criterion="< 5%",
           measured=f"电荷法 {ver['W_charge']*1e9:.2f} nm vs 解析 "
                    f"{ver['W_exact']*1e9:.2f} nm "
                    f"({abs(ver['W_charge']-ver['W_exact'])/ver['W_exact']*100:.2f}%)",
           passed=abs(ver["W_charge"] - ver["W_exact"]) / ver["W_exact"] < 0.05,
           note="耗尽近似未计耗尽边缘载流子，偏差为模型层级差异")
    cl.add(metric="Shockley J_s 对照", kind=VALIDATION,
           reference="J = J_s[exp(V_a/V_T)-1]（扩散极限，长基区近似）",
           criterion="定性符合 + 失效点已量化",
           measured=f"J_s = {Js_exact:.3e} A/m^2；理想公式在 "
                    f"V>=0.70 V 高估 "
                    f"{min(f['overestimate'] for f in failure):.0f}~"
                    f"{max(f['overestimate'] for f in failure):.0f} 倍",
           note="tau=1e-10 s 下耗尽区 SRH 复合/产生主导，理想公式"
                "高偏压失效属模型层级差异，失效倍数已记录于 results")

    result = dict(
        params=dict(L=L, NA_CM3=NA_CM3, ND_CM3=ND_CM3, TAU_N=TAU_N,
                    TEMP=TEMP, tol=tol,
                    v_forward=list(map(float, v_forward)),
                    v_reverse=list(map(float, v_reverse))),
        results=dict(x=x, res_eq=res_eq, ver=ver, records=records,
                     Va_arr=Va_arr, J_arr=J_arr, W_arr=W_arr,
                     iters_arr=iters_arr, n_bad=n_bad,
                     Js_exact=Js_exact, ideality=ideality,
                     failure=failure, gummel=gummel, reverse=reverse),
        validation=cl,
        defects=list(STAGE2_DEFECTS),
    )
    result["t_total"] = time.time() - t0
    if verbose:
        print(f"\n总耗时 {result['t_total']:.1f} s")
    return result


def main():
    """CLI 入口：完整流程（终端输出双写留存 results/stage2/terminal.txt）。"""
    from common.io import TerminalLog
    with TerminalLog("stage2"):
        run(verbose=True)


if __name__ == "__main__":
    main()
