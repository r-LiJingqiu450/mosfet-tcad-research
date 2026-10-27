# -*- coding: utf-8 -*-
"""
tests/test_stage_closure.py —— 阶段闭环测试（validation_spec.md 第 3 节）。

    1. S2 平衡退化：V_a = 0 时数值解给出 V_bi、np = n_i^2、J = 0；
    2. S3 V_D = 0 退化 S1：表格 phi_n = 0 列与第一阶段 MOS C-V 求解器
       在同一网格上的 Q_s 一致（逐点退化到机器精度）；生产口径
       （fig_cv_loop：240 点表 vs 120 点扫描）全程最大等效栅压差
       0.550 mV @ V_G=1.40 V，判据 < 1 mV。
"""

import numpy as np
import pytest

from semiconductor_sim.common import constants as C
from semiconductor_sim.common.grid import make_surface_grid
from semiconductor_sim.stage1_mos.solver import solve_poisson
from semiconductor_sim.stage1_mos.postprocess import postprocess as s1_postprocess
from semiconductor_sim.stage2_pn.run import run as s2_run
from semiconductor_sim.stage3_mosfet.solver import solve_slice, surface_charge
from semiconductor_sim.stage3_mosfet.table import SliceTable


# ============================================================================
# 1. S2 平衡退化
# ============================================================================

@pytest.fixture(scope="module")
def s2_result():
    """小域 + 两个偏压点的 S2 完整流程（够快，验证条目齐全）。"""
    return s2_run(verbose=False, L=2.0e-6, v_forward=[0.5], v_reverse=[-1.0])


def test_s2_equilibrium_degeneration(s2_result):
    """V_a=0 -> V_bi、np=n_i^2、J=0（四标准键与验证条目齐全）。"""
    assert set(["params", "results", "validation", "defects"]) <= set(
        s2_result)
    cl = s2_result["validation"]
    assert cl.stage == "stage2"
    # 全部有判据的检查通过
    assert cl.all_passed(), [c.metric for c in cl.failures()]
    metrics = {c.metric: c for c in cl.checks}
    # 平衡态三项：V_bi 边界一致性（自洽）、np 恒等（格式自洽）、J=0（自洽）
    vbi = [c for c in cl.checks if "V_bi" in c.metric][0]
    assert vbi.kind == "自洽"
    mass = [c for c in cl.checks if "质量作用" in c.metric][0]
    assert mass.kind == "格式自洽"
    j0 = [c for c in cl.checks if "J = 0" in c.metric][0]
    assert j0.kind == "自洽"
    # 缺陷登记
    assert s2_result["defects"] == ["DEF-S2-01", "DEF-S2-02", "DEF-S2-03"]


# ============================================================================
# 2. S3 V_D = 0 退化 S1 MOS C-V
# ============================================================================

@pytest.fixture(scope="module")
def cv_closure():
    """S3 表格 phi_n=0 列 vs S1 求解器（同一 x 网格）。"""
    x = make_surface_grid(1.2e-6, d0=0.125e-9, h_max=4.0e-9)
    tab = SliceTable(np.linspace(-0.15, 1.6, 120), np.linspace(0.0, 2.0, 11), x)
    tab.build()
    # S1：同一网格按 phi_s 延续求解
    qs1, vg1 = [], []
    prev = None
    for ps in np.linspace(-0.15, 1.6, 120):
        phi, _, _ = solve_poisson(ps, x, C.NA, phi_init=prev)
        prev = phi
        qs, _ = s1_postprocess(phi, x)
        qs1.append(qs)
        vg1.append(C.VFB + ps - qs / C.COX)
    qs1 = np.array(qs1)
    vg1 = np.array(vg1)
    # S3：表格 phi_n=0 列
    ipn0 = int(np.argmin(np.abs(tab.phi_n_grid - 0.0)))
    vg3 = tab.vg[:, ipn0]
    qs3 = tab.qs[:, ipn0]
    m3 = np.isfinite(vg3) & (vg3 > -1.2) & (vg3 < 2.8)
    # 对齐到共同 V_G 网格（仅重叠区）
    vg_common = np.linspace(-1.0, 2.5, 240)
    qs3_i = np.interp(vg_common, vg3[m3], qs3[m3])
    qs1_i = np.interp(vg_common, vg1, qs1)
    dvg = np.abs(qs3_i - qs1_i) / C.COX          # V
    return dict(vg=vg_common, dvg_mV=dvg * 1e3, max_mV=float(dvg.max() * 1e3))


def test_s3_vd0_degenerates_to_s1(cv_closure):
    """V_D=0 时 S3 与 S1 的 Q_s-V_G 一致：等效栅压差 < 1 mV。

    本测试用 120x11 粗表，且表格 phi_s 采样网格与 S1 扫描网格相同
    （linspace(-0.15, 1.6, 120)）——phi_n=0 时两阶段解同一组方程
    （S3 准平衡 Poisson 逐项退化为 S1），逐点一致到机器精度
    （实测 max < 1e-10 mV）。判据取 < 1 mV：覆盖生产配置的插值
    口径对比（见 test_s3_vd0_production_interp_closure：240 点表 vs
    120 点 S1 扫描，全程最大 0.550 mV）。
    """
    assert cv_closure["max_mV"] < 1.0, (
        f"C-V 闭环失效: 最大等效栅压差 {cv_closure['max_mV']:.3f} mV")
    # 中压区（|V_G|<1）应显著更准（表插值误差随 phi_s 灵敏度上升）
    mid = np.abs(cv_closure["vg"]) < 1.0
    assert cv_closure["dvg_mV"][mid].max() < 0.5


@pytest.fixture(scope="module")
def cv_closure_production():
    """生产口径的 C-V 闭环（fig_cv_loop 复算，约 2 s）。

    S3 侧：生产 phi_s 网格（240 点、起点 -0.3）上只建 phi_n=0 列
    的真实 SliceTable（与生产 240x110 表同网格，列间独立故单列等价）；
    S1 侧：fig_cv_loop 同款 120 点扫描（起点 -0.15）。两阶段在不同
    采样网格上独立求解，插值到公共 V_G 网格比较——检验的是完整
    生产链路（含网格差异引入的插值误差），而非逐点退化。
    """
    x = make_surface_grid(1.2e-6, d0=0.125e-9, h_max=4.0e-9)
    tab = SliceTable(np.linspace(-0.3, 1.6, 240), np.linspace(0.0, 2.0, 1), x)
    tab.build()
    vg3, qs3 = tab.vg[:, 0], tab.qs[:, 0]
    qs1, vg1 = [], []
    prev = None
    for ps in np.linspace(-0.15, 1.6, 120):
        phi, _, _ = solve_poisson(ps, x, C.NA, phi_init=prev)
        prev = phi
        qs, _ = s1_postprocess(phi, x)
        qs1.append(qs)
        vg1.append(C.VFB + ps - qs / C.COX)
    qs1, vg1 = np.array(qs1), np.array(vg1)
    m3 = np.isfinite(vg3) & (vg3 > -1.2) & (vg3 < 2.8)
    vg_common = np.linspace(-1.0, 2.5, 240)
    qs3_i = np.interp(vg_common, vg3[m3], qs3[m3])
    qs1_i = np.interp(vg_common, vg1, qs1)
    dvg_mV = np.abs(qs3_i - qs1_i) / C.COX * 1e3
    return dict(vg=vg_common, dvg_mV=dvg_mV,
                max_mV=float(dvg_mV.max()),
                vg_at_max=float(vg_common[np.argmax(dvg_mV)]),
                m3_vg=(float(vg3[m3].min()), float(vg3[m3].max())))


def test_s3_vd0_production_interp_closure(cv_closure_production):
    """生产口径 C-V 闭环（S3 报告 8.10 节 / 图 S3-F10 的数值锁定）。

    有效比较区 V_G in [-1.09, 2.60] V（phi_n=0 列 SKIP_OV 上限
    phi_s=0.9 对应 V_G=2.60 V；深积累角区被 m3 掩码剔除）。
    判据与实测（2024 系列化复算）：
      - 全程（V_G in [-1, 2.5]）最大等效栅压差 < 1 mV
        （实测 0.550 mV @ V_G=1.40 V，强反型起始段斜率最大处）；
      - V_G = 0 / 1 / 2 V 检查点 < 0.01 / 0.15 / 0.5 mV
        （实测 0.002 / 0.087 / 0.390 mV）。
    偏差非单调：两网格插值误差在 0.05~0.55 mV 间波动，不随栅压单调。
    """
    r = cv_closure_production
    # 有效比较区覆盖绘图范围（无平外推）
    assert r["m3_vg"][0] < -1.0 and r["m3_vg"][1] > 2.5, (
        f"有效比较区 {r['m3_vg']} 未覆盖 [-1, 2.5] V")
    assert r["max_mV"] < 1.0, (
        f"C-V 闭环（生产口径）失效: 最大 {r['max_mV']:.3f} mV "
        f"@ V_G={r['vg_at_max']:.2f} V")
    for v, tol in [(0.0, 0.01), (1.0, 0.15), (2.0, 0.5)]:
        i = np.argmin(np.abs(r["vg"] - v))
        assert r["dvg_mV"][i] < tol, (
            f"V_G={v} 检查点超差: {r['dvg_mV'][i]:.3f} mV > {tol} mV")
