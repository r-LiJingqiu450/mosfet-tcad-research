# -*- coding: utf-8 -*-
"""
stage3_mosfet/tests/test_stage3.py —— MOSFET GCA 单元测试。

由原件 mosfet/tests/test_mosfet.py 迁移（导入改为系列包路径，
_square_law_vt 已公开为 square_law_id_vt）。

覆盖：
    1. 切片平带态 np = n_i²、V_G=0；
    2. 表面电荷三点差分与电荷积分一致；
    3. 栅极 phi_s(V_G) 单调递增；
    4. GCA 线性区还原 Square-law（分区间严格判据 + 符号判据）；
    5. 亚阈值摆幅 SS（固定窗口法）+ 局部 SS 物理下限 + 衬偏效应理论对照；
    6. 弱反型不施加 CLM（亚阈饱和由指数拖尾自然给出）；
    7. 电流连续性：线性区及饱和区独立差分检查 < 1e-3；
    8. 深反型角点跳过不污染工作区（V_G<=2 全内插、无钳位）。
"""

import warnings

import numpy as np
import pytest

from semiconductor_sim.common import constants as C
from semiconductor_sim.common.grid import make_surface_grid
from semiconductor_sim.stage3_mosfet.solver import (
    solve_slice, surface_charge, electron_sheet_density, total_charge_density)
from semiconductor_sim.stage3_mosfet.table import gate_voltage, SliceTable
from semiconductor_sim.stage3_mosfet.gca import gca_solve, id_vg, id_vd
from semiconductor_sim.stage3_mosfet import postprocess as pp
from semiconductor_sim.stage3_mosfet.run import square_law_id_vt



X = make_surface_grid(1.2e-6, d0=0.125e-9, h_max=4.0e-9)


@pytest.fixture(scope="module")
def table():
    from semiconductor_sim.stage3_mosfet.phase1 import production_table
    return production_table()[0]


@pytest.fixture(scope="module")
def x():
    return X


def test_flatband_equilibrium():
    """平带态：np = n_i²，V_G = 0，Q_s = 0。"""
    res = solve_slice(0.0, 0.0, X)
    npe = np.max(np.abs(res["n"] * res["p"] - C.NI ** 2)) / C.NI ** 2
    assert npe < 1e-10
    qs, _ = surface_charge(res, X)
    assert abs(qs) < 1e-10
    assert abs(gate_voltage(0.0, qs)) < 1e-10


def test_surface_charge_consistency():
    """表面电荷：三点差分与电荷积分一致（弱反型/耗尽区）。"""
    for ps in [0.3, 0.6, 0.85]:
        res = solve_slice(ps, 0.0, X)
        qs3, _ = surface_charge(res, X)
        qs_int = total_charge_density(res, X)
        assert abs(qs3 - qs_int) / abs(qs_int) < 0.02


def test_phi_s_monotonic(table):
    """phi_s(V_G) 单调递增（V_D=0 退化 MOS 电容）。"""
    vg = np.linspace(0.0, 2.0, 41)
    ps = np.array([table.solve_vg(v, 0.0)[0] for v in vg])
    assert np.all(np.diff(ps) >= 0)


def test_gca_linear_square_law(table):
    """
    GCA 线性区与 Square-law（提取 V_T）偏差：分区间严格判据。

    深强反型（V_G >= 1.75）：|偏差| <= 18%（V_G=2 处实测 13.6%）；
    近阈强反型（1.3 <= V_G < 1.75）：|偏差| <= 30%（软开启拖尾，
    V_G≈1.25 处最大 31.4% 出自该区间下沿，不在本判据内）。
    符号判据：GCA 系统性地低于方定律（电荷片物理，见论文 8.6 节）。
    """
    vg_all = np.linspace(0, 2.0, 81)
    id_tmp = id_vg(table, vg_all, 0.05)
    vt, _ = pp.extract_vt(vg_all, id_tmp)
    for lo, hi, bound in [(1.75, 2.0, 18.0), (1.3, 1.75, 30.0)]:
        vg = np.linspace(lo, hi, 11)
        id_gca = id_vg(table, vg, 0.05)
        sq_fit = square_law_id_vt(vg, 0.05, vt)
        dev = (id_gca - sq_fit) / sq_fit * 100
        assert np.all(dev < 0), "GCA 应系统性低于方定律"
        assert np.max(np.abs(dev)) < bound


def test_subthreshold_slope(table):
    """亚阈值摆幅：固定窗口法 65~90 mV/dec；局部最小值与衬偏效应理论一致。"""
    vg = np.linspace(0.0, 2.0, 121)
    id_ = id_vg(table, vg, 0.05)
    ss = pp.extract_ss(vg, id_)
    assert 65.0 < ss < 90.0
    # 局部 SS 曲线：物理下限 ln10*V_T = 59.6 mV/dec 不得被突破
    #（CLM 虚假施加曾导致 V_G≈0.9 处伪 50 mV/dec，回归防护）
    ss_loc = pp.local_ss(vg, id_)
    sub = (vg > 0.5) & (vg < 1.3)
    assert np.nanmin(ss_loc[sub]) > 60.0
    # 近阈值局部最小 SS 与理论 ln10*VT*(1+Cdep(W_max)/Cox) 相差 < 6 mV/dec
    ss_min = np.nanmin(ss_loc[sub])
    ss_th_min = np.log(10) * C.VT * (
        1 + C.EPS_SI / (np.sqrt(2 * C.EPS_SI * C.TWICE_PSI_F /
                               (C.Q_E * C.NA)) * C.COX)) * 1000
    assert abs(ss_min - ss_th_min) < 6.0


def test_clm_not_in_weak_inversion(table):
    """弱反型（V_G=0.9, V_D=0.05）不得被判饱和、不得施加 CLM。"""
    r = gca_solve(0.9, 0.05, table)
    assert not r["saturated"]
    assert r["y_p"] == C.L_CH


def test_continuity_linear(table):
    """电流连续性：线性区 I_D(y) 波动 < 1e-3（验收指标）。"""
    for vg in [1.5, 2.0]:
        r = gca_solve(vg, 0.05, table)
        rel, _, _ = pp.continuity_check(r)
        assert rel < 1e-3


def test_continuity_saturation_finite_difference(table):
    """Independent point-gradient truncation error must pass the original 1e-3 target."""
    r = gca_solve(2.0, 1.5, table, y=np.linspace(0,C.L_CH,641))
    rel, _, _ = pp.continuity_check(r)
    assert rel < 1e-3


def test_saturation_plateau(table):
    """饱和区：V_D 增大时 I_D 变化缓慢（平台，含一阶 CLM）。"""
    ids, vdsats, _ = id_vd(table, 2.0, np.linspace(1.5, 2.0, 11), clm=True)
    # 深饱和段（V_D > V_Dsat）内 I_D 相对变化 < 20%（一阶 CLM）
    rel = (ids[-1] - ids[0]) / ids[0]
    assert 0.0 < rel < 0.20


def test_pinch_off_moves_to_source(table):
    """夹断点随 V_D 增大向源端移动。"""
    _, _, yps = id_vd(table, 2.0, np.array([1.0, 1.5, 2.0]), clm=True)
    assert yps[-1] < yps[0]     # 向源端（y 减小）


def test_table_column_vs_direct_solve(table, x):
    """表格 phi_n=0 列与同网格直接求解一致（建表自洽性）。"""
    ipn0 = int(np.argmin(np.abs(table.phi_n_grid - 0.0)))
    for ips in [40, 80, 120, 150]:
        ps = table.phi_s_grid[ips]
        if not np.isfinite(table.vg[ips, ipn0]):
            continue
        res = solve_slice(ps, 0.0, x)
        qs, _ = surface_charge(res, x)
        vg = C.VFB + ps - qs / C.COX
        ns = electron_sheet_density(res, x)
        assert abs(vg - table.vg[ips, ipn0]) < 1e-6 * max(1.0, abs(vg))
        assert abs(qs - table.qs[ips, ipn0]) < 1e-6 * max(1e-4, abs(qs))
        assert abs(ns - table.ns[ips, ipn0]) < 1e-4 * max(1e10, ns)


def test_skip_corner_does_not_clamp_workspace(table):
    """
    深反型角点（phi_s-phi_n>0.9）跳过后，V_G<=2 V 工作区内所有
    phi_n 列必须为内插（有效 V_G 上限 > 2 V），不得发生钳位外推。
    """
    for ipn in range(table.npn):
        colv = np.isfinite(table.vg[:, ipn])
        assert table.vg[colv, ipn].max() > 2.0


@pytest.mark.parametrize('vd,ith,branch',[(.05,1e-7,'saturation'),(.05,1e-5,'linear')])
def test_threshold_inverse_correct_branch(vd,ith,branch):
    from semiconductor_sim.stage3_mosfet.reference import square_law_constant_current
    r=square_law_constant_current(vd,ith)
    assert r['branch']==branch
    assert np.isclose(square_law_id_vt(r['vg'],vd,C.VT_LONG),ith,rtol=1e-12)


def test_clm_disabled_profile_finite_at_high_drain(table):
    from semiconductor_sim.stage3_mosfet.legacy_gca import gca_solve as old
    with np.errstate(all='raise'):
        for solve in [gca_solve,old]:
            r=solve(2,1.5,table,clm=False)
            assert all(np.all(np.isfinite(r[k])) for k in ['phi_n','ns','phi_s','id_local'])
            assert r['phi_n'][-1] <= 1.5


@pytest.mark.parametrize('vg,vd',[(.5,.05),(1,.05),(1.5,1.5),(2,1.5)])
def test_independent_continuum_reference(table,vg,vd):
    from semiconductor_sim.stage3_mosfet.reference import reference_current
    # The fixture uses the same production table as study and golden.
    ref=reference_current(vg,vd)['ID']
    assert abs(gca_solve(vg,vd,table)['ID']/ref-1)<.01


def test_charge_budget_is_classical_identity(table):
    from semiconductor_sim.stage3_mosfet.reference import charge_budget
    ps,ns,qs=table.solve_vg(2,0)
    r=charge_budget(2,ps,ns,qs)
    assert abs(r['closure_V'])<1e-12
    assert .14<r['surface_term_V']<.17


def test_conservative_edges_and_independent_gradient_are_distinct(table):
    from semiconductor_sim.stage3_mosfet.postprocess import edge_continuity_check
    results=[]
    for n in [161,321,641]:
        r=gca_solve(2,1.5,table,y=np.linspace(0,C.L_CH,n))
        assert edge_continuity_check(r)<1e-9
        results.append(pp.continuity_check(r)[0])
    assert results[-1]<results[0]/3
    assert results[-1]<1e-3


def test_output_grid_does_not_change_terminal_current(table):
    currents=[gca_solve(2,1.5,table,y=np.linspace(0,C.L_CH,n))['ID']
              for n in [81,161,641]]
    assert max(currents)==min(currents)


def test_fine_quadrature_knots_leave_exact_interpolant_integral_unchanged(table):
    a=gca_solve(1.5,1.5,table)['ID']
    b=gca_solve(1.5,1.5,table,nphi=1001)['ID']
    assert abs(a/b-1)<1e-12


@pytest.mark.parametrize('vg,pn',[(5.,0.),(1.,-.1),(1.,2.1)])
def test_table_outside_domain_raises(table,vg,pn):
    with pytest.raises(ValueError):table.solve_vg(vg,pn)


def test_nonconvergence_is_not_a_valid_table_entry(monkeypatch):
    import semiconductor_sim.stage3_mosfet.table as module
    monkeypatch.setattr(module,'solve_slice',lambda *a,**kw:
                        dict(phi=np.zeros_like(X),converged=False,status='injected_failure'))
    with pytest.raises(RuntimeError,match='Invalid slice'):
        module.SliceTable([0,.1],[0,.1],X).build()


def test_small_iteration_budget_cannot_claim_convergence():
    from semiconductor_sim.stage3_mosfet.solver import solve_poisson_slice
    r=solve_poisson_slice(.85,0,X,max_iter=1,warn=False)
    assert not r['converged']
    assert r['correction_inf']>1e-9


@pytest.mark.parametrize('mode',['linear','log'])
def test_transport_against_closed_form_sheets(mode):
    class AnalyticTable:
        phi_n_grid=np.linspace(0,1,9)
        sheet_interpolation=mode
        def ns_vs_phi_n(self,vg,pn):
            return (1e16*np.exp(-4*pn) if mode=='log' else 1e16*(2-pn))
        def solve_vg(self,vg,pn):
            return (vg,float(self.ns_vs_phi_n(vg,np.array([pn]))[0]),0.)
    t=AnalyticTable();r=gca_solve(1,1,t,y=np.linspace(0,C.L_CH,641))
    integral=1e16*((1-np.exp(-4))/4 if mode=='log' else 1.5)
    expected=C.Q_E*C.W_CH*C.MU_N/C.L_CH*integral
    assert abs(r['ID']/expected-1)<1e-13
    yf=r['y']/C.L_CH
    expected_pn=(-np.log1p(-yf*(1-np.exp(-4)))/4 if mode=='log'
                 else 2-np.sqrt(4-3*yf))
    assert np.max(abs(r['phi_n']-expected_pn))<1e-12


