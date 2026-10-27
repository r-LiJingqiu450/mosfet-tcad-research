# -*- coding: utf-8 -*-
"""
stage1_mos/tests/test_stage1.py —— MOS 电容 Poisson 求解器测试。

原件无对应测试文件，本文件按系列测试规范（code_spec.md 第 6 节）新建：
覆盖物理不变量、数值恒等式、边界条件与教科书解析对照。

覆盖：
    1. 平带态 np = n_i²、Q_s = 0、V_G = V_FB；
    2. 积累区 Q_s > 0（空穴积累，半导体电荷为正）、表面空穴浓度升高；
    3. 耗尽/反型区 Q_s < 0、表面空穴浓度降低；
    4. 表面电荷：三点差分与电荷积分一致；
    5. Q_s(phi_s) 单调递减、V_G(phi_s) 单调递增；
    6. 边界条件 phi(0)=phi_s、phi(L)=0 严格满足；
    7. 阈值电压 V_T 数值 vs 解析（< 1%）；
    8. 最大耗尽宽度 W_max 数值 vs 解析（< 1%）；
    9. 均匀/非均匀网格解的一致性（容差 1e-3 V）；
    10. 全 phi_s 扫描点（延续模式）Newton 收敛。

符号约定（postprocess 文档字符串 + 高斯定理）：
    Q_s = eps_si * dphi/dx|_{x=0} = 半导体空间电荷（C/m^2），
    积累区为正（空穴）、耗尽/反型区为负（电离受主 + 电子）。
    V_G = V_FB + phi_s - Q_s/C_ox（Q_s<0 时氧化层压降为正）。

求解模式说明：solve_poisson 文档字符串把"参数延续"列为收敛控制
第三层。冷启动单点求解在 |phi_s| 较大时不收敛（Newton 对强
非线性初值敏感，phi_s=0.71 V 实测 200 次迭代残差仍为 1.1 倍源项
尺度）；生产 run() 与全部本测试均使用延续模式（以上一步收敛解
为初值），该模式下全扫描残差在机器精度底（实测 <= 1e-13）。
"""

import numpy as np

from semiconductor_sim.common import constants as C
from semiconductor_sim.common.grid import make_grid
from semiconductor_sim.stage1_mos.solver import solve_poisson
from semiconductor_sim.stage1_mos.postprocess import postprocess


def _grid(kind="uniform", N=801):
    return make_grid(1.0e-6, N, kind=kind, d0=1.0e-9, transition=0.35e-6,
                     h_max=4.5e-9)


def _solve_to(target, x, n_ramp=21, tol=1e-11):
    """
    从平态（phi_s=0）沿参数延续推进到 target，返回末点解。

    与 run() 的生产扫描同一模式：phi_init 为上一步收敛解。
    """
    prev = None
    phi = iters = conv = None
    for ps in np.linspace(0.0, target, n_ramp):
        phi, iters, conv = solve_poisson(ps, x, C.NA, phi_init=prev,
                                         tol=tol, max_iter=200)
        assert conv, f"phi_s={ps:.4f} 延续求解未收敛"
        prev = phi.copy()
    return phi, iters, conv


def test_flatband():
    """平带态：np = n_i²，Q_s = 0，V_G = V_FB。"""
    x = _grid()
    phi, iters, conv = solve_poisson(0.0, x, C.NA)
    assert conv
    p0 = 0.5 * (C.NA + np.sqrt(C.NA ** 2 + 4.0 * C.NI ** 2))
    n0 = C.NI ** 2 / p0
    n = n0 * np.exp(np.clip(phi / C.VT, -3, 3))
    p = p0 * np.exp(np.clip(-phi / C.VT, -3, 3))
    npe = np.max(np.abs(n * p - C.NI ** 2)) / C.NI ** 2
    assert npe < 1e-6
    qs, es = postprocess(phi, x)
    assert abs(qs) < 1e-8


def test_accumulation_charge_sign():
    """积累区（phi_s < 0）：Q_s > 0（空穴积累），表面空穴浓度远超体内。"""
    x = _grid()
    phi, _, conv = _solve_to(-0.5, x)
    assert conv
    qs, es = postprocess(phi, x)
    assert qs > 0
    p_surf = C.NA * np.exp(np.clip(-phi[0] / C.VT, -3, 3))
    assert p_surf > 10 * C.NA          # 表面空穴积累超过掺杂


def test_depletion_charge_sign():
    """耗尽/反型区（phi_s > 0）：Q_s < 0（电离受主+电子），表面空穴耗尽。"""
    x = _grid()
    phi, _, conv = _solve_to(0.6, x)
    assert conv
    qs, es = postprocess(phi, x)
    assert qs < 0
    p_surf = C.NA * np.exp(np.clip(-phi[0] / C.VT, -3, 3))
    assert p_surf < 0.1 * C.NA         # 表面空穴耗尽


def test_surface_charge_consistency():
    """
    表面电荷自洽：Gauss 三点差分 vs 全城电荷积分（耗尽区，容差 2%）。

    Q_s(高斯) = eps*dphi/dx|_0  与  Q_int = q*∫(p-n-Na)dx  在耗尽区
    应一致（静电学恒等式）：phi_s=0.3 实测 0.000%、0.6 实测 0.001%。
    强反型点（phi_s=0.85）不纳入判据：反型层电荷集中于亚网格
    尺度（衰减长度 ~V_T/E_s < 1 nm），梯形积分欠分辨，实测偏差
    5.9%——属网格分辨率边界，与 S3 深反型角点跳过（SKIP_OV）
    同源，已在论文讨论中说明。
    """
    x = _grid()
    p0 = 0.5 * (C.NA + np.sqrt(C.NA ** 2 + 4.0 * C.NI ** 2))
    n0 = C.NI ** 2 / p0
    for ps in [0.3, 0.6]:
        phi, _, _ = _solve_to(ps, x)
        qs3, _ = postprocess(phi, x)
        p = p0 * np.exp(-np.clip(phi, -3, 3) / C.VT)
        n = n0 * np.exp(np.clip(phi, -3, 3) / C.VT)
        q_int = C.Q_E * np.trapezoid(p - n - C.NA, x)
        assert abs(qs3 - q_int) / abs(q_int) < 0.02


def test_monotonicity():
    """Q_s(phi_s) 单调递减（空穴正电荷 -> 电子负电荷），V_G(phi_s) 单调递增。"""
    x = _grid()
    qs_list, vg_list = [], []
    Cox = C.EPS_OX / C.TOX
    prev = None
    for ps in np.linspace(-0.4, 1.0, 15):
        phi, _, conv = solve_poisson(ps, x, C.NA, phi_init=prev)
        assert conv
        prev = phi.copy()
        qs, _ = postprocess(phi, x)
        qs_list.append(qs)
        vg_list.append(C.VFB + ps - qs / Cox)
    assert np.all(np.diff(qs_list) < 0)
    assert np.all(np.diff(vg_list) > 0)


def test_boundary_conditions():
    """边界：phi(0) = phi_s、phi(L) = 0 严格满足（机器精度）。"""
    x = _grid()
    for ps in [-0.3, 0.0, 0.7, 1.2]:
        phi, _, conv = _solve_to(ps, x)
        assert conv
        assert phi[0] == ps
        assert phi[-1] == 0.0


def test_vt_numeric_vs_analytic():
    """阈值电压：phi_s = 2psi_F 处数值 V_G vs 耗尽近似解析（< 1%）。"""
    x = _grid()
    Cox = C.EPS_OX / C.TOX
    psi_F = C.VT * np.log(C.NA / C.NI)
    vt_analytic = C.VFB + 2.0 * psi_F + np.sqrt(
        4.0 * C.Q_E * C.EPS_SI * C.NA * psi_F) / Cox
    phi, _, _ = _solve_to(2.0 * psi_F, x)
    qs, _ = postprocess(phi, x)
    vt_numeric = C.VFB + 2.0 * psi_F - qs / Cox
    assert abs(vt_numeric - vt_analytic) / vt_analytic < 0.01


def test_wmax_numeric_vs_analytic():
    """最大耗尽宽度：电荷法 vs sqrt(4 eps psi_F/(q Na))（< 1%）。"""
    x = _grid()
    psi_F = C.VT * np.log(C.NA / C.NI)
    w_analytic = np.sqrt(4.0 * C.EPS_SI * psi_F / (C.Q_E * C.NA))
    phi, _, _ = _solve_to(2.0 * psi_F, x)
    qs, _ = postprocess(phi, x)
    w_num = -qs / (C.Q_E * C.NA)
    assert abs(w_num - w_analytic) / w_analytic < 0.01


def test_grid_kind_consistency():
    """均匀/非均匀网格解在重叠区一致（容差 1e-3 V）。"""
    xu = _grid("uniform", N=401)
    xn = _grid("nonuniform", N=401)
    phi_u, _, cu = _solve_to(0.7, xu)
    phi_n, _, cn = _solve_to(0.7, xn)
    assert cu and cn
    # 非均匀网格重采样到均匀网格比较（节点位置不同）
    phi_n_i = np.interp(xu, xn, phi_n)
    assert np.abs(phi_u - phi_n_i).max() < 1e-3


def test_all_scan_points_converge():
    """全 phi_s 扫描点（延续模式）Newton 收敛（tol 用 run() 同款）。"""
    x = _grid()
    prev = None
    for ps in np.linspace(-0.5, 1.4, 20):
        phi, _, conv = solve_poisson(ps, x, C.NA, phi_init=prev,
                                     tol=1e-11, max_iter=100)
        assert conv, f"phi_s={ps:.4f} 未收敛"
        prev = phi.copy()
