# -*- coding: utf-8 -*-
"""
common/ 七模块单元测试（tests/test_common.py）。

覆盖：
    1. constants  派生量与解析函数（VT、psi_F、gamma、V_T、Square-law 极限）
    2. grid       四类网格的结构性质（单调/端点/表面段 h0=h1/对称/参数校验）
    3. thomas     标量三对角 vs 稠密解；带状入口与标量一致；块 Thomas vs 稠密
    4. bernoulli  B(0)=1、B'e 恒等式、大值极限、导数 vs 数值微分
    5. qfp        列变换链式法则（跨 14 个数量级载流子浓度）
    6. validation Check/Checklist 四分类记录与表格行生成
    7. io         JSON/NPZ roundtrip、golden 比对、TerminalLog 双写
"""

import os

import numpy as np
import pytest

from semiconductor_sim.common import constants as C
from semiconductor_sim.common.grid import (make_grid, make_surface_grid,
                                           make_junction_grid,
                                           make_channel_grid, grid_stats)
from semiconductor_sim.common.thomas import (thomas_solve, thomas_solve_banded,
                                             block_thomas)
from semiconductor_sim.common.bernoulli import bernoulli, bernoulli_deriv
from semiconductor_sim.common.qfp import qfp_columns
from semiconductor_sim.common.validation import (Check, Checklist, VALIDATION,
                                                 SELF_CHECK, FORMAT_CHECK,
                                                 ARTIFACT, rel_dev,
                                                 max_rel_dev, format_sci,
                                                 format_pct)
from semiconductor_sim.common.io import (save_json, load_json, compare_golden,
                                         save_npz, load_npz, TerminalLog,
                                         ensure_results_dir, GOLDEN)


# ============================================================================
# 1. constants
# ============================================================================

def test_thermal_voltage():
    """V_T = kT/q = 25.852 mV @300K。"""
    assert abs(C.VT - 0.025852) < 1e-6


def test_fermi_potential():
    """psi_F = V_T ln(Na/n_i)。"""
    assert abs(C.PSI_F - C.VT * np.log(C.NA / C.NI)) < 1e-12
    assert abs(C.PSI_F - 0.3572) < 1e-3


def test_body_factor_and_vt_long():
    """gamma 与长沟道阈值：派生量 == analytic_vt() 独立计算。"""
    gamma = np.sqrt(2 * C.EPS_SI * C.Q_E * C.NA) / C.COX
    assert abs(C.GAMMA - gamma) < 1e-12
    vt = C.VFB + 2 * C.PSI_F + C.GAMMA * np.sqrt(2 * C.PSI_F)
    assert abs(C.VT_LONG - vt) < 1e-12
    assert abs(C.analytic_vt() - C.VT_LONG) < 1e-12
    assert abs(C.VT_LONG - 0.9963) < 1e-3      # 论文引用值


def test_square_law_limits():
    """Square-law：亚阈为 0；线性区/饱和区公式对拍。"""
    vg = np.array([0.5, 1.5, 2.0])
    vt = C.analytic_vt()
    # 亚阈（V_G < V_T）必须为 0
    assert C.square_law_id(0.5, 0.05) == 0.0
    # 线性区（V_D < V_G - V_T）
    id_lin = C.square_law_id(1.5, 0.05)
    ov = 1.5 - vt
    ref = C.MU_N * C.COX * (C.W_CH / C.L_CH) * (ov * 0.05 - 0.5 * 0.05 ** 2)
    assert abs(id_lin - ref) / ref < 1e-12
    # 饱和区（V_D > V_G - V_T）
    id_sat = C.square_law_id(1.5, 1.5)
    ref_sat = 0.5 * C.MU_N * C.COX * (C.W_CH / C.L_CH) * ov ** 2
    assert abs(id_sat - ref_sat) / ref_sat < 1e-12
    # 连续性：V_D = V_G - V_T 处两段相等
    id_edge = C.square_law_id(1.5, ov)
    assert abs(id_edge - id_sat) / id_sat < 1e-12


def test_vbi_and_depletion_analytic():
    """S2 解析：V_bi 与耗尽宽度（对拍手算公式）。"""
    vbi = C.V_bi_analytic(C.NA, C.ND)
    assert abs(vbi - C.VT * np.log(C.NA * C.ND / C.NI ** 2)) < 1e-15
    w0 = C.W_analytic(C.NA, C.ND, 0.0)
    ref = np.sqrt(2 * C.EPS_SI / C.Q_E * (1 / C.NA + 1 / C.ND) * vbi)
    assert abs(w0 - ref) / ref < 1e-12
    # 反偏展宽、正偏收窄
    assert C.W_analytic(C.NA, C.ND, -1.0) > w0
    assert C.W_analytic(C.NA, C.ND, +0.5) < w0


def test_equilibrium_carriers_mass_action():
    """平衡载流子满足质量作用 np = n_i²（机器精度）。"""
    p_p, n_p, p_n, n_n = C.equilibrium_carriers()
    assert abs(n_p * p_p / C.NI ** 2 - 1.0) < 1e-12
    assert abs(p_n * n_n / C.NI ** 2 - 1.0) < 1e-12
    assert abs(p_p - C.NA) / C.NA < 1e-2      # p 侧多子 ≈ Na


# ============================================================================
# 2. grid
# ============================================================================

def test_surface_grid_structure():
    """S3 纵向网格：端点、表面均匀段 h0=h1（三点差分要求）、单调、h_max。"""
    d0, h_max = 0.125e-9, 4.0e-9
    x = make_surface_grid(1.2e-6, d0=d0, h_max=h_max)
    h = np.diff(x)
    assert x[0] == 0.0 and abs(x[-1] - 1.2e-6) < 1e-18
    assert np.all(np.diff(x) > 0)
    assert abs(h[0] - d0) < 1e-18 and abs(h[1] - d0) < 1e-18
    assert h.max() <= h_max * (1 + 1e-12)
    # 表面段覆盖 transition=0.08 um
    n_surf = int(round(0.08e-6 / d0))
    assert abs(h[:n_surf] - d0).max() < 1e-18


def test_surface_grid_rejects_bad_params():
    with pytest.raises(ValueError):
        make_surface_grid(1.0e-6, d0=-1.0)
    with pytest.raises(ValueError):
        make_surface_grid(1.0e-6, transition=2.0e-6)   # transition > x_max
    with pytest.raises(ValueError):
        make_surface_grid(1.0e-6, d0=10e-9, h_max=1e-9)  # h_max < d0


def test_junction_grid_symmetry():
    """S2 结网格：x=0 为节点、左右镜像、单调。"""
    x = make_junction_grid(5.0e-6)
    assert abs(x[np.argmin(np.abs(x))]) < 1e-18      # 0 恰为节点
    assert abs(x[0] + 2.5e-6) < 1e-18 and abs(x[-1] - 2.5e-6) < 1e-18
    assert np.all(np.diff(x) > 0)
    # 镜像对称（节点数奇偶均容忍 1e-18 舍入）
    xr = x[x >= 0]
    xl = -x[x <= 0][::-1]
    n = min(len(xr), len(xl))
    assert np.abs(xr[:n] - xl[:n]).max() < 1e-18


def test_channel_grid():
    y = make_channel_grid(1.0e-6, Ny=41)
    assert len(y) == 41 and y[0] == 0.0 and y[-1] == 1.0e-6
    with pytest.raises(ValueError):
        make_channel_grid(1.0e-6, Ny=2)


def test_make_grid_uniform_and_nonuniform():
    xu = make_grid(1.0e-6, 801, kind="uniform")
    assert len(xu) == 801 and xu[0] == 0.0 and xu[-1] == 1.0e-6
    assert np.allclose(np.diff(xu), 1.0e-6 / 800)
    xn = make_grid(1.0e-6, 801, kind="nonuniform", d0=1.0e-9)
    assert xn[0] == 0.0 and abs(xn[-1] - 1.0e-6) < 1e-18
    assert np.all(np.diff(xn) > 0)
    # 表面密、体内疏
    assert np.diff(xn)[0] < np.diff(xn)[-1]


def test_grid_stats():
    x = make_surface_grid(1.2e-6, d0=0.125e-9)
    st = grid_stats(x)
    assert st["N"] == len(x) and st["L"] == x[-1]
    assert st["h0"] == np.diff(x)[0]
    assert st["ratio"] == st["h_max"] / st["h0"]


# ============================================================================
# 3. thomas
# ============================================================================

def test_thomas_vs_dense():
    """标量 Thomas：随机严格对角占优系统 vs 稠密解。"""
    rng = np.random.default_rng(7)
    n = 300
    a = rng.normal(size=n) * 0.1
    b = rng.normal(size=n) + 5.0          # 对角占优
    c = rng.normal(size=n) * 0.1
    d = rng.normal(size=n)
    a[0] = 0.0
    c[-1] = 0.0
    A = np.diag(b) + np.diag(c[:-1], 1) + np.diag(a[1:], -1)
    x_ref = np.linalg.solve(A, d)
    assert np.abs(thomas_solve(a, b, c, d) - x_ref).max() < 1e-10


def test_thomas_banded_matches_scalar():
    """带状入口（scipy 优先）与标量追赶法逐位一致。"""
    rng = np.random.default_rng(11)
    n = 200
    a = rng.normal(size=n) * 0.1
    b = rng.normal(size=n) + 5.0
    c = rng.normal(size=n) * 0.1
    d = rng.normal(size=n)
    a[0] = 0.0
    c[-1] = 0.0
    x1 = thomas_solve(a, b, c, d)
    x2 = thomas_solve_banded(a, b, c, d)
    assert np.abs(x1 - x2).max() < 1e-12


def test_block_thomas_vs_dense():
    """块 Thomas（3x3）：随机/对角占优系统 vs 稠密解（迁移自 S2 测试）。"""
    rng = np.random.default_rng(42)
    M = 200
    x_true = rng.normal(size=(M, 3))
    results = []
    # 1) 随机块三对角（对角占优保证非奇异）
    Ad = rng.normal(size=(M, 3, 3)) + 3 * np.eye(3)
    Al = rng.normal(size=(M, 3, 3)) * 0.3
    Au = rng.normal(size=(M, 3, 3)) * 0.3
    Al[0] = 0
    Au[-1] = 0
    J = np.zeros((3 * M, 3 * M))
    for j in range(M):
        J[3 * j:3 * j + 3, 3 * j:3 * j + 3] = Ad[j]
        if j > 0:
            J[3 * j:3 * j + 3, 3 * (j - 1):3 * (j - 1) + 3] = Al[j]
        if j < M - 1:
            J[3 * j:3 * j + 3, 3 * (j + 1):3 * (j + 1) + 3] = Au[j]
    rhs = (J @ x_true.ravel()).reshape(M, 3)
    results.append(("随机块三对角",
                    np.abs(block_thomas(Ad, Al, Au, rhs) - x_true).max()))
    # 2) 对角占优系统（接近实际结构）
    Ad2 = np.zeros((M, 3, 3))
    Ad2[:, 0, 0] = 4
    Ad2[:, 1, 1] = 5
    Ad2[:, 2, 2] = 6
    Al2 = np.zeros((M, 3, 3))
    Al2[:, 0, 0] = -1
    Al2[:, 1, 1] = -2
    Al2[:, 2, 2] = -3
    Au2 = np.zeros((M, 3, 3))
    Au2[:, 0, 0] = -1
    Au2[:, 1, 1] = -2
    Au2[:, 2, 2] = -3
    Al2[0] = 0
    Au2[-1] = 0
    J2 = np.zeros((3 * M, 3 * M))
    for j in range(M):
        J2[3 * j:3 * j + 3, 3 * j:3 * j + 3] = Ad2[j]
        if j > 0:
            J2[3 * j:3 * j + 3, 3 * (j - 1):3 * (j - 1) + 3] = Al2[j]
        if j < M - 1:
            J2[3 * j:3 * j + 3, 3 * (j + 1):3 * (j + 1) + 3] = Au2[j]
    rhs2 = (J2 @ x_true.ravel()).reshape(M, 3)
    results.append(("对角占优",
                    np.abs(block_thomas(Ad2, Al2, Au2, rhs2) - x_true).max()))
    for name, err in results:
        assert err < 1e-10, f"{name}: 块 Thomas 与稠密解不一致 {err:.3e}"


# ============================================================================
# 4. bernoulli
# ============================================================================

def test_bernoulli_values_and_limits():
    """B(0)=1；大值极限 B(+inf)=x、B(-inf)->0；恒等式 B(x)e^{-x}=B(-x)。"""
    assert bernoulli(0.0) == 1.0
    assert abs(bernoulli(800.0) - 800.0) < 1e-12          # e^{-x} 下溢
    assert abs(bernoulli(-800.0)) < 1e-12                  # 负极限 ~0
    xs = np.array([-30.0, -2.0, -0.1, 0.0, 0.1, 2.0, 30.0])
    B = bernoulli(xs)
    Bm = bernoulli(-xs)
    assert np.abs(B * np.exp(-xs) - Bm).max() < 1e-12 * max(1.0, np.abs(B).max())


def test_bernoulli_small_taylor():
    """
    小 |x| 分支精度（阈值 1e-8，与原件 pn_junction.py 一致）：

        |x| < 1e-8 : 泰勒 B = 1 + x/2（x^2/12 项在双精度下不可见，
                    ulp(1)=2.2e-16 > 8.3e-18），与 1+x/2 逐位一致；
        其余区间   : 直接式 x/(1-e^{-x})，与 expm1 高精度参照在
                    x >= 1e-4 时相对误差 <= 3e-14（实测）。

    已知边界（不修复，如实记录）：1e-8 < |x| < ~1e-4 时间接式
    1-e^{-x} 存在抵消 cancellation，绝对误差 ~ eps/|x|（阈值处实测
    2.2e-9）。求解器工作点的相邻节点 x 或低于阈值（泰勒精确）、
    或高于 1e-4（误差 <3e-14）；该 gap 不触及收敛判据 1e-9 V。
    本测试锁定上述分区行为，包括阈值处 2e-9 量级的已知跳变。
    """
    # 泰勒分支：与 1+x/2 逐位一致（x^2/12 不可见）
    for x in (1e-12, 1e-10, 1e-9):
        assert float(bernoulli(x)) == float(1.0 + x / 2.0)
    # 阈值边界：x=1e-8 起进入直接式，跳变量 ~2e-9（已知 cancellation）
    assert abs(float(bernoulli(1e-8)) - float(bernoulli(9.999e-9))) < 1e-8
    # 直接式区间：与 expm1 参照（无抵消）比对，x >= 1e-4 相对误差 <= 1e-13
    for x in (1e-4, 1e-3, 0.1, 0.5, 2.0, 10.0, -0.3, -2.0, -10.0):
        ref = x / (-np.expm1(-x))
        assert abs(float(bernoulli(x)) - ref) <= 1e-13 * abs(ref)
    # 大 |x| 极限
    assert float(bernoulli(800.0)) == 800.0
    assert float(bernoulli(-800.0)) == 0.0


def test_bernoulli_deriv_vs_fd():
    """B'(x) vs 中心有限差分（中等 |x| 区间，避开泰勒/极限分支）。"""
    xs = np.array([-5.0, -1.0, -0.2, 0.2, 1.0, 5.0])
    h = 1e-6
    fd = (bernoulli(xs + h) - bernoulli(xs - h)) / (2 * h)
    an = bernoulli_deriv(xs)
    assert np.abs(fd - an).max() < 1e-6
    # 端点极限
    assert abs(bernoulli_deriv(0.0) - 0.5) < 1e-12
    assert abs(bernoulli_deriv(800.0) - 1.0) < 1e-12
    assert abs(bernoulli_deriv(-800.0)) < 1e-12


# ============================================================================
# 5. qfp
# ============================================================================

def test_qfp_columns_chain_rule():
    """qfp_columns(A, cn, cp) == A @ du/dv（链式法则代数单测，迁移自 S2）。"""
    rng = np.random.default_rng(42)
    A = rng.normal(size=(50, 3, 3)) * rng.choice(
        [1e-20, 1.0, 1e10], size=(50, 1, 1))
    cn = 10.0 ** rng.uniform(9, 23, size=50)     # n/V_T
    cp = 10.0 ** rng.uniform(9, 23, size=50)     # p/V_T
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
# 6. validation
# ============================================================================

def test_checklist_records_and_rows():
    """Checklist：四分类记录、all_passed 语义、to_rows 表行。"""
    cl = Checklist("test")
    cl.add(metric="a", kind=VALIDATION, reference="ref", criterion="< 1%",
           measured="0.1%", passed=True)
    cl.add(metric="b", kind=SELF_CHECK, reference="ref", criterion="自洽",
           measured="ok", passed=True)
    cl.add(metric="c", kind=ARTIFACT, reference="ref", criterion="无",
           artifact="ART-T-01")
    cl.add(metric="d", kind=VALIDATION, reference="ref", criterion="< 5%",
           measured="9%", passed=False)
    assert len(cl.checks) == 4
    assert not cl.all_passed()                    # 有 False
    assert len(cl.failures()) == 1
    rows = cl.to_rows()
    assert len(rows) == 4
    assert rows[0][-1].startswith("通过")
    assert rows[3][-1].startswith("未达标")
    assert rows[2][-1] == "—"                     # artifact 无判据
    # 移除失败项后全部通过
    cl.checks = [c for c in cl.checks if c.passed is not False]
    assert cl.all_passed()


def test_checklist_numpy_bool_normalization():
    """
    DEF-C-02 回归：numpy 布尔判据的归一化。

    调用方的判据多是 numpy 标量比较结果（np.bool_）。若原样存入，
    `is True/is False` 身份判断会同时逃过二者：失败项被 all_passed/
    failures 漏判、通过项在 to_rows 显示"—"。add() 必须归一化，
    三个判断方法亦须对 np.bool_ 稳健。
    """
    cl = Checklist("test")
    cl.add(metric="np-pass", kind=VALIDATION, reference="ref",
           criterion="< 6%", measured="4.6%", passed=np.bool_(True))
    cl.add(metric="np-fail", kind=VALIDATION, reference="ref",
           criterion="< 2%", measured="4.6%", passed=np.bool_(False))
    cl.add(metric="py-none", kind=ARTIFACT, reference="ref",
           criterion="无", artifact="ART-T-02")
    assert isinstance(cl.checks[0].passed, bool)
    assert not cl.all_passed()
    assert [c.metric for c in cl.failures()] == ["np-fail"]
    rows = cl.to_rows()
    assert rows[0][-1].startswith("通过")
    assert rows[1][-1].startswith("未达标")
    assert rows[2][-1] == "—"
    # 直接构造 Check（绕过 add）时判断方法同样稳健
    cl2 = Checklist("test")
    cl2.checks.append(Check(metric="raw", kind=VALIDATION, reference="r",
                            criterion="< 1%", measured="2%",
                            passed=np.bool_(False)))
    assert not cl2.all_passed()
    assert len(cl2.failures()) == 1
    assert cl2.to_rows()[0][-1].startswith("未达标")


def test_validation_helpers():
    """rel_dev / max_rel_dev / format_* 数值工具。"""
    a = np.array([1.0, 2.0, 3.0])
    b = np.array([1.0, 2.5, 3.0])
    assert np.allclose(rel_dev(a, b), [0.0, -0.2, 0.0])
    m, i = max_rel_dev(a, b)
    assert abs(m - 0.2) < 1e-15 and i == 1
    assert format_sci(1.234e-5) == "1.234e-05"
    assert format_sci(None) == "—"
    assert format_pct(0.0123) == "1.23%"


# ============================================================================
# 7. io
# ============================================================================

def test_json_roundtrip_and_golden(tmp_path):
    """JSON 存取 roundtrip；compare_golden 通过/失败/缺失跳过。"""
    data = {"vt": 1.0422, "ss": 77.8, "arr": [1, 2, 3]}
    p = save_json("_test_common_json", data, dirpath=str(tmp_path))
    assert load_json("_test_common_json", dirpath=str(tmp_path)) == data
    cur = {"vt": 1.0422 * (1 + 1e-9), "ss": 77.8}
    res = compare_golden("_test_common_json", cur, rtol=1e-6,
                         dirpath=str(tmp_path))
    assert res["vt"][0] and res["ss"][0]           # 均在容差内
    cur2 = {"vt": 1.10}
    res2 = compare_golden("_test_common_json", cur2, rtol=1e-6,
                          dirpath=str(tmp_path))
    assert not res2["vt"][0]                       # 超容差判失败
    # golden 缺失 -> None（跳过）
    assert compare_golden("_no_such_golden", cur, dirpath=str(tmp_path)) is None


def test_npz_roundtrip(tmp_path):
    x = np.linspace(0, 1, 11)
    save_npz("_test_common_npz", x=x, phi=np.ones(5))
    z = load_npz("_test_common_npz")
    assert np.array_equal(z["x"], x)
    assert np.array_equal(z["phi"], np.ones(5))


def test_terminal_log_dual_write(tmp_path):
    """TerminalLog：屏幕回显 + 文件留存一致。"""
    d = ensure_results_dir("_test_stage")
    log = TerminalLog("_test_stage", echo=False)
    print("hello series")
    log.close()
    txt = open(d + "/terminal.txt", encoding="utf-8").read()
    assert "hello series" in txt


def test_ensure_results_dir(tmp_path):
    d = ensure_results_dir("_test_stage2")
    assert d.endswith("_test_stage2")
    assert os.path.isdir(d)
