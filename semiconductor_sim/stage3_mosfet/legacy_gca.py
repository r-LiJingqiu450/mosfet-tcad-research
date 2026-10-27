# -*- coding: utf-8 -*-
"""
渐变沟道近似（GCA）求解器与偏压扫描 —— 沿沟道积分电流。

系列化说明：本模块由原 gca_solver.py 与 sweep.py 合并（扫描函数
是 gca_solve 的直接包装，无独立实现）；物理常数取自
common.constants。

核心公式（精确 GCA，漂移+扩散合并后的准平衡形式）：

    沟道电子电流密度 J_n = -q mu_n n d(phi_n)/dy
    面电流  I_D = q W mu_n N_s(y) d(phi_n)/dy   （与 y 无关）

分离变量并从源(y=0, phi_n=0)积分到漏(y=L, phi_n=V_D)：

    I_D = (q W mu_n / L) * ∫_0^{V_D} N_s(phi_n) d(phi_n)

N_s(phi_n) 由 gate_boundary.GateTable 在固定 V_G 下双线性插值得到。
该积分形式：
    1) 线性区自动还原 Square-law（N_s ≈ C_ox(V_G-V_T-phi_n)/q 线性衰减）；
    2) 饱和平台自然涌现（phi_n > V_G-V_T 后 N_s 指数拖尾，积分收敛）；
    3) 亚阈值斜率自动含衬偏效应因子（N_s ∝ exp((phi_s-phi_n)/V_T)）。

饱和区沟道长度调制（CLM，一阶）：
    夹断区长度 Delta_L = sqrt(2 eps_si (V_D - V_Dsat)/(q N_a))，
    I_D = I_Dsat / (1 - Delta_L/L)，夹断点 y_p = L - Delta_L 随 V_D 增大
    向源端移动。

    CLM 仅当源端处于强反型（N_s(0) >= N_MIN_STRONG）时施加。弱反型/
    深亚阈区不存在"沟道夹断"：N_s(phi_n) 的指数拖尾使积分
    ∫_0^{V_D} N_s dphi_n 在 V_D > 若干 V_T 后自然收敛，给出平滑的
    亚阈饱和平台；此时若仍套用耗尽近似的 CLM，会在 I_D-V_G 曲线上
    产生虚假的"饱和"扭结（开发期实测：V_G≈0.9、V_D=0.05 时局部
    SS 被压低到 50 mV/dec，低于 ln10·V_T 物理下限）。
"""

import numpy as np

from common import constants as C

N_MIN_STRONG = 1.0e15   # 强反型起始面密度 (m^-2) = 1e11 cm^-2；
# 只有源端 N_s 超过该值才认为存在可被夹断的强反型沟道，CLM 才适用


def _cumtrapz(y, x):
    out = np.zeros_like(x)
    if len(x) > 1:
        dx = np.diff(x)
        out[1:] = np.cumsum(0.5 * (y[:-1] + y[1:]) * dx)
    return out


def _clm_length(vd, vdsat, na=C.NA):
    if vd <= vdsat:
        return 0.0
    return float(np.sqrt(2.0 * C.EPS_SI * (vd - vdsat) / (C.Q_E * na)))


def gca_solve(vg, vd, table, y=None, nphi=401, nmin=1e12, clm=True):
    """
    固定 (V_G, V_D) 解 GCA，返回电流与沿沟道剖面。

    table : gate_boundary.GateTable（已 build）

    返回 dict：
        ID, vdsat, saturated, y_p, y, phi_n, ns, phi_s, id_local
    """
    if y is None:
        y = np.linspace(0.0, C.L_CH, 81)
    y = np.asarray(y, dtype=float)
    L = C.L_CH

    vdsat_cross = table.pinch_off_potential(vg, nmin=nmin)

    # 源端反型强度决定饱和/CLM 的处理方式：
    #   强反型（N_s(0) >= N_MIN_STRONG）：存在可夹断的沟道，积分上限
    #       取 vdsat（N_s 跌破 nmin 处），V_D > vdsat 后施加一阶 CLM；
    #   弱反型/深亚阈：不存在沟道夹断，N_s(phi_n) 指数拖尾对 V_D 的
    #       全积分自然收敛（V_D > 若干 V_T 后即平台），严禁截断与 CLM
    #       （开发期实测：截断会丢掉拖尾承载的大部分电流，V_G≈0.87 处
    #        I_D 出现 2.4 倍伪跌落；CLM 则制造伪饱和扭结）。
    ns0 = float(table.solve_vg(vg, 0.0)[1])
    strong = ns0 >= N_MIN_STRONG
    if ns0 < nmin:
        # 深亚阈：N_s(0) 已低于 nmin，纯扩散尾电流
        v_eff = vd
        saturated = False
        vdsat = vd
    elif strong:
        v_eff = min(vd, vdsat_cross)
        vdsat = vdsat_cross
        saturated = bool(vd > v_eff)
    else:
        v_eff = vd
        vdsat = vd            # 语义：弱反型无夹断点
        saturated = False

    # ---- N_s(phi_n) 曲线（表格插值）并积分 ----
    phi_n_fine = np.linspace(0.0, v_eff, nphi)
    ns_fine = table.ns_vs_phi_n(vg, phi_n_fine)
    ns_fine = np.maximum(ns_fine, 0.0)
    integ = float(np.trapezoid(ns_fine, phi_n_fine))
    IDsat = (C.Q_E * C.W_CH * C.MU_N / L) * integ

    if clm and saturated:
        dL = min(_clm_length(vd, vdsat), 0.95 * L)
        ID = IDsat / (1.0 - dL / L)
        y_p = L - dL
    else:
        dL = 0.0
        ID = IDsat
        y_p = L

    # ---- 沟道剖面 y(phi_n) = (q W mu/ID) ∫_0^{phi_n} N_s dphi_n ----
    if ID > 0:
        yy = (C.Q_E * C.W_CH * C.MU_N / ID) * _cumtrapz(ns_fine, phi_n_fine)
    else:
        yy = np.zeros_like(phi_n_fine)

    y_chan = np.linspace(0.0, y_p, 201)
    phi_n_chan = np.interp(y_chan, yy, phi_n_fine)
    ns_chan = np.interp(y_chan, yy, ns_fine)

    if saturated and dL > 0.0:
        y_pinch = np.linspace(y_p, L, 51)
        phi_n_pinch = vdsat + (vd - vdsat) * (y_pinch - y_p) / (L - y_p)
        ns_pinch = np.full_like(y_pinch, nmin)
        y_prof = np.concatenate([y_chan, y_pinch[1:]])
        phi_n_prof = np.concatenate([phi_n_chan, phi_n_pinch[1:]])
        ns_prof = np.concatenate([ns_chan, ns_pinch[1:]])
    else:
        y_prof = y_chan
        phi_n_prof = phi_n_chan
        ns_prof = ns_chan

    phi_s_prof = np.array([table.solve_vg(vg, pn)[0] for pn in phi_n_prof])

    # ---- 重采样到用户 y 网格 ----
    phi_n_on_y = np.interp(y, y_prof, phi_n_prof)
    ns_on_y = np.interp(y, y_prof, ns_prof)
    phi_s_on_y = np.interp(y, y_prof, phi_s_prof)

    # ---- 局部电流（连续性校验，沟道区）----
    dphi_dy = np.gradient(phi_n_on_y, y)
    id_local = C.Q_E * C.W_CH * C.MU_N * ns_on_y * dphi_dy

    return dict(ID=float(ID), vdsat=vdsat, saturated=saturated, y_p=y_p,
                y=y, phi_n=phi_n_on_y, ns=ns_on_y, phi_s=phi_s_on_y,
                id_local=id_local)


# ==============================================================================
# 偏压扫描（原 sweep.py，并入本模块）
# ==============================================================================

# -*- coding: utf-8 -*-
"""
偏压扫描（I_D-V_G、I_D-V_D）与器件工作点求解（并入 gca 模块）。

复用 table.SliceTable（一次构建）与本模块 gca_solve，避免逐点重复求解。
"""

def id_vg(table, vg_list, vd, **kw):
    """固定 V_D，扫描 V_G，返回 I_D 数组。"""
    return np.array([gca_solve(vg, vd, table, **kw)["ID"] for vg in vg_list])


def id_vd(table, vg, vd_list, **kw):
    """固定 V_G，扫描 V_D，返回 (I_D, vdsat, y_p) 数组。"""
    ids, vdsats, yps = [], [], []
    for vd in vd_list:
        r = gca_solve(vg, vd, table, **kw)
        ids.append(r["ID"])
        vdsats.append(r["vdsat"])
        yps.append(r["y_p"])
    return np.array(ids), np.array(vdsats), np.array(yps)


def transfer_characteristics(table, vg_list, vd_lin=0.05, vd_sat=1.5, **kw):
    """
    转移特性：线性区（小 V_D）与饱和区（大 V_D）的 I_D-V_G。
    返回 dict(vg, id_lin, id_sat)。
    """
    return dict(vg=np.asarray(vg_list, dtype=float),
                id_lin=id_vg(table, vg_list, vd_lin, **kw),
                id_sat=id_vg(table, vg_list, vd_sat, **kw))


def output_characteristics(table, vg_list, vd_list, **kw):
    """
    输出特性：多个 V_G 的 I_D-V_D 曲线族。
    返回 dict(vd, ids, vdsats, yps)，ids 形状 (len(vg_list), len(vd_list))。
    """
    vd = np.asarray(vd_list, dtype=float)
    ids = np.empty((len(vg_list), len(vd)))
    vdsats = np.empty_like(ids)
    yps = np.empty_like(ids)
    for i, vg in enumerate(vg_list):
        ids[i], vdsats[i], yps[i] = id_vd(table, vg, vd, **kw)
    return dict(vd=vd, vg=vg_list, ids=ids, vdsats=vdsats, yps=yps)
