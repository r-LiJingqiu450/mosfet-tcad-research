# -*- coding: utf-8 -*-
"""
参数提取与验证（第三阶段）：
    阈值电压 V_T（恒定电流判据）、亚阈值摆幅 SS、跨导 g_m、
    输出电阻 r_o、沟道电流连续性校验。

    系列化说明：常数取自 common.constants（全系列唯一参数源）。
"""

import numpy as np

from common import constants as C


def extract_vt(vg, id_lin, i_th=None):
    """
    恒定电流判据提取 V_T。
    默认 I_th = 0.1 μA（W/L=1 的行业常用判据，V_D 取线性区小偏压）。
    返回 V_T 及所用判据电流。
    """
    if i_th is None:
        i_th = 1e-7                 # 0.1 μA（W=L=1 μm）
    vg = np.asarray(vg, dtype=float)
    id_lin = np.asarray(id_lin, dtype=float)
    above = id_lin >= i_th
    if not above.any():
        return float(np.nan), i_th
    idx = int(np.argmax(above))
    if idx == 0:
        return float(vg[0]), i_th
    # 对数插值
    v0, v1 = vg[idx - 1], vg[idx]
    i0 = max(id_lin[idx - 1], 1e-300)
    i1 = max(id_lin[idx], 1e-300)
    vt = v0 + (np.log(i_th) - np.log(i0)) / (np.log(i1) - np.log(i0)) * (v1 - v0)
    return float(vt), i_th


def extract_ss(vg, id_log, vd=None):
    """
    亚阈值摆幅 SS（mV/dec）：在固定电流窗口对 log10(I_D)-V_G 线性拟合。

    窗口逐级下试（1e-14~1e-10 / 1e-13~1e-9 / 1e-12~1e-8 A），取第一个
    样本数 >= 5 的窗口。窗口位于深亚阈区（约 4 个数量级低于阈值电流），
    该处耗尽层较薄、C_dep/C_ox 较大，SS 高于阈值点值——这正是衬偏效应
    的表现，与理论 SS = ln10·V_T·(1+C_dep/C_ox) 逐点对照见 ss_theory。
    注意：本函数不是"最陡拟合"，早期论文稿中"滑动窗口最陡值"的表述
    与实现不符，已更正。
    """
    vg = np.asarray(vg, dtype=float)
    id_log = np.asarray(id_log, dtype=float)
    lg = np.log10(np.maximum(id_log, 1e-300))
    for lo, hi in [(1e-14, 1e-10), (1e-13, 1e-9), (1e-12, 1e-8)]:
        mask = (id_log >= lo) & (id_log <= hi)
        if mask.sum() >= 5:
            slope = np.polyfit(vg[mask], lg[mask], 1)[0]
            if slope > 1e-6:
                return float(1000.0 / slope)
    return float("nan")


def local_ss(vg, id_):
    """
    局部亚阈值摆幅 SS(V_G) = dV_G/dlog10(I_D)（mV/dec）曲线。
    用于检验 SS 随 V_G 的演变（深亚阈大、近阈值小、强反型再上升），
    其近阈值最小值应与衬偏效应理论 ss_theory(phi_s≈2psi_F) 一致。
    """
    vg = np.asarray(vg, dtype=float)
    id_ = np.asarray(id_, dtype=float)
    lg = np.log10(np.maximum(id_, 1e-300))
    slope = np.gradient(lg, vg)
    out = np.where(slope > 1e-9, 1000.0 / slope, np.nan)
    return out


def ss_theory(phi_s, cox=C.COX, eps=C.EPS_SI, na=C.NA):
    """
    衬偏效应理论亚阈值摆幅 SS = ln10·V_T·(1 + C_dep/C_ox)（mV/dec），
    C_dep = eps_si/W，W = sqrt(2·eps_si·phi_s/(q·Na))（耗尽近似）。
    phi_s 较小时 W 薄、C_dep 大 → SS 大；phi_s→2psi_F 时 W→W_max，
    SS → ln10·V_T·(1+eps_si/(W_max·C_ox))（近阈值最小值）。
    """
    phi_s = np.asarray(phi_s, dtype=float)
    w = np.sqrt(2.0 * eps * np.maximum(phi_s, 1e-9) / (C.Q_E * na))
    cdep = eps / w
    return np.log(10.0) * C.VT * (1.0 + cdep / cox) * 1000.0


def extract_gm(vg, id_):
    """跨导 g_m = dI_D/dV_G（S），返回 (gm, gm_max)。"""
    vg = np.asarray(vg, dtype=float)
    id_ = np.asarray(id_, dtype=float)
    gm = np.gradient(id_, vg)
    return gm, float(gm.max())


def extract_ro(vd, id_):
    """输出电阻 r_o = (dI_D/dV_D)^{-1}（Ω）。"""
    vd = np.asarray(vd, dtype=float)
    id_ = np.asarray(id_, dtype=float)
    did = np.gradient(id_, vd)
    did = np.where(np.abs(did) > 1e-18, did, 1e-18)
    return 1.0 / did


def continuity_check(res, tol=1e-3):
    """
    沟道电流连续性校验：I_D(y) = q W mu N_s dphi_n/dy 沿沟道的相对波动。
    只统计强沟道区（N_s > 20% N_s,max）并剔除两端各 5% 的边界伪影；
    此掩码用于限制独立差分求导截断误差；不意味着准费米电流遗漏扩散。
    返回 (rel_fluctuation, mean_id, mask)。
    """
    id_local = np.asarray(res["id_local"], dtype=float)
    ns = np.asarray(res["ns"], dtype=float)
    y = np.asarray(res["y"], dtype=float)
    mask = ns > 0.2 * ns.max()          # 强沟道区
    n = len(y)
    edge = max(1, int(0.05 * n))        # 剔除两端
    if edge < n:
        mask[:edge] = False
        mask[-edge:] = False
    if mask.sum() < 3:
        return np.nan, np.nan, mask
    id_chan = id_local[mask]
    mean_id = float(np.mean(id_chan))
    rel = float(np.max(np.abs(id_chan - mean_id)) / max(abs(mean_id), 1e-30))
    return rel, mean_id, mask


def edge_continuity_check(res):
    """Internal integral consistency, not independent physical validation."""
    if 'id_edges' not in res:
        raise ValueError('No conservative edges in empirical legacy mode')
    ids=np.asarray(res['id_edges'])
    return float(np.max(np.abs(ids-res['ID']))/max(abs(res['ID']),1e-300))
