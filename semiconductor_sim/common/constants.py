# -*- coding: utf-8 -*-
"""
系列统一物理常数与器件参数（S1/S2/S3 共用）。

数值约定（全系列一致）：
    - 内部全部 SI 单位（m、m^-3、V、A/m^2、s、F/m^2），输出层转换 cm；
    - 电势 phi 以衬底/体区参考为基准（S1/S3: phi_bulk = 0；
      S2: phi_p = 0，phi_n = V_a）；
    - 载流子 Boltzmann 统计：n = n_i exp((phi-phi_n)/V_T)，
      p = n_i exp((phi_p-phi)/V_T)；
    - 300 K 热电压 V_T = kT/q = 25.852 mV。

本文件是全系列唯一参数来源：任何阶段不得自定义重复常数。
器件参数按"阶段块"组织；跨阶段一致的参数放在通用块，
标注其被哪些阶段使用。
"""

import numpy as np

# ================================================================================
# 1. 通用物理常数（三阶段完全一致）
# ================================================================================
Q_E = 1.602176634e-19      # 元电荷 q (C)
K_B = 1.380649e-23         # 玻尔兹曼常数 k (J/K)
EPS0 = 8.8541878128e-12    # 真空介电常数 (F/m)
TEMP = 300.0               # 温度 (K)

VT = K_B * TEMP / Q_E      # 热电压 kT/q (V)，300 K 下 0.025852 V

EPS_SI = 11.7 * EPS0       # 硅介电常数 (F/m)
EPS_OX = 3.9 * EPS0        # 二氧化硅介电常数 (F/m)

NI = 1.0e16                # 硅本征载流子浓度 (m^-3)，即 1e10 cm^-3
NI_CM3 = 1.0e10            # 硅本征载流子浓度 (cm^-3)

# 硅能带参数（S1 能带图绘制用），300 K
EC_MINUS_EI = 0.56          # Ec - Ei (eV)
EI_MINUS_EV = 0.56          # Ei - Ev (eV)
EG = 1.12                   # 硅禁带宽度 (eV)（S2 能带图用）

# ================================================================================
# 2. 迁移率 / 扩散系数 / 复合寿命（S2、S3 使用）
# ================================================================================
MU_N = 0.135               # 电子迁移率 (m^2/Vs) = 1350 cm^2/Vs
MU_P = 0.045               # 空穴迁移率 (m^2/Vs) = 450  cm^2/Vs
D_N = MU_N * VT            # 电子扩散系数 (m^2/s)，Einstein 关系
D_P = MU_P * VT            # 空穴扩散系数 (m^2/s)

TAU_N = 1.0e-10            # 电子寿命 (s)
TAU_P = 1.0e-10            # 空穴寿命 (s)

# ================================================================================
# 3. 器件参数（跨阶段一致部分）
# ================================================================================
NA_CM3 = 1.0e16            # p 型区受主浓度 (cm^-3)：S1 衬底 / S2 p 侧 / S3 衬底
NA = NA_CM3 * 1e6          # (m^-3)

ND_CM3 = 1.0e16            # S2 n 侧施主浓度 (cm^-3)
ND = ND_CM3 * 1e6          # (m^-3)

TOX_NM = 20.0              # 栅氧化层厚度 (nm)：S1、S3 一致
TOX = TOX_NM * 1e-9        # (m)
COX = EPS_OX / TOX         # 栅氧单位面积电容 (F/m^2) = 1.727e-3

VFB = 0.0                  # 平带电压 (V)：S1、S3 一致（第一版未含 Phi_MS）

# ---- S3 沟道几何 ----
L_CH = 1.0e-6              # 沟道长度 (m)
W_CH = 1.0e-6              # 沟道宽度 (m)（电流按 A/μm 归一）

# ================================================================================
# 4. 派生参考量（由上述参数计算，供解析对照）
# ================================================================================
PSI_F = VT * np.log(NA / NI)        # 费米势 V_T ln(Na/n_i) ≈ 0.357 V
PHI_BULK = -PSI_F                    # 衬底体内电势（phi_p=0 参考）
TWICE_PSI_F = 2.0 * PSI_F            # 反型判据 2·psi_F ≈ 0.714 V

GAMMA = np.sqrt(2.0 * EPS_SI * Q_E * NA) / COX          # 衬偏效应因子 (V^0.5)
VT_LONG = VFB + TWICE_PSI_F + GAMMA * np.sqrt(TWICE_PSI_F)   # S3 长沟道阈值 (V)

L_D = np.sqrt(EPS_SI * VT / (Q_E * NA))   # 德拜长度 (m) ≈ 41 nm（网格设计依据）


# ================================================================================
# 5. 解析参考公式（各阶段"验证"的独立参照）
# ================================================================================
def analytic_psi_f(Na=NA):
    """费米势 psi_F = V_T ln(Na/n_i)。"""
    return VT * np.log(Na / NI)


def analytic_vt(vfb=VFB, na=NA, cox=COX):
    """S1/S3 长沟道阈值电压（耗尽近似 + 体电荷模型）。"""
    psif = analytic_psi_f(na)
    gamma = np.sqrt(2.0 * EPS_SI * Q_E * na) / cox
    return vfb + 2.0 * psif + gamma * np.sqrt(2.0 * psif)


def square_law_id(vg, vd, vfb=VFB, na=NA, cox=COX, mu_n=MU_N,
                  L=L_CH, W=W_CH):
    """
    S3 长沟道 Square-law 漏电流（理想，忽略沟道长度调制）：
        线性区 (vd <= vg - vt):  I_D = mu_n Cox (W/L) [(vg-vt)vd - vd^2/2]
        饱和区 (vd >  vg - vt):  I_D = (1/2) mu_n Cox (W/L) (vg-vt)^2
    返回电流 (A)。
    """
    vt = analytic_vt(vfb, na, cox)
    vg = np.asarray(vg, dtype=float)
    vd = np.asarray(vd, dtype=float)
    ov = vg - vt                       # 过驱动电压
    vdsat = np.maximum(ov, 0.0)
    id_lin = mu_n * cox * (W / L) * (ov * vd - 0.5 * vd ** 2)
    id_sat = 0.5 * mu_n * cox * (W / L) * vdsat ** 2
    return np.where(ov <= 0.0, 0.0, np.where(vd < vdsat, id_lin, id_sat))


def V_bi_analytic(Na=NA, Nd=ND):
    """S2 pn 结内建电势 V_bi = V_T ln(Na·Nd/n_i^2)。"""
    return VT * np.log(Na * Nd / NI ** 2)


def W_analytic(Na, Nd, V_a=0.0):
    """
    S2 耗尽层总宽度 W = sqrt(2 eps/q (1/Na + 1/Nd) (V_bi - V_a))。
    V_a > 0 为正偏（势垒降低）。
    """
    vbi = V_bi_analytic(Na, Nd)
    return np.sqrt(2.0 * EPS_SI / Q_E * (1.0 / Na + 1.0 / Nd) * (vbi - V_a))


def Js_analytic(Na, Nd):
    """S2 反偏饱和电流密度 J_s = q n_i^2 (D_p/(L_p N_d) + D_n/(L_n N_a))
    （长基区 Shockley 近似，少子寿命 tau 取 TAU_N/TAU_P）。"""
    L_p = np.sqrt(D_P * TAU_P)
    L_n = np.sqrt(D_N * TAU_N)
    return Q_E * NI ** 2 * (D_P / (L_p * Nd) + D_N / (L_n * Na))


def equilibrium_carriers(Na=NA, Nd=ND):
    """S2 平衡载流子（p 侧 / n 侧），返回 (p_p, n_p, p_n, n_n)（m^-3）。"""
    p_p = 0.5 * (Na + np.sqrt(Na ** 2 + 4.0 * NI ** 2))
    n_p = NI ** 2 / p_p
    n_n = 0.5 * (Nd + np.sqrt(Nd ** 2 + 4.0 * NI ** 2))
    p_n = NI ** 2 / n_n
    return p_p, n_p, p_n, n_n
