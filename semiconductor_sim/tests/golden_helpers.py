# -*- coding: utf-8 -*-
"""
tests/golden_helpers.py —— golden 基准数字提取（scripts/make_golden.py 与
tests/test_golden.py 共用同一降配配置）。

设计（docs/spec/reproduce_spec.md 第 3 节、code_spec.md 第 6 节）：
    1. golden 必须秒级可复现：三阶段均采用"快速降配"（远小于生产的
       网格点数 / 表格规模 / 扫描点数），保证 pytest 全程在分钟级；
    2. 完整生产数字以 results/stageX/terminal.txt 为准，golden 只锁
       "代码重构 / 环境迁移不改变数值"这一条回归线；
    3. 提取量全部为标量（float/int），由 common.io.compare_golden 逐键
       rtol=1e-6 比对（机器确定性求解链的预期复现精度）；
    4. 三个 stage_numbers() 函数被 make_golden.py 与 test_golden.py
       同时调用，配置只有这一份，杜绝"基准一套配置、测试另一套配置"。

降配配置与生产配置对应：
    stage1 : NGRID=201（生产 801，网格收敛性另由 analyze_grid_convergence 覆盖）
    stage2 : L=2 um、正偏 1 点 + 反偏 1 点（生产 L=5 um、16+50 点）
    stage3 : 与生产共用400x220 log配置，81点栅压扫描
"""

import numpy as np

from semiconductor_sim.common import constants as C
from semiconductor_sim.common.grid import make_surface_grid
from semiconductor_sim.stage1_mos.run import run as s1_run
from semiconductor_sim.stage2_pn.run import run as s2_run
from semiconductor_sim.stage3_mosfet.table import SliceTable
from semiconductor_sim.stage3_mosfet.gca import id_vg
from semiconductor_sim.stage3_mosfet import postprocess as pp

# ---- 三阶段降配配置（唯一来源；make_golden 与 test_golden 共享）----
S1_CFG = dict(NGRID=201, verbose=False)              # 生产 NGRID=801
S2_CFG = dict(verbose=False, L=2.0e-6,               # 生产 L=5.0e-6
              v_forward=[0.5], v_reverse=[-1.0])     # 生产 16 + 50 点
from semiconductor_sim.stage3_mosfet.config import PRODUCTION
S3_TABLE_CFG = PRODUCTION.to_dict()
S3_VG = np.linspace(0.0, 2.0, 81)                    # 与 test_stage3 同网格
S3_VD_LIN, S3_VD_SAT = 0.05, 1.5                     # 与生产 run() 一致


# ============================================================================
# stage1：MOS 电容 Q-V / C-V 标量
# ============================================================================

def stage1_numbers():
    """
    S1 降配运行（NGRID=201），提取阈值/耗尽/平带三类验证标量。

    对标 run() 验证清单前三项（VALIDATION）+ Newton 迭代统计（SELF_CHECK）。
    """
    r = s1_run(**S1_CFG)
    return dict(
        vt_numeric=float(r["Vt_numeric"]),
        vt_analytic=float(r["Vt_analytic"]),
        w_num=float(r["W_num"]),
        w_max_analytic=float(r["W_max_analytic"]),
        cfb_numeric=float(r["Cfb_numeric"]),
        cfb_analytic=float(r["Cfb_analytic"]),
        psi_F=float(r["psi_F"]),
        L_D=float(r["L_D"]),
        Cox=float(r["Cox"]),
        iter_mean=float(r["iter_mean"]),
        iter_max=int(r["iter_max"]),
    )


# ============================================================================
# stage2：PN 结平衡态 + 偏压点标量
# ============================================================================

def stage2_numbers():
    """
    S2 降配运行（L=2 um，V_a=+0.5/-1.0 两个偏压点），提取：

        ver 八项   ：平衡态四重验证（V_bi / np / J / W / iters）
        Js_exact   ：Shockley 解析饱和电流
        J/W@±V_a   ：两个偏压点的电流与耗尽宽度（records 逐点核对）

    与 test_stage_closure 的 s2 fixture 同配置，闭环测试与 golden
    互为印证。
    """
    r = s2_run(**S2_CFG)
    ver = r["results"]["ver"]
    recs = {round(float(rec["V_a"]), 6): rec for rec in r["results"]["records"]}
    rec_f = recs[0.5]      # 正偏注入点
    rec_r = recs[-1.0]     # 反偏产生点
    out = dict(
        V_bi_num=float(ver["V_bi_num"]),
        V_bi_exact=float(ver["V_bi_exact"]),
        np_err=float(ver["np_err"]),
        J_eq=float(ver["J_eq"]),
        J_scale=float(ver["J_scale"]),
        W_charge=float(ver["W_charge"]),
        W_exact=float(ver["W_exact"]),
        eq_iters=int(ver["iters"]),
        Js_exact=float(r["results"]["Js_exact"]),
        J_fwd=float(rec_f["J"]),
        W_fwd=float(rec_f["W_charge"]),
        iters_fwd=int(rec_f["iters"]),
        J_rev=float(rec_r["J"]),
        W_rev=float(rec_r["W_charge"]),
        iters_rev=int(rec_r["iters"]),
    )
    return out


# ============================================================================
# stage3：MOSFET GCA 降配表格 + 转移特性标量
# ============================================================================

def _s3_table():
    """Compatibility helper; production config is centrally defined."""
    from semiconductor_sim.stage3_mosfet.phase1 import production_table
    return production_table()[0]


def stage3_numbers():
    from semiconductor_sim.stage3_mosfet.golden import production_numbers
    return production_numbers()

STAGE_NUMBERS=dict(stage1=stage1_numbers,stage2=stage2_numbers,stage3=stage3_numbers)
