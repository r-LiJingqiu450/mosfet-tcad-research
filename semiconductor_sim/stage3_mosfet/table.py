# -*- coding: utf-8 -*-
"""
栅极混合边界与 (phi_s, phi_n) -> (V_G, N_s, Q_s) 二维表格。

系列化说明：单截面求解与电荷后处理取自 stage3_mosfet.solver
（其内部复用 common.thomas / common.constants）。

栅极边界条件（理想氧化层，无界面态/固定电荷）：

    eps_ox (V_G - phi_s - V_FB) / t_ox = eps_si E_s = -Q_s

第一阶段已证明等价写法（把 phi_s 当自变量、反算 V_G）：

    V_G = V_FB + phi_s - Q_s/C_ox           （phi_s 相对衬底体内）

本模块沿用这一"直接扫 phi_s、反算 V_G"的方案（栅极边界残差按定义为零），
构建 (phi_s, phi_n) -> (V_G, N_s, Q_s) 二维表格，之后所有 V_G/V_D 扫描
与 GCA 电流积分都走廉价双线性插值 + 单调反演，不再重复求解。

深反型角点（phi_s - phi_n > SKIP_OV）不在器件工作范围内，且反型层
厚度低于网格分辨率、数值病态，构建时跳过（置 NaN）。该角点对应的
V_G 下限约 2.6 V（phi_n=0 列）至 3.2 V（大 phi_n 列，phi_s 已顶到
网格上限 1.6 V），高于本文 V_G <= 2 V 的全部工作点：实测 110 个
phi_n 列的有效 V_G 上限最小为 2.019 V，故 V_G=2 V 处所有列均为
内插而非钳位外推。跳过原因（直接求解实测）：(phi_s=1.45, phi_n=0)
截面 Newton 虽"收敛"，但 n_surf=2.3e28 cm^-3（超过硅原子密度）、
E_surf=3.0e7 V/cm，反型层厚度 ~V_T/E_s ~ 0.01 nm 远低于 0.125 nm
表面网格——非简并统计与网格分辨率双重失效，解无物理意义。
"""

import warnings

import numpy as np

from common import constants as C
from .solver import solve_slice, surface_charge, electron_sheet_density

SKIP_OV = 0.90   # phi_s - phi_n 超过该值视为深反型角点（对应 V_G >~ 2.6 V），跳过



def gate_voltage(phi_s, qs, vfb=C.VFB):
    """由表面势 phi_s（相对衬底体内）与半导体面电荷 Q_s 反算栅压 V_G。"""
    return vfb + phi_s - qs / C.COX


def gate_residual(phi_s, qs, vg, vfb=C.VFB):
    """栅极边界残差（诊断用）：eps_si E_s - C_ox (V_G - V_FB - phi_s)。"""
    es = -qs / C.EPS_SI
    return C.EPS_SI * es - C.COX * (vg - vfb - phi_s)


class SliceTable:
    """
    (phi_s, phi_n) -> (V_G, N_s, Q_s) 二维规则网格表格（直接扫 phi_s）。

    属性
    ----
    phi_s_grid : (Nps,) 表面势采样（相对衬底体内）
    phi_n_grid : (Npn,) 沟道电子准费米势采样
    vg, ns, qs : (Nps, Npn) 对应输出（深反型角点为 NaN）
    x          : 纵向网格
    """

    def __init__(self, phi_s_grid, phi_n_grid, x, na=C.NA, sheet_interpolation="log"):
        if (len(phi_s_grid) < 2 or len(phi_n_grid) < 1
                or np.any(np.diff(phi_s_grid) <= 0)
                or np.any(np.diff(phi_n_grid) <= 0)
                or not np.all(np.isfinite(phi_s_grid))
                or not np.all(np.isfinite(phi_n_grid))):
            raise ValueError('Table axes must be finite and strictly increasing')
        self.phi_s_grid = np.asarray(phi_s_grid, dtype=float)
        self.phi_n_grid = np.asarray(phi_n_grid, dtype=float)
        self.x = x
        if sheet_interpolation not in ('linear','log'):
            raise ValueError('sheet_interpolation must be linear or log')
        self.sheet_interpolation = sheet_interpolation
        self.na = na
        self.nps = len(self.phi_s_grid)
        self.npn = len(self.phi_n_grid)
        self.vg = np.full((self.nps, self.npn), np.nan)
        self.ns = np.full((self.nps, self.npn), np.nan)
        self.qs = np.full((self.nps, self.npn), np.nan)

    # ------------------------------------------------------------------
    def build(self, verbose=False, tol=1e-9, hot=True):
        """
        构建表格：phi_s 外层（从平带点向两端）、phi_n 内层。

        phi_s 外层热启动（hot=True）：相邻 phi_s 行间隔仅 ~0.008 V，
        同 phi_n 列直接沿用上一行解作初值。收敛判据放宽到 1e-9 V 后
        （见 poisson_dd_1d 文档字符串），热启动 2~4 步收敛，且消除了
        "真解附近因舍入噪声无法收敛"的停顿；任一列失败即回退冷启动
        （耗尽近似初值，全部工作点 2~7 步收敛、永不失败）。
        phi_n 内层不延续（陈旧 phi_n 剖面与新区界冲突，经验会失败）。

        hot=False 时全部冷启动（用于校验热启动等价性）。
        """
        self.build_stats=dict(potential_min=float('inf'),potential_max=float('-inf'),
                              clip_evaluations=0,solved_slices=0,skipped_slices=0,attempted_solves=0,failed_hot_attempts=0)
        def record_attempt(result):
            self.build_stats['attempted_solves']+=1
            for key in ['potential_min','potential_max']:
                op=min if key=='potential_min' else max
                self.build_stats[key]=op(self.build_stats[key],result.get(key,0.))
            self.build_stats['clip_evaluations']+=result.get('clip_evaluations',0)
        order = np.argsort(np.abs(self.phi_s_grid))       # 平带点起步
        n_fail = 0
        n_skip = 0
        n_hot = 0
        prev_phi = [None] * self.npn                       # 每列上一行解

        for rank, ips in enumerate(order):
            phi_s = self.phi_s_grid[ips]
            for ipn in range(self.npn):
                phi_n = self.phi_n_grid[ipn]
                # 深反型角点：跳过（不在工作范围、反型层低于网格分辨率）
                if phi_s - phi_n > SKIP_OV:
                    n_skip += 1
                    self.build_stats["skipped_slices"]+=1
                    continue
                res = None
                if hot and prev_phi[ipn] is not None:
                    res = solve_slice(phi_s, phi_n, self.x, Na=self.na,
                                      state=(prev_phi[ipn], None, None),
                                      tol=tol, warn=False)
                    record_attempt(res)
                    if res["converged"]:
                        n_hot += 1
                    else:
                        self.build_stats["failed_hot_attempts"]+=1
                        res = None
                if res is None:
                    res = solve_slice(phi_s, phi_n, self.x, Na=self.na,
                                      tol=tol, warn=False)
                    record_attempt(res)
                if not res['converged']:
                    raise RuntimeError(f'Invalid slice at {phi_s=}, {phi_n=}: '
                                       f"{res.get('status','not_converged')}")
                self.build_stats['solved_slices']+=1
                prev_phi[ipn] = res["phi"]
                qs_val, _ = surface_charge(res, self.x)
                self.vg[ips, ipn] = gate_voltage(phi_s, qs_val, vfb=C.VFB)
                self.qs[ips, ipn] = qs_val
                self.ns[ips, ipn] = electron_sheet_density(res, self.x)
                if not res["converged"]:
                    n_fail += 1
            if verbose:
                print(f"  phi_s = {phi_s:+.4f} V 行完成")
        if n_fail:
            warnings.warn(f"SliceTable.build: {n_fail} 个截面未收敛。",
                          RuntimeWarning, stacklevel=2)
        self.n_hot = n_hot
        return self

    # ------------------------------------------------------------------
    @staticmethod
    def _bilinear(g1, g2, data, v1, v2):
        v1 = float(np.clip(v1, g1[0], g1[-1]))
        v2 = float(np.clip(v2, g2[0], g2[-1]))
        i1 = int(np.clip(np.searchsorted(g1, v1) - 1, 0, len(g1) - 2))
        i2 = int(np.clip(np.searchsorted(g2, v2) - 1, 0, len(g2) - 2))
        t1 = (v1 - g1[i1]) / (g1[i1 + 1] - g1[i1])
        t2 = (v2 - g2[i2]) / (g2[i2 + 1] - g2[i2])
        c00 = data[i1, i2]
        c10 = data[i1 + 1, i2]
        c01 = data[i1, i2 + 1]
        c11 = data[i1 + 1, i2 + 1]
        return (c00 * (1 - t1) * (1 - t2) + c10 * t1 * (1 - t2)
                + c01 * (1 - t1) * t2 + c11 * t1 * t2)

    # ------------------------------------------------------------------
    def solve_vg(self, vg, phi_n):
        """
        给定 V_G 与 phi_n，反解 phi_s 使 gate_voltage = vg，返回
        (phi_s, N_s, Q_s)。V_G 随 phi_s 单调递增，故沿 phi_s 轴线性反演。
        """
        # phi_n 方向：定位括号区间 [ipn, ipn+1]，t2 ∈ [0,1] 为内插系数
        # （早期版本用 searchsorted 的上节点作基准、t2 为负，等于向下
        #  外推半个单元，会放大 phi_n 网格敏感性并引入网格依赖误差）
        if (not np.isfinite(vg+phi_n) or phi_n < self.phi_n_grid[0]
                or phi_n > self.phi_n_grid[-1] or self.npn < 2):
            raise ValueError('Query outside supported table phi_n domain')
        ipn = int(np.clip(np.searchsorted(self.phi_n_grid, phi_n) - 1, 0,
                          self.npn - 2))
        ipn_lo = ipn
        ipn_hi = ipn + 1
        t2 = 0.0
        if self.phi_n_grid[ipn_hi] > self.phi_n_grid[ipn_lo]:
            t2 = (phi_n - self.phi_n_grid[ipn_lo]) / \
                 (self.phi_n_grid[ipn_hi] - self.phi_n_grid[ipn_lo])

        def invert_1d(vg_col_1d, data_col_1d, logarithmic=False):
            valid = np.isfinite(vg_col_1d)
            if valid.sum() < 2:
                return np.nan, np.nan
            g = self.phi_s_grid[valid]
            vcol = vg_col_1d[valid]
            dcol = data_col_1d[valid]
            if np.any(np.diff(vcol) <= 0):
                raise ValueError('Nonmonotone gate-voltage table column')
            if vg < vcol[0] or vg > vcol[-1]:
                raise ValueError('VG query outside valid table domain')
            if vg == vcol[0]:
                return g[0], dcol[0]
            if vg == vcol[-1]:
                return g[-1], dcol[-1]
            i = int(np.searchsorted(vcol, vg)) - 1
            i = min(max(i, 0), len(vcol) - 2)
            tt = (vg - vcol[i]) / (vcol[i + 1] - vcol[i])
            value = (np.exp((1-tt)*np.log(dcol[i])+tt*np.log(dcol[i+1]))
                     if logarithmic else dcol[i]+tt*(dcol[i+1]-dcol[i]))
            return g[i]+tt*(g[i+1]-g[i]), value

        phi_s_lo, ns_lo = invert_1d(self.vg[:, ipn_lo], self.ns[:, ipn_lo], self.sheet_interpolation=="log")
        phi_s_hi, ns_hi = invert_1d(self.vg[:, ipn_hi], self.ns[:, ipn_hi], self.sheet_interpolation=="log")
        qs_lo = invert_1d(self.vg[:, ipn_lo], self.qs[:, ipn_lo])[1]
        qs_hi = invert_1d(self.vg[:, ipn_hi], self.qs[:, ipn_hi])[1]
        phi_s = phi_s_lo + t2 * (phi_s_hi - phi_s_lo)
        ns = (np.exp((1-t2)*np.log(ns_lo)+t2*np.log(ns_hi))
              if self.sheet_interpolation=='log' else ns_lo+t2*(ns_hi-ns_lo))
        qs = qs_lo + t2 * (qs_hi - qs_lo)
        return float(phi_s), float(ns), float(qs)

    def ns_vs_phi_n(self, vg, phi_n_arr):
        """固定 V_G，输出 N_s(phi_n) 数组（GCA 积分用）。"""
        out = np.empty(len(phi_n_arr))
        for k, pn in enumerate(phi_n_arr):
            _, ns, _ = self.solve_vg(vg, pn)
            out[k] = ns
        return out

    def pinch_off_potential(self, vg, nmin=1e12):
        """夹断电势 phi_n,p：N_s(phi_n) 首次低于 nmin 的 phi_n（二分反演）。"""
        ns = self.ns_vs_phi_n(vg, self.phi_n_grid)
        if ns[-1] > nmin:
            return float(self.phi_n_grid[-1])
        idx = int(np.where(ns < nmin)[0][0])
        if idx == 0:
            return float(self.phi_n_grid[0])
        pn0, pn1 = self.phi_n_grid[idx - 1], self.phi_n_grid[idx]
        ns0, ns1 = ns[idx - 1], ns[idx]
        return float(pn0 + (nmin - ns0) * (pn1 - pn0) / (ns1 - ns0))
