# -*- coding: utf-8 -*-
"""
准费米势（QFP）列变换（S2 全 Newton 系统使用，S3/S4 预留）。

把雅可比块三对角矩阵的未知量从 (phi, n, p) 变换到 (phi, phi_n, phi_p)，
使所有未知量都是电势尺度（V），消除指数耦合造成的条件数灾难。
"""

import numpy as np


def qfp_columns(A, cn, cp):
    """
    雅可比的列变换：(phi, n, p) 变量 -> (phi, phi_n, phi_p) 准费米势变量。

    由载流子统计 n = n_i exp((phi-phi_n)/V_T)、p = n_i exp((phi_p-phi)/V_T)
    （对同一节点 k）推出变量代换的链式法则：
        d/dphi_k   = d/dphi_k + (n_k/V_T) d/dn_k - (p_k/V_T) d/dp_k
        d/dphi_n,k = -(n_k/V_T) d/dn_k
        d/dphi_p,k = +(p_k/V_T) d/dp_k

    为什么必须做这个变换：(phi, n, p) 直接作未知量时，连续性方程对
    phi 的耦合 ~ n/V_T（耗尽区 1e10 ~ 准中性区 1e22，跨 12 个数量级），
    行缩放无法归一各列，块矩阵条件数实测 ~1e36（double 只有 1e16），
    Newton 步完全失真。换成准费米势后所有未知量都是电势尺度（V），
    n、p 由指数映射自动保持正，指数耦合被解析吸收。

    参数
    ----
    A  : (M,3,3) 块三对角雅可比的某一层（对角/次/超），列为 [phi, n, p]
    cn : (M,) 该层每个块所对应节点的 n_k/V_T
    cp : (M,) 该层每个块所对应节点的 p_k/V_T

    返回变换后的 (M,3,3)，列为 [phi, phi_n, phi_p]（单位 V）。
    """
    New = np.empty_like(A)
    New[:, :, 0] = A[:, :, 0] + cn[:, None] * A[:, :, 1] \
        - cp[:, None] * A[:, :, 2]
    New[:, :, 1] = -cn[:, None] * A[:, :, 1]
    New[:, :, 2] = cp[:, None] * A[:, :, 2]
    return New
