# -*- coding: utf-8 -*-
"""
系列统一网格生成。

收录三阶段实际使用的全部网格构造：
    make_grid            S1 均匀 / 非均匀（二分求解几何增长比）
    make_surface_grid    S3 表面密铺三段式（栅界面 x=0 起）
    make_junction_grid   S2 对称三段式（结区 x=0 居中）
    make_channel_grid    S3 沟道截面均匀网格
    grid_stats           网格统计（点数、首层/末层间距、增长比）

设计约束（全系列一致）：
    - 表面/结区首层间距 d0 由最陡特征尺度决定（反型层衰减长度
      V_T/E_s 或德拜长度 L_D ≈ 41 nm）；
    - 最大间距 h_max 必须小于特征屏蔽长度（经验 h_max <= L_D/8），
      否则过渡区空间电荷变化无法解析；
    - 相邻间距增长比 <= 1.2（S1/S2）或 1.15（S3 表面段），保证格式精度。
"""

import numpy as np


# ================================================================================
# S1：均匀 / 非均匀网格
# ================================================================================
def _solve_growth_ratio(h0, total, n_segments):
    """
    数值解求几何增长比 r：
        h0 * (r^n - 1) / (r - 1) = total
    n_segments 段，首层间距 h0，总跨度 total。
    左侧关于 r 单调递增，二分求解；r=1 时退化为 h0*n。
    """
    if n_segments < 1:
        return 1.0

    def sum_of_ratio(r):
        if abs(r - 1.0) < 1e-12:
            return h0 * n_segments
        return h0 * (r ** n_segments - 1.0) / (r - 1.0)

    lo, hi = 1.0, 1.0
    for _ in range(200):
        if sum_of_ratio(hi) >= total:
            break
        hi *= 2.0
    if sum_of_ratio(lo) >= total:
        return 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if sum_of_ratio(mid) >= total:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def make_grid(L, N, kind="uniform", d0=1.0e-9, transition=0.35e-6,
              h_max=4.5e-9, max_growth=1.2):
    """
    生成一维网格（S1 MOS 电容）。

    参数
    ----
    L          : 域长度 (m)
    N          : 总点数（均匀网格严格等于 N；非均匀网格近似）
    kind       : 'uniform' | 'nonuniform'
    d0         : 非均匀网格表面首层间距 (m)
    transition : 细密段 [0, transition] 的末端位置 (m)
    h_max      : 非均匀网格最大间距 (m)，须 < 德拜长度/8
    max_growth : 相邻间距增长比上限（通常 <= 1.2，保证格式精度）

    非均匀网格结构（三段拼接）：
        段 A  [0, transition]          均匀细网格 (d0)
        段 B  (transition, x_B]        几何增长过渡 (d0 -> h_max, 逐层 ≤1.2)
        段 C  (x_B, L]                 均匀粗网格 (h_max)

    返回严格单调递增数组。
    """
    import warnings

    if kind == "uniform":
        return np.linspace(0.0, L, N)

    # ---- 段 A：表面细网格 ----
    N1 = max(2, int(round(transition / d0)) + 1)
    if d0 * (N1 - 1) >= L:
        warnings.warn("make_grid: transition/d0 组合超出域长，退化为均匀网格",
                      RuntimeWarning, stacklevel=2)
        return np.linspace(0.0, L, N1)
    segA = np.linspace(0.0, transition, N1)
    hA = segA[1] - segA[0]

    # ---- 段 B：几何过渡（间距从 hA 增长到 h_max）----
    ratio = h_max / hA
    if ratio <= 1.02:
        segB = np.array([transition + hA])
    else:
        kB = max(2, int(np.ceil(np.log(ratio) / np.log(max_growth))))
        r = ratio ** (1.0 / kB)
        segB = transition + hA * (r ** np.arange(1, kB + 1) - 1.0) / (r - 1.0)
        segB = segB[segB < L - h_max]         # 保证段 C 至少有两层
        if len(segB) == 0:
            segB = np.array([transition + hA])

    # ---- 段 C：粗均匀网格至远场 ----
    N3 = max(3, int(round((L - segB[-1]) / h_max)) + 1)
    segC = np.linspace(segB[-1], L, N3)

    x = np.concatenate([segA, segB, segC])
    return np.unique(x)


# ================================================================================
# S3：栅界面表面密铺三段式
# ================================================================================
def make_surface_grid(x_max, d0=0.5e-9, h_max=4.0e-9, transition=0.08e-6,
                      max_growth=1.15):
    """
    S3 纵向 x 网格：x[0]=0（栅界面），x[-1]=x_max（衬底欧姆接触）。

    结构：表面 [0, transition] 内均匀步长 d0（反型层厚度仅 ~5-10 nm，
    且表面电荷用三点单侧差分要求前两段步长相等），随后按 max_growth
    几何增长至 h_max，再保持 h_max 直至 x_max。

    返回严格单调递增数组。
    """
    if not (d0 > 0 and h_max >= d0 and x_max > transition > 0):
        raise ValueError("网格参数非法（需 0 < d0 <= h_max, 0 < transition < x_max）")

    xs = [0.0]
    # ---- 表面均匀段 ----
    while xs[-1] + d0 < transition - 1e-18:
        xs.append(xs[-1] + d0)
    xs.append(transition)
    # ---- 几何增长段 ----
    h = d0
    while xs[-1] < x_max - 1e-18:
        h = min(h * max_growth, h_max)
        if xs[-1] + h >= x_max:
            h = x_max - xs[-1]
        xs.append(xs[-1] + h)
    return np.asarray(xs, dtype=float)


# ================================================================================
# S2：pn 结对称三段式
# ================================================================================
def make_junction_grid(L, d0=1.5e-9, transition=0.6e-6, h_max=8e-9,
                       max_growth=1.2):
    """
    S2 以结区 (x=0) 为中心的对称三段式非均匀网格，域 [-L/2, L/2]。

    设计依据：耗尽层宽度 W_0 ~ 0.43 um（平衡）~1.2 um（反偏 5 V），
    德拜长度 L_D = sqrt(eps V_T/(q N)) ~ 40 nm。网格须解析结区
    空间电荷的陡峭分布：
        段 A  [-transition, transition]   均匀细网格 (d0 = 1.5 nm)
        段 B  过渡段                      几何增长 (d0 -> h_max, 逐层 <=1.2)
        段 C  远场                       均匀粗网格 (h_max = 8 nm < L_D/5)

    返回 x 数组 (m)，保证 x=0 恰为节点（结位于节点上）。
    """
    half = 0.5 * L

    # ---- 右半边 [0, half] ----
    N1 = max(2, int(round(transition / d0)) + 1)
    segA = np.linspace(0.0, transition, N1)
    hA = segA[1] - segA[0]

    # 段 B：几何过渡（间距从 hA 增长到 h_max）
    ratio = h_max / hA
    if ratio <= 1.02:
        segB = np.array([transition + hA])
    else:
        kB = max(2, int(np.ceil(np.log(ratio) / np.log(max_growth))))
        r = ratio ** (1.0 / kB)
        segB = transition + hA * (r ** np.arange(1, kB + 1) - 1.0) / (r - 1.0)
        segB = segB[segB < half - h_max]
        if len(segB) == 0:
            segB = np.array([transition + hA])

    # 段 C：远场粗网格
    N3 = max(3, int(round((half - segB[-1]) / h_max)) + 1)
    segC = np.linspace(segB[-1], half, N3)

    xr = np.unique(np.concatenate([segA, segB, segC]))
    # 镜像拼接：结处节点共享（不重复）
    x = np.concatenate([-xr[::-1], xr[1:]])
    return x


# ================================================================================
# S3：沟道截面网格
# ================================================================================
def make_channel_grid(L, Ny=41):
    """S3 沟道截面位置：y ∈ [0, L] 均匀 Ny 个点（含两端源/漏）。"""
    if Ny < 3:
        raise ValueError("Ny 至少为 3")
    return np.linspace(0.0, L, Ny)


# ================================================================================
# 统计
# ================================================================================
def grid_stats(x):
    """网格统计：点数、表面/体内步长、覆盖率。"""
    h = np.diff(x)
    return dict(N=len(x), h0=h[0], h_min=h.min(), h_max=h.max(),
                L=x[-1], ratio=h.max() / h[0])
