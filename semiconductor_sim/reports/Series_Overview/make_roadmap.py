#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_roadmap.py —— 生成 Series_Overview 的 S0-F1 系列路线图。

四节点横向递进（S1 MOS 电容 → S2 PN 结 → S3 MOSFET(GCA) → S4 全 DD 预留），
箭头标注阶段闭环退化关系，底部横带为 common/ 公共层。输出：
reports/Series_Overview/figures/fig1_roadmap.png
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(HERE, "figures", "fig1_roadmap.png")


def _register_cjk_font():
    """注册 Noto Sans CJK SC（与各分册 make_figures 相同的逻辑）。"""
    import matplotlib.font_manager as fm
    import tempfile
    from fontTools.ttLib import TTCollection

    for cand in fm.fontManager.ttflist:
        if "Noto Sans CJK SC" in cand.name:
            plt.rcParams["font.sans-serif"] = ["Noto Sans CJK SC", "DejaVu Sans"]
            plt.rcParams["axes.unicode_minus"] = False
            return True
    for path in ("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
                 "/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc"):
        if os.path.exists(path):
            try:
                coll = TTCollection(path)
                for face in coll.fonts:
                    name = face["name"].getDebugName(1)
                    if name and "Noto Sans CJK SC" in name:
                        fd, tmp = tempfile.mkstemp(suffix=".otf")
                        face.save(tmp)
                        os.close(fd)
                        fm.fontManager.addfont(tmp)
                        plt.rcParams["font.sans-serif"] = ["Noto Sans CJK SC", "DejaVu Sans"]
                        plt.rcParams["axes.unicode_minus"] = False
                        return True
            except Exception:
                pass
    return False


def main():
    _register_cjk_font()
    fig, ax = plt.subplots(figsize=(9.2, 5.4))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 62)
    ax.axis("off")

    # 标题
    ax.text(50, 59.5, "一维半导体器件数值仿真系列：复杂度逐级递进路线", ha="center",
            va="center", fontsize=13, fontweight="bold", color="#1a1a1a")

    nodes = [
        # (x, 标题, 副题1, 副题2, 颜色, 是否实线框)
        (13, "S1　MOS 电容", "非线性 Poisson\n（Boltzmann 统计）", "非均匀网格 + 阻尼 Newton\n+ 残差底噪豁免", "#2c3e50", True),
        (37, "S2　PN 结", "Poisson + 漂移-扩散\n全耦合（SG + QFP）", "Gummel 外层迭代\n块 Thomas + SG 离散", "#2c3e50", True),
        (63, "S3　MOSFET（GCA）", "单截面准平衡 Poisson\n+ 渐变沟道近似电流积分", "预解表格 + 热启动\nLAPACK 带状求解", "#2c3e50", True),
        (88, "S4　全 DD（预留）", "二维/一维全漂移-扩散\n（含量子修正预留）", "接口已定义\n实现待续", "#9aa5b1", False),
    ]
    bw, bh = 21, 22
    for x, title, l1, l2, color, solid in nodes:
        box = FancyBboxPatch((x - bw / 2, 26), bw, bh,
                             boxstyle="round,pad=0.6,rounding_size=1.6",
                             linewidth=2.0 if solid else 1.4,
                             edgecolor=color if solid else "#9aa5b1",
                             facecolor="#f4f7fa" if solid else "#f7f7f7",
                             linestyle="-" if solid else "--")
        ax.add_patch(box)
        ax.text(x, 44.2, title, ha="center", va="center", fontsize=10.5,
                fontweight="bold", color=color)
        ax.text(x, 37.6, l1, ha="center", va="center", fontsize=8.4, color="#333")
        ax.text(x, 30.6, l2, ha="center", va="center", fontsize=8.0, color="#555")

    # 递进箭头
    for x0, x1 in [(24.5, 25.5), (48.5, 51.5), (74.5, 77.5)]:
        pass
    arrow_pairs = [(24.3, 26.2), (48.3, 52.2), (74.3, 77.2)]
    for x0, x1 in arrow_pairs:
        ax.add_patch(FancyArrowPatch((x0, 37), (x1, 37), arrowstyle="-|>",
                                     mutation_scale=16, linewidth=1.8,
                                     color="#2c3e50"))

    # 闭环退化标注（箭头下方）
    loops = [
        (38, "阶段闭环：S2 平衡退化\n（V_a=0 → V_bi、J=0）"),
        (63, "阶段闭环：S3 V_D=0\n退化 S1 MOS C-V"),
        (88, "（预留）S4 V_G=0\n退化 S2"),
    ]
    for x, txt in loops:
        ax.text(x, 19.5, txt, ha="center", va="center", fontsize=7.8, color="#7a6410")

    # 公共层横带
    band = FancyBboxPatch((6, 4), 88, 8, boxstyle="round,pad=0.4,rounding_size=1.2",
                          linewidth=1.4, edgecolor="#7f8c8d", facecolor="#eef2f6")
    ax.add_patch(band)
    ax.text(50, 8.9, "common/ 公共层：constants（唯一参数源）· grid · thomas（三对角/带状/块）· bernoulli（SG 核）· qfp · validation（四分类）· io",
            ha="center", va="center", fontsize=8.0, color="#2c3e50")
    ax.text(50, 5.6, "统一 run() 契约：params / results / validation / defects；golden 基准 rtol = 1e-6；pytest 66 项全部通过",
            ha="center", va="center", fontsize=7.6, color="#555")

    # 从各节点到公共层的虚线
    for x in (13, 37, 63, 88):
        ax.plot([x, x], [26, 12.6], linestyle=":", linewidth=1.0, color="#95a5a6")

    fig.savefig(OUT, dpi=200, bbox_inches="tight", facecolor="white")
    print("roadmap ->", OUT)


if __name__ == "__main__":
    main()
