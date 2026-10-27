# -*- coding: utf-8 -*-
"""
一维半导体器件数值仿真系列（semiconductor_sim）。

统一仓库入口：三阶段共享 common/ 公共模块，各阶段独立成册。

    python -m semiconductor_sim.stage1_mos.run        # S1 MOS 电容
    python -m semiconductor_sim.stage2_pn.run         # S2 PN 结
    python -m semiconductor_sim.stage3_mosfet.run     # S3 MOSFET GCA

目录约定（docs/spec/code_spec.md）：
    common/       七模块公共层（constants/grid/thomas/bernoulli/qfp/validation/io）
    stageX_*/     各阶段求解器、run() 主驱动、make_figures() 绘图、tests/
    docs/spec/    六类系列规范（M0 冻结）
    results/      阶段结果（stageX/）与 golden 基准（golden/）
    reports/      各分册报告与系列总览
"""

import os as _os
import sys as _sys

# 导入兼容层：各阶段模块沿用原件的仓库内导入风格（from common.xxx import
# ...），而系列测试/脚本使用包绝对导入（from semiconductor_sim.xxx import
# ...）。包被导入时把仓库根目录插入 sys.path，两种风格均可解析——
# 不修改任何阶段源码的 import 行，保持"拆分保等价"的行级忠实。
_ROOT = _os.path.dirname(_os.path.abspath(__file__))
if _ROOT not in _sys.path:
    _sys.path.insert(0, _ROOT)

__version__ = "1.2.0-phase1"
