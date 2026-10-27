# -*- coding: utf-8 -*-
"""pytest 路径引导：使 semiconductor_sim 包可从仓库父目录导入。

运行方式（在仓库父目录）：
    pytest semiconductor_sim/tests -q
"""

import os
import sys

_PARENT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _PARENT not in sys.path:
    sys.path.insert(0, _PARENT)
