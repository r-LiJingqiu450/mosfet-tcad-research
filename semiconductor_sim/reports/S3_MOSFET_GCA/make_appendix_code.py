#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_appendix_code.py —— 生成 S3 分册附录 F（完整源码）HTML 片段。

按依赖顺序读取仓库源码，HTML 转义后写入同目录 appendix_code.html；
build_paper.py 构建 S3 时自动调用本脚本，保证 PDF 中的代码与仓库严格一致。

收录范围（相对路径相对于仓库根 semiconductor_sim/）：
  common/        公共层 7 个模块（bernoulli/qfp 为 S2 传承、由 test_common.py 守护，一并收录）
  stage3_mosfet/ 第三阶段实现 6 个模块 + 单元测试
  tests/         系列测试 4 个文件（golden 基准、跨阶段闭环、公共层测试）
各包 __init__.py 为空文件，不收录。
"""
import html
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))          # reports/S3_MOSFET_GCA -> 仓库根
OUT = os.path.join(HERE, "appendix_code.html")

FILES = [
    # ---- common/ 公共层（S1/S2/S3 共享）----
    "common/constants.py",
    "common/grid.py",
    "common/thomas.py",
    "common/bernoulli.py",
    "common/qfp.py",
    "common/validation.py",
    "common/io.py",
    # ---- stage3_mosfet/ 第三阶段实现 ----
    "stage3_mosfet/solver.py",
    "stage3_mosfet/table.py",
    "stage3_mosfet/gca.py",
    "stage3_mosfet/postprocess.py",
    "stage3_mosfet/run.py",
    "stage3_mosfet/make_figures.py",
    # ---- 测试 ----
    "stage3_mosfet/tests/test_stage3.py",
    "tests/golden_helpers.py",
    "tests/test_golden.py",
    "tests/test_stage_closure.py",
    "tests/test_common.py",
]

HEAD = """<!-- ================= 附录 F 完整源码（文末补充，构建时自动生成） ================= -->
<h2 class="appendix-title">附录 F　完整源码（全文收录）</h2>
<p>应读者要求，本附录按依赖顺序收录本分册数值链路与测试套件的全部源文件
全文：<code>common/</code> 公共层（constants、grid、thomas、bernoulli、qfp、
validation、io；其中 bernoulli 与 qfp 为第二阶段传承、由
<code>tests/test_common.py</code> 守护的公共组件，一并收录）、
<code>stage3_mosfet/</code> 第三阶段实现（solver、table、gca、postprocess、
run、make_figures）及其单元测试、<code>tests/</code> 系列测试（golden 基准
比对、跨阶段闭环、公共层测试），共 {nfiles} 个文件、{nlines} 行（相对路径
相对于仓库根目录 semiconductor_sim/；各包 <code>__init__.py</code> 为空文件，
未收录）。本附录由 <code>reports/S3_MOSFET_GCA/make_appendix_code.py</code>
在报告构建时从仓库源码自动生成，与附录 C 的行号指向、附录 D 的复现命令
一一对应；源码的任何改动都会在下一次构建时同步到本附录。</p>
"""


def main():
    parts = []
    total = 0
    for rel in FILES:
        path = os.path.join(ROOT, rel)
        if not os.path.exists(path):
            raise SystemExit(f"缺少源文件: {rel}")
        with open(path, encoding="utf-8") as f:
            text = f.read()
        n = text.count("\n") + (0 if text.endswith("\n") else 1)
        total += n
        parts.append(
            f'<div class="filecap">{rel}（全文 {n} 行）</div>\n'
            f'<pre class="code">{html.escape(text.rstrip(chr(10)))}</pre>'
        )
    doc = HEAD.format(nfiles=len(FILES), nlines=total) + "\n".join(parts) + "\n"
    with open(OUT, "w", encoding="utf-8") as f:
        f.write(doc)
    print(f"附录 F 已生成: {os.path.relpath(OUT, ROOT)}"
          f"（{len(FILES)} 个文件，{total} 行）")


if __name__ == "__main__":
    main()
