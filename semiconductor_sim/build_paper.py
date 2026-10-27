#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_paper.py —— 系列报告统一构建入口（doc_template_spec.md 第 16 节）。

输入：reports/Sx_*/body.html（16 节正文片段）+ reports/Sx_*/figures/（图）
输出：reports/Sx_*/<名称>.pdf

    一维半导体器件数值仿真系列
    ├── reports/S1_MOS_Capacitor/body.html      -> S1_MOS_Capacitor.pdf
    ├── reports/S2_PN_Junction/body.html        -> S2_PN_Junction.pdf
    ├── reports/S3_MOSFET_GCA/body.html         -> S3_MOSFET_GCA.pdf
    └── reports/Series_Overview/body.html       -> Series_Overview.pdf

用法（仓库根目录）：
    python build_paper.py                # 全部四份
    python build_paper.py S1             # 只建 S1（S1/S2/S3/Overview）

约束（reproduce_spec.md 第 4 节）：报告构建脚本不得内联重算任何数字；
正文数字必须与 results/stageX/terminal.txt 同源。
"""

import os
import subprocess
import sys

from weasyprint import HTML

BASE = os.path.dirname(os.path.abspath(__file__))
REPORTS = os.path.join(BASE, "reports")

# ============================================================================
# 系列统一版式（源自三份原件报告 CSS，按 figure_spec/doc_template_spec 统一）
# ============================================================================
CSS = """
@page {
    size: A4;
    margin: 2.4cm 2.2cm 2.2cm 2.2cm;
    @bottom-center { content: counter(page); font-family: "Noto Serif CJK SC"; font-size: 9pt; color: #555; }
}
@page cover { @bottom-center { content: none; } }

body { font-family: "Noto Serif CJK SC", serif; font-size: 10.5pt; line-height: 1.75; color: #1a1a1a; }
.cover { page: cover; text-align: center; page-break-after: always; padding-top: 4.0cm; }
.cover .series { font-size: 12pt; letter-spacing: 3px; color: #2c3e50; margin-bottom: 0.8cm;
                 font-family: "Noto Sans CJK SC", sans-serif; }
.cover .book-id { display: inline-block; font-family: "Noto Sans CJK SC", sans-serif;
                  font-size: 11pt; color: #fff; background: #2c3e50; padding: 2px 14px;
                  border-radius: 3px; margin-bottom: 1.0cm; letter-spacing: 2px; }
.cover h1 { font-family: "Noto Sans CJK SC", sans-serif; font-size: 21pt; line-height: 1.5; margin: 0.6cm 0 0.3cm; }
.cover .subtitle { font-size: 12pt; color: #444; margin-bottom: 2.2cm; }
.cover .meta { font-size: 11pt; line-height: 2.1; }
.cover .meta .blank { display: inline-block; border-bottom: 1px solid #333; min-width: 6.5cm; height: 0.75em; }
.cover .footer-note { margin-top: 2.8cm; font-size: 9.5pt; color: #777; }

h2 { font-family: "Noto Sans CJK SC", sans-serif; font-size: 14pt; border-left: 4px solid #2c3e50;
     padding-left: 8px; margin: 1.05em 0 0.5em; page-break-after: avoid; }
h3 { font-family: "Noto Sans CJK SC", sans-serif; font-size: 11.5pt; color: #2c3e50;
     margin: 0.9em 0 0.4em; page-break-after: avoid; }
p { margin: 0.45em 0; text-align: justify; }

.eq { text-align: center; margin: 0.65em 0; font-size: 10.5pt; page-break-inside: avoid; }
.eq .num { position: absolute; right: 0; color: #555; }
.eqwrap { position: relative; }

figure { margin: 1.0em 0; text-align: center; page-break-inside: avoid; }
figure img { width: 92%; max-width: 15.5cm; }
figcaption { font-size: 9.5pt; color: #333; margin-top: 0.35em; }
figcaption b { color: #000; }

table { border-collapse: collapse; margin: 0.8em auto; width: 96%; font-size: 9.5pt; page-break-inside: avoid; }
th { background: #2c3e50; color: #fff; padding: 5px 8px; font-family: "Noto Sans CJK SC", sans-serif; }
td { border: 1px solid #bbb; padding: 4px 8px; text-align: center; }
td.l { text-align: left; }
td.nw { white-space: nowrap; }
tr:nth-child(even) td { background: #f4f6f8; }

.abstract { margin: 0.2em 0 0.8em 0; }
.abstract .kw { margin-top: 0.5em; }
.abstract-title { font-family: "Noto Sans CJK SC", sans-serif; font-weight: bold; font-size: 11pt; }

.ref ol { padding-left: 1.6em; font-size: 9.8pt; }
.ref li { margin-bottom: 0.35em; }

.appendix-title { page-break-before: always; }

pre.code, pre { font-family: "DejaVu Sans Mono", monospace; font-size: 7.6pt; line-height: 1.32;
    background: #fafafa; border: 1px solid #d5d5d5; border-radius: 3px;
    padding: 7px 9px; white-space: pre-wrap; overflow-wrap: break-word;
    page-break-inside: auto; margin: 0.5em 0; }
code { font-family: "DejaVu Sans Mono", monospace; font-size: 9pt; background: #f0f0f0; padding: 0 3px; border-radius: 2px; }

.callout { border-left: 4px solid #b8860b; background: #fdf6e3; padding: 6px 12px;
           margin: 0.7em 0; font-size: 9.8pt; }
.callout.blue { border-left-color: #2c3e50; background: #eef2f6; }
.series-note { border: 1px solid #c9d4de; background: #f4f7fa; border-radius: 4px;
               padding: 8px 14px; margin: 0.8em 0; font-size: 10pt; }
.filecap { font-family: "Noto Sans CJK SC", sans-serif; font-size: 9pt; font-weight: bold;
           color: #2c3e50; margin: 0.9em 0 0.25em; border-bottom: 1px solid #c9d4de;
           padding-bottom: 2px; }
"""

REPORT_DIRS = {
    "S1": "S1_MOS_Capacitor",
    "S2": "S2_PN_Junction",
    "S3": "S3_MOSFET_GCA",
    "Overview": "Series_Overview",
}


def build_one(key):
    """构建一份报告：body.html + figures/ -> <目录名>.pdf。缺失则跳过。"""
    d = os.path.join(REPORTS, REPORT_DIRS[key])
    body_path = os.path.join(d, "body.html")
    if not os.path.exists(body_path):
        print(f"[{key}] 跳过：缺少正文 {os.path.relpath(body_path, BASE)}")
        return None
    with open(body_path, encoding="utf-8") as f:
        body = f.read()
    # 可选补充片段：报告目录带 make_appendix_code.py 时（如 S3 的附录 F 完整源码），
    # 先在构建时从仓库源码重新生成，再拼接到正文之后（见该脚本说明）。
    gen = os.path.join(d, "make_appendix_code.py")
    if os.path.exists(gen):
        subprocess.run([sys.executable, gen], check=True, cwd=BASE)
        extra_path = os.path.join(d, "appendix_code.html")
        if os.path.exists(extra_path):
            with open(extra_path, encoding="utf-8") as f:
                body += "\n" + f.read()
    html = (f'<!DOCTYPE html>\n<html lang="zh-CN">\n<head><meta charset="utf-8">'
            f'<style>{CSS}</style></head>\n<body>\n{body}\n</body>\n</html>')
    out = os.path.join(d, REPORT_DIRS[key] + ".pdf")
    HTML(string=html, base_url=d).write_pdf(out)
    return out


def main(argv):
    keys = argv[1:] or list(REPORT_DIRS)
    for k in keys:
        if k not in REPORT_DIRS:
            print(f"未知报告: {k}；可选: {list(REPORT_DIRS)}")
            return 2
        out = build_one(k)
        if out:
            print(f"[{k}] -> {os.path.relpath(out, BASE)}")
    print("完成。")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
