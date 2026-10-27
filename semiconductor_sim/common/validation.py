# -*- coding: utf-8 -*-
"""
系列统一验证 / 自洽检查工具。

规范（docs/spec/validation_spec.md）：
    验证        有独立解析/独立估计参照，必须给出参照公式、判据、实测、结论
    自洽检查    无独立参照，仅内部一致性（网格/表格收敛等），必须标注"自洽"
    格式自洽    质量作用定律等恒等式（机器精度），属于自洽的子类
    artifact    数值构造导致、非物理，标注 ART-Sx-yy，不得包装成物理效应

Checklist 收集各阶段全部指标，run() 将其放入返回 dict 的
"validation" 键；to_rows() 直接生成论文验证表的三线表行。
"""

from dataclasses import dataclass, field, asdict

import numpy as np

# 指标类型（规范 3.2）
VALIDATION = "验证"
SELF_CHECK = "自洽"
FORMAT_CHECK = "格式自洽"
ARTIFACT = "artifact"


@dataclass
class Check:
    """单条验证/自洽检查记录。"""
    metric: str                 # 指标名（如 "C-V 解析对比"）
    kind: str                   # VALIDATION / SELF_CHECK / FORMAT_CHECK / ARTIFACT
    reference: str              # 参照公式或基准描述
    criterion: str              # 判据（如 "< 0.1%"）
    measured: str = ""          # 实测值（格式化后的字符串）
    passed: object = None       # True / False / None（artifact 无判据）
    note: str = ""              # 备注（物理归属、已知边界说明等）
    artifact: str = ""          # ART-Sx-yy 编号（仅 artifact 类）


@dataclass
class Checklist:
    """一个阶段的全部检查记录。"""
    stage: str
    checks: list = field(default_factory=list)

    def add(self, **kw):
        # passed 归一化为 Python bool：调用方传入的判据多是 numpy 标量
        # 比较结果（np.bool_），而 all_passed/failures/to_rows 用
        # `is True / is False` 身份判断——np.bool_ 会同时逃过两者，
        # 导致失败项被漏判、通过项显示"—"（docs/defects DEF-C-02）。
        p = kw.get("passed")
        if p is not None and not isinstance(p, bool):
            kw["passed"] = bool(p)
        self.checks.append(Check(**kw))
        return self.checks[-1]

    def all_passed(self):
        """全部有判据的检查是否通过（artifact 与 None 不参与）。"""
        return all(c.passed is None or bool(c.passed) for c in self.checks)

    def failures(self):
        return [c for c in self.checks
                if c.passed is not None and not bool(c.passed)]

    def to_rows(self):
        """生成论文验证表行：[指标, 类型, 参照, 判据, 实测, 结论]。"""
        rows = []
        for c in self.checks:
            if c.passed is None:
                concl = "—"
            elif bool(c.passed):
                concl = "通过"
            else:
                concl = "未达标"
            rows.append([c.metric, c.kind, c.reference, c.criterion,
                         c.measured, concl + (("；" + c.note) if c.note else "")])
        return rows


# ================================================================================
# 数值工具
# ================================================================================
def rel_dev(a, b):
    """相对偏差 (a-b)/b（b 为参照）。"""
    return (np.asarray(a, dtype=float) - np.asarray(b, dtype=float)) \
        / np.asarray(b, dtype=float)


def max_rel_dev(a, b, mask=None):
    """
    数组 a 相对参照 b 的最大相对偏差（返回 (最大值, 位置索引)）。
    mask 为 True 的位置参与统计。
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    d = np.abs(rel_dev(a, b))
    if mask is not None:
        d = np.where(mask, d, -np.inf)
    i = int(np.nanargmax(d))
    return float(d[i]), i


def format_sci(x, digits=3):
    """科学计数法格式化（论文表格用）。"""
    if isinstance(x, str):
        return x
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "—"
    return f"{x:.{digits}e}"


def format_pct(x, digits=2):
    """百分比格式化。"""
    if isinstance(x, str) or x is None:
        return "—"
    return f"{x * 100:.{digits}f}%"
