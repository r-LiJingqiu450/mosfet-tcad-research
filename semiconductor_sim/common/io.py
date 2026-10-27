# -*- coding: utf-8 -*-
"""
系列统一 IO：结果目录、终端输出留存、golden 基准对比。

规范（docs/spec/reproduce_spec.md）：
    - 每阶段运行输出固定写入 results/stageX/terminal.txt；
    - 关键数字以 JSON/NPZ 存 results/golden/，tests/test_golden.py 比对；
    - 图固定文件名写入 results/stageX/。
"""

import json
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(BASE, "results")
GOLDEN = os.path.join(RESULTS, "golden")


def ensure_results_dir(stage):
    """确保 results/stageX/ 存在并返回路径。stage 为 'stage1' 等。"""
    d = os.path.join(RESULTS, stage)
    os.makedirs(d, exist_ok=True)
    return d


def terminal_path(stage):
    return os.path.join(ensure_results_dir(stage), "terminal.txt")


class TerminalLog:
    """
    终端输出双写：print 到屏幕的同时留存 results/stageX/terminal.txt。

    用法：
        log = TerminalLog("stage1")
        print(...)              # 正常 print 即被捕获
        log.close()
    或 with TerminalLog("stage1") as log: ...

    echo=False 时只写文件不回显（供测试/静默批处理使用）。
    """

    def __init__(self, stage, echo=True):
        self.echo = echo
        self._file = open(terminal_path(stage), "w", encoding="utf-8")
        self._stdout = sys.stdout
        # 构造即接管 sys.stdout：文档承诺的两种用法（裸调用 + with）
        # 都必须生效。早期实现只在 __enter__ 里切换，裸调用时 print
        # 全部漏记（文件为空）——见 docs/defects  DEF-C-01。
        sys.stdout = self

    def write(self, s):
        if self.echo:
            self._stdout.write(s)
        self._file.write(s)

    def flush(self):
        self._stdout.flush()
        self._file.flush()

    def close(self):
        self._file.close()
        sys.stdout = self._stdout

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False


def save_json(name, obj, dirpath=GOLDEN):
    """存 JSON（golden 基准）。name 不含扩展名。"""
    os.makedirs(dirpath, exist_ok=True)
    path = os.path.join(dirpath, name + ".json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=1, default=_jsonable)
    return path


def load_json(name, dirpath=GOLDEN):
    with open(os.path.join(dirpath, name + ".json"), encoding="utf-8") as f:
        return json.load(f)


def _jsonable(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(f"不可 JSON 化: {type(o)}")


def compare_golden(name, current, rtol=1e-6, atol=0.0, dirpath=GOLDEN):
    """
    当前结果与 golden 基准比对。

    current : {key: value} 数值字典（与 save_json 存的结构一致）
    返回 dict：{key: (passed, rel_dev)}；golden 缺失时返回 None 表示跳过。
    """
    try:
        ref = load_json(name, dirpath)
    except FileNotFoundError:
        return None
    out = {}
    for k, v in current.items():
        if k not in ref:
            continue
        rv = float(ref[k])
        cv = float(v)
        if rv == 0.0:
            ok = abs(cv - rv) <= atol
            dev = abs(cv - rv)
        else:
            dev = abs(cv - rv) / abs(rv)
            ok = dev <= rtol
        out[k] = (ok, dev)
    return out


def save_npz(name, **arrays):
    """存 NPZ（golden 数组基准）。"""
    os.makedirs(GOLDEN, exist_ok=True)
    path = os.path.join(GOLDEN, name + ".npz")
    np.savez(path, **arrays)
    return path


def load_npz(name):
    return np.load(os.path.join(GOLDEN, name + ".npz"))
