# -*- coding: utf-8 -*-
"""
Scharfetter–Gummel 电流离散的 Bernoulli 插值核（S2 使用，S3 预留）。

B(x) = x / (1 - e^{-x})，B(0) = 1。
数值稳定实现：小 |x| 泰勒展开，大 |x| 取极限。
"""

import numpy as np


def bernoulli(x):
    """
    Bernoulli 函数 B(x) = x / (1 - e^{-x})，B(0) = 1。

    SG 电流格式的核心插值核。数值稳定实现：
        |x| < 1e-8  : 泰勒展开 B = 1 + x/2 + x^2/12
        大正数 x    : e^{-x} 下溢为 0 -> B = x（正确极限）
        大负数 x    : e^{|x|} 可达 inf -> B = x/(1-inf) = -0.0（正确极限 ~0）
    """
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = np.abs(x) < 1e-8
    xs = x[small]
    out[small] = 1.0 + 0.5 * xs + xs * xs / 12.0
    xl = x[~small]
    with np.errstate(over="ignore"):
        out[~small] = xl / (1.0 - np.exp(-xl))
    return out


def bernoulli_deriv(x):
    """
    Bernoulli 函数导数 B'(x)，用于全 Newton 格式中电流对电势的耦合项。

    恒等式（推导自 B(x) = x/(1-e^{-x})）：
        B'(x) = B(x) (1 - B(-x)) / x
    验证：x=0 时 B'(0)=1/2；x->+inf 时 B'->1；x->-inf 时 B'->0。
    小 |x| 用泰勒 B' = 1/2 + x/6 + x^2/24。
    """
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = np.abs(x) < 1e-6
    xs = x[small]
    out[small] = 0.5 + xs / 6.0 + xs * xs / 24.0
    xl = x[~small]
    out[~small] = bernoulli(xl) * (1.0 - bernoulli(-xl)) / xl
    return out
