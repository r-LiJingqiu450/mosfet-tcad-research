"""C-10 numerical issue E10-1: the zero-bias network conductance must be taken from the near-zero (source) side.

Test network: a series chain with huge conductances next to both contacts and one tiny link in the middle, mimicking
deep subthreshold.
- The exact conductance is the series sum.
- The drain-side current, sum c (1 - phi) with phi ~ 1 - 1e-40, is lost to cancellation.
- The source-side current, sum c phi with phi ~ 1e-40, keeps full relative accuracy.
"""
from __future__ import annotations

import numpy as np
import pytest


def chain(cs):
    n = len(cs) + 1
    return n, np.arange(n - 1), np.arange(1, n), np.asarray(cs, float)


def test_source_side_conductance_is_cancellation_free():
    pytest.importorskip('devsim')
    from phase3.c6.wprobe import solve_net
    from phase3.c10.probe10 import solve_net_robust
    cs = [1e30] * 5 + [1e-10] + [1e30] * 5
    n, a, b, c = chain(cs)
    exact = 1.0 / np.sum(1.0 / c)
    G_drain, I_src, _, _ = solve_net(n, a, b, c, np.array([0]), np.array([n - 1]))
    G_robust, _, _, _ = solve_net_robust(n, a, b, c, np.array([0]), np.array([n - 1]))
    assert abs(G_robust / exact - 1) < 1e-12
    assert abs(-I_src / exact - 1) < 1e-12
    assert not abs(G_drain / exact - 1) < 1e-3          # the drain-side value is wrong (documents E10-1)


def test_benign_network_both_sides_agree():
    pytest.importorskip('devsim')
    from phase3.c6.wprobe import solve_net
    from phase3.c10.probe10 import solve_net_robust
    n, a, b, c = chain([1.0, 2.0, 3.0, 4.0])
    exact = 1.0 / np.sum(1.0 / c)
    G_drain, I_src, _, _ = solve_net(n, a, b, c, np.array([0]), np.array([n - 1]))
    G_robust, _, _, _ = solve_net_robust(n, a, b, c, np.array([0]), np.array([n - 1]))
    for g in (G_drain, -I_src, G_robust):
        assert abs(g / exact - 1) < 1e-12
