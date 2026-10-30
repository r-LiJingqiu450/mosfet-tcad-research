"""Unit tests for the C-3 analysis helper (R04 7.4 classification with actual ratios)."""
import math

from phase3.c3.t1_analysis import classify


def _seq(c, a, p, h):
    return [c + a * x ** p for x in h]


def test_asymptotic_exact_power_law_nonuniform_ratios():
    h = [1.0, 1.4147, 1.4147 * 1.4137, 1.4147 * 1.4137 * 1.4124]     # finest first, actual-style ratios
    f = _seq(2.0, 1e-3, 2.0, h)
    out = classify(f, h)
    assert out['cls'] == 'asymptotic'
    assert abs(out['p123'] - 2.0) < 1e-9 and abs(out['p234'] - 2.0) < 1e-9
    assert abs(out['richardson'] - 2.0) < 1e-12
    # U = 1.25 |d12| / (r12^p - 1) equals 1.25 x the true error of the finest value here
    assert math.isclose(out['U'], 1.25 * 1e-3, rel_tol=1e-9)


def test_nominal_ratio_would_bias_order():
    h = [1.0, 1.4147, 1.4147 * 1.4137, 1.4147 * 1.4137 * 1.4124]
    f = _seq(0.0, 1.0, 2.0, h)
    nominal = [math.sqrt(2) ** k for k in range(4)]
    assert abs(classify(f, h)['p123'] - 2.0) < 1e-9
    assert abs(classify(f, nominal)['p123'] - 2.0) > 1e-3


def test_oscillatory_and_roundoff():
    h = [1.0, 2.0, 4.0, 8.0]
    assert classify([1.0, 1.1, 0.95, 1.2], h)['cls'] == 'oscillatory'
    assert classify([1.0, 1.0, 1.0, 1.3], h)['cls'] == 'converged_to_roundoff'


def test_monotone_nonasymptotic():
    h = [1.0, 2.0, 4.0, 8.0]
    f = [1.0, 1.0 + 1e-3, 1.0 + 1e-3 + 2.5e-3, 1.0 + 1e-3 + 2.5e-3 + 1e-1]   # p123 ~ 1.3, p234 ~ 5.3
    out = classify(f, h)
    assert out['cls'] == 'monotone_nonasymptotic'
    assert out['U'] == 3 * max(abs(x - f[0]) for x in f)
