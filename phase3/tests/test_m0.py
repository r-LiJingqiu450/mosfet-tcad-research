"""Package B-2 regression tests for the parameter-injection layer (phase3/m0.py)."""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import pytest

from phase3 import m0
from semiconductor_sim.stage3_mosfet import reference as prod

P = m0.M0Params.f0()


def test_f0_matches_constants():
    from common import constants as C
    assert P.cox == C.COX and P.na == C.NA and P.L == C.L_CH and P.mu_n == C.MU_N
    assert abs(P.psi_b0 - (-0.357159)) < 5e-7          # R05 PM checklist value


@pytest.mark.parametrize('vg,pn', [(.05, 0.), (.5, .05), (1.25, .5), (2., 1.5)])
def test_slice_bitwise(vg, pn):
    a = prod.reference_slice(vg, pn); b = m0.reference_slice(vg, pn, P)
    for k in ('phi_s', 'qs', 'ns', 'ns_excess', 'ns_background'):
        assert a[k] == b[k], k


@pytest.mark.parametrize('vg,vd', [(.5, .05), (1.25, .2), (2., 1.5)])
def test_current_b2_tolerance(vg, vd):
    a = prod.reference_current(vg, vd)['ID']; b = m0.reference_current(vg, vd, P)['ID']
    assert abs(a - b) <= 1e-13 * abs(a)


def test_injection_changes_result():
    p2 = m0.M0Params(**{**{k: v for k, v in P.to_dict().items()
                           if k not in ('cox', 'p0', 'psi_b0')}, 'tox': 10e-9})
    assert m0.reference_current(1.25, .05, p2)['ID'] > m0.reference_current(1.25, .05, P)['ID']


def test_conductance_is_vd_to_zero_limit():
    g = m0.conductance0(1.25, P)
    vd = 1e-5
    i = m0.current_per_um(1.25, vd, P)
    assert abs(i / (vd * g) - 1) < 2 * vd / P.vt     # O(VD/(2VT)) term


def test_threshold_finders_agree():
    vt, st = m0.vt_cc_direct(.05, P)
    assert st == 'ok' and 1.0 < vt < 1.1
    import numpy as np
    vg = np.arange(vt - .02, vt + .021, .005)
    i = [m0.current_per_um(v, .05, P) for v in vg]
    vt2, st2 = m0.vt_cc_tabulated(vg, i, m0.ith_per_um(P))
    assert st2 == 'ok' and abs(vt2 - vt) < 1e-5


def test_tabulated_reports_status():
    assert m0.vt_cc_tabulated([0, 1], [1e-9, 1e-8], 1e-7)[1] == 'not_bracketed'
    assert math.isnan(m0.vt_cc_tabulated([0, 1], [1e-9, 1e-8], 1e-7)[0])
