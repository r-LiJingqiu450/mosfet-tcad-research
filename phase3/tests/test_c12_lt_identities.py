"""C-12: identities and limits inside the Liang-Taur double-gate model (IEEE T-ED 51(8) 1385, 2004; Liang, UCSD 2006).

- lambda_1 from LT04 eq. (8) reproduces the values quoted in the dissertation (Figs. 4.2(c), 4.5(a)-(c)).
- The published conjugate function g1 equals eps(x) u_L1(x,0)/eps_si (eps-weighted eigenfunction).
- The closed form b1 (LT04 eq. 12) equals the projection eq. (11) of the LT boundary data.
- The long-channel, V_ds -> 0 end correction of eqs. (13)+(16) tends to (lam1/pi) Ein(b1/V_T) in the thin-film limit,
  and to (lam1/pi)[ln(b1/V_T) + C(s)] with C(0) = Euler gamma.
- Rows (Thomson) <= slice (Dirichlet = LT04 eq. 16) for the frozen LT potential.
"""
from __future__ import annotations

import math

import pytest

from phase3.c12 import lt_endcorr as lt
from phase3.c12 import lt_vt_rolloff as vt


@pytest.mark.parametrize('name', ['G2_tsi5_ti1.5_SiO2', 'G3_tsi25_ti1.5_SiO2', 'G4_tsi25_ti3_k7.6', 'G5_tsi25_ti13.5_k35.1'])
def test_lambda1_reproduces_quoted_values(name):
    tsi, ti, ei = lt.GEOMS[name]
    assert abs(lt.lam1(tsi, ti, ei) - lt.QUOTED_LAMBDA1[name]) < 0.1


@pytest.mark.parametrize('name', list(lt.GEOMS))
def test_conjugate_function_is_eps_weighted_eigenfunction(name):
    assert lt.conjugate_identity_residual(*lt.GEOMS[name]) < 1e-12


@pytest.mark.parametrize('name', list(lt.GEOMS))
def test_b1_closed_form_equals_projection(name):
    tsi, ti, ei = lt.GEOMS[name]
    cf = lt.b1_closed_form(tsi, ti, ei)
    assert abs(lt.project_lt_data(math.pi / lt.lam1(tsi, ti, ei), tsi, ti, ei) / cf - 1) < 1e-12
    assert abs(lt.b1_quadrature_eq11(tsi, ti, ei) / cf - 1) < 1e-10


def test_thin_film_limit_is_ein():
    for beta in (5.0, 15.0, 30.0):
        assert abs(lt.F_slice(beta, 1e-3) / lt.ein(beta) - 1) < 1e-6
        assert abs(lt.F_rows(beta, 1e-3) / lt.ein(beta) - 1) < 1e-6


def test_asymptotic_constant():
    assert abs(lt.C_const(1e-4) - lt.GAMMA) < 1e-6
    cs = [lt.C_const(s) for s in (0.2, 0.5, 0.8, 1.2)]
    assert all(a > b for a, b in zip(cs, cs[1:]))            # finite film thickness lowers C(s) below gamma
    s, beta = 0.886, 400.0
    assert abs(lt.F_slice(beta, s) - (math.log(beta) + lt.C_const(s))) < 1e-3


def test_rows_below_slice_and_finite_L_limit():
    s, beta = 0.886, 15.0
    assert lt.F_rows(beta, s) < lt.F_slice(beta, s)
    assert abs(lt.dL_finite_L(beta, s, 20.0) / lt.F_slice(beta, s) - 1) < 1e-6
    l1 = lt.lam1(10.0, 1.5, lt.EPS_OX)
    for L in (24.0, 60.0, 200.0):
        r = vt.rows_over_slice(L, 0.43, 1.56, l1, math.pi * 10.0 / (2 * l1))
        assert 0.98 < r <= 1.0


def test_higher_modes_negligible_for_end_correction():
    tsi, ti, ei = lt.GEOMS['G1_tsi10_ti1.5_SiO2']
    l1 = lt.lam1(tsi, ti, ei)
    s = math.pi * tsi / (2 * l1)
    A0 = 0.4
    one = l1 / math.pi * lt.F_slice(lt.b1_closed_form(tsi, ti, ei, A0) / lt.VT, s)
    full = lt.dL_full_series(tsi, ti, ei, A0, nmodes=150)['dL_nm']
    assert abs(full / one - 1) < 1e-3
