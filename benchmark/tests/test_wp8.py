"""WP8 instruments: double-double arithmetic, double-double residual, iterative refinement, grid-convergence procedure.

They check the instruments on manufactured inputs against exact rational arithmetic (fractions.Fraction). They are
not results.
"""
from __future__ import annotations

import os
import sys
from fractions import Fraction

import numpy as np
import pytest
import scipy.sparse as sp
import scipy.sparse.linalg as spla

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, 'benchmark', 'wp8'))
import ddlinalg as DD  # noqa: E402


def F(x):
    return Fraction(float(x))


def test_error_free_transformations_are_exact():
    rng = np.random.default_rng(1)
    a = rng.standard_normal(2000) * 10.0 ** rng.integers(-150, 150, 2000)
    b = rng.standard_normal(2000) * 10.0 ** rng.integers(-150, 150, 2000)
    s, e = DD.two_sum(a, b)
    p, f = DD.two_prod(a, b)
    for i in range(0, 2000, 7):
        assert F(s[i]) + F(e[i]) == F(a[i]) + F(b[i])
        assert F(p[i]) + F(f[i]) == F(a[i]) * F(b[i])


def _random_sparse(n, density, rng, spread=12):
    A = sp.random(n, n, density=density, random_state=np.random.RandomState(7), format='csr')
    A.data = rng.standard_normal(A.nnz) * 10.0 ** rng.uniform(-spread, spread, A.nnz)
    A = A + sp.eye(n, format='csr') * 1e3
    return A.tocsr()


def test_dd_residual_matches_exact_rational_residual():
    rng = np.random.default_rng(2)
    n = 40
    A = _random_sparse(n, 0.15, rng)
    xh = rng.standard_normal(n) * 10.0 ** rng.uniform(-8, 8, n)
    xl = xh * rng.uniform(-1, 1, n) * 2.0 ** -60
    b = A @ xh + rng.standard_normal(n) * 1e-20
    rh, rl = DD.DDResidual(A)(b, xh, xl)
    Ad = A.toarray()
    for i in range(n):
        exact = F(b[i]) - sum(F(Ad[i, j]) * (F(xh[j]) + F(xl[j])) for j in range(n) if Ad[i, j] != 0)
        scale = sum(abs(F(Ad[i, j]) * F(xh[j])) for j in range(n)) + abs(F(b[i]))
        assert abs((F(rh[i]) + F(rl[i])) - exact) <= Fraction(1, 10 ** 29) * scale


def _exact_solve(Ad, b):
    """Gaussian elimination in rational arithmetic (small n only)."""
    n = len(b)
    M = [[F(Ad[i, j]) for j in range(n)] + [F(b[i])] for i in range(n)]
    for c in range(n):
        piv = max(range(c, n), key=lambda r: abs(M[r][c]))
        M[c], M[piv] = M[piv], M[c]
        for r in range(c + 1, n):
            if M[r][c] != 0:
                f = M[r][c] / M[c][c]
                M[r] = [M[r][k] - f * M[c][k] for k in range(n + 1)]
    x = [Fraction(0)] * n
    for i in range(n - 1, -1, -1):
        x[i] = (M[i][n] - sum(M[i][k] * x[k] for k in range(i + 1, n))) / M[i][i]
    return x


def test_refinement_recovers_a_tiny_output_component_that_double_lu_misses():
    """A tiny solution component (the analogue of a low current) is lost by a double-precision solve and recovered
    by refinement with double-double residuals, against the exact solution of the same double-precision system."""
    rng = np.random.default_rng(3)
    n = 24
    Q1, _ = np.linalg.qr(rng.standard_normal((n, n)))
    Q2, _ = np.linalg.qr(rng.standard_normal((n, n)))
    A = Q1 @ np.diag(np.logspace(0, -9, n)) @ Q2           # condition number 1e9
    x_true = rng.standard_normal(n)
    x_true[0] = 1e-11                                       # tiny output component
    b = A @ x_true
    lu = spla.splu(sp.csc_matrix(A))
    x0 = lu.solve(b)
    exact = _exact_solve(A, b)
    G_exact = float(exact[0])
    out = DD.refine(sp.csr_matrix(A), b, lu.solve, x0, out_index=0, sign=1.0)
    assert abs(x0[0] / G_exact - 1) > 1e-3                  # the double solve misses the tiny component
    assert out['converged']
    assert abs(out['G'] / G_exact - 1) < 1e-10              # refinement recovers it


def test_three_grid_procedure_recovers_a_manufactured_order_and_limit():
    import gci_wp8 as GC
    f_exact, C, p = 3.5e-12, 2.0e-9, 2.0
    h1, r21, r32 = 1.0e-3, 1.41, 1.43                       # unequal ratios exercise q(p)
    h2, h3 = h1 * r21, h1 * r21 * r32
    f = [f_exact + C * h ** p for h in (h1, h2, h3)]
    g = GC.three_grid(f[0], f[1], f[2], r21, r32)
    assert g['monotone'] and abs(g['p'] - p) < 1e-9 and abs(g['f_ext'] / f_exact - 1) < 1e-12
    assert abs(g['asymptotic_ratio'] - 1) <= 2 * g['e_a21']      # equals f1/f2 for an exact power law
    osc = GC.three_grid(1.0, 1.1, 0.95, 1.41, 1.41)
    assert not osc['monotone']


def test_refinement_reports_non_convergence_when_the_factorization_is_too_poor():
    rng = np.random.default_rng(4)
    n = 12
    A = rng.standard_normal((n, n))
    b = rng.standard_normal(n)
    bad = lambda r: 1e3 * r                                  # noqa: E731  (a useless inner solver)
    out = DD.refine(sp.csr_matrix(A), b, bad, np.zeros(n), out_index=0, sign=1.0, max_it=8)
    assert not out['converged']


def test_vector_refinement_and_adjoint_reach_the_exact_solutions():
    """xsens.refine_to (all components, stopped on the componentwise backward error) and its transposed use give the
    exact solutions of A x = b and A^T y = e of the same double-precision system, and y[k] = x[k] when b = e = e_k."""
    os.environ.setdefault('DEVSIM_MATH_LIBS', '/usr/lib/x86_64-linux-gnu/libblas.so.3:/usr/lib/x86_64-linux-gnu/liblapack.so.3')
    try:
        import xsens as XS                                  # imports DEVSIM (the deck modules) as well
    except (ImportError, RuntimeError) as exc:              # pragma: no cover - environment without DEVSIM
        pytest.skip(f'xsens needs DEVSIM: {exc}')
    rng = np.random.default_rng(5)
    n = 20
    Q1, _ = np.linalg.qr(rng.standard_normal((n, n)))
    Q2, _ = np.linalg.qr(rng.standard_normal((n, n)))
    A = Q1 @ np.diag(np.logspace(0, -8, n)) @ Q2
    k = 3
    e = np.zeros(n)
    e[k] = 1.0
    lu = spla.splu(sp.csc_matrix(A))
    xh, xl, hx = XS.refine_to(DD.DDResidual(sp.csr_matrix(A)), e, lu.solve, lu.solve(e))
    yh, yl, hy = XS.refine_to(DD.DDResidual(sp.csr_matrix(A.T)), e, lambda r: lu.solve(r, trans='T'), lu.solve(e, trans='T'))
    assert hx[-1] <= XS.OMEGA_TOL and hy[-1] <= XS.OMEGA_TOL
    x_ex, y_ex = _exact_solve(A, e), _exact_solve(A.T, e)
    for v_h, v_l, ex in ((xh, xl, x_ex), (yh, yl, y_ex)):
        err = max(abs(float(F(v_h[i]) + F(v_l[i]) - ex[i])) / max(abs(float(ex[i])), 1e-300) for i in range(n))
        assert err < 1e-12
    assert abs((xh[k] + xl[k]) / (yh[k] + yl[k]) - 1) < 1e-14
