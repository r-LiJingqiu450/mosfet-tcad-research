"""Unit tests for the C-4 numerical helpers: the equilibrated solver used by T1 v2 (L1v2) and the
checkpoint test of the bounded job runner."""
import json

import numpy as np
import pytest
import scipy.sparse as sp

from phase3.c4.jobrunner import complete


def _badly_scaled(n=60, seed=1):
    rng = np.random.default_rng(seed)
    A = sp.random(n, n, density=0.15, random_state=rng, format='csr') + sp.eye(n) * 3.0
    rows = sp.diags(10.0 ** rng.uniform(-12, 12, n))          # 24 decades of row scaling
    cols = sp.diags(10.0 ** rng.uniform(-6, 6, n))
    return (rows @ A @ cols).tocsc()


def test_equilibrated_solver_consistent_transformations():
    pytest.importorskip('devsim')
    from phase3.c4.t1v2 import EqSolver
    A = _badly_scaled()
    x_true = np.linspace(1.0, 2.0, A.shape[0])
    b = A @ x_true
    s = EqSolver(A)
    x = s.solve(b)
    # backward stability in the equilibrated metric (C-1(b) B3 metric)
    assert s.eta(A, x, b) < 1e-14
    # forward error of the scaled unknowns y = x / D_c is consistent with kappa(A_s) u
    y, y_true = x / s.dc, x_true / s.dc
    kappa = np.linalg.cond(s.As.toarray(), np.inf)
    fwd = np.max(np.abs(y - y_true)) / np.max(np.abs(y_true))
    assert fwd <= 10 * kappa * 2.0 ** -53
    c = np.zeros(A.shape[0]); c[3] = -1.0
    lam = s.solve_T(c)
    assert np.max(np.abs(A.T @ lam - c)) <= 1e-8 * np.max(np.abs(A.T) @ np.abs(lam))
    # dual-weighted residual is exact in exact arithmetic for a linear output:
    #   c^T(x_hat - x) = lam^T (A x_hat - b);
    # in floating point the residual carries rounding of size ~u(|A||x_hat| + |b|), so the identity
    # holds to within the worst-case rounding term E_cond (the conservative T1 v2 term)
    x_hat = x_true * (1 + 1e-6 * np.sin(np.arange(A.shape[0])))
    true_err = c @ (x_hat - x_true)
    dwr = lam @ (A @ x_hat - b)
    e_cond = 2.0 ** -53 * np.sum(np.abs(lam) * (abs(A) @ np.abs(x_hat) + np.abs(b)))
    assert abs(true_err - dwr) <= 10 * e_cond
    assert abs(true_err - dwr) < 1e-3 * abs(true_err)


def test_jobrunner_checkpoint(tmp_path):
    p = tmp_path / 'r.json'
    assert not complete(str(p))
    p.write_text('{"points": {}}')
    assert not complete(str(p))                                 # incomplete record
    p.write_text(json.dumps({'end_utc': '2026-10-31T00:00:00Z'}))
    assert complete(str(p))
    p.write_text('{"end_utc": ')                                # truncated file
    assert not complete(str(p))
