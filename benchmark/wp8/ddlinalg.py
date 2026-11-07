"""Double-double residuals and mixed-precision iterative refinement for sparse linear systems (WP8).

A double-double number is the unevaluated sum hi + lo of two doubles, |lo| <= ulp(hi)/2, worth about 32 significant
digits. The error-free transformations below (TwoSum, Veltkamp splitting, Dekker's TwoProduct) are the classical ones
(Higham, Accuracy and Stability of Numerical Algorithms, 2nd ed., SIAM 2002, sections 4.3 and 4.4).

Only the residual r = b - A x and the accumulation of x are carried in double-double. Each correction is solved in
double precision with a factorization the caller supplies; this is mixed-precision iterative refinement. It converges
when that factorization is accurate enough for the scaled matrix, and the attainable accuracy is then set by the
double-double residual, not by the factorization.

All functions work elementwise on NumPy float64 arrays; nothing here depends on DEVSIM.
"""
from __future__ import annotations

import numpy as np

SPLITTER = 134217729.0          # 2**27 + 1 (Veltkamp splitting constant for binary64)


def two_sum(a, b):
    """s + e == a + b exactly, with s = fl(a + b)."""
    s = a + b
    bb = s - a
    e = (a - (s - bb)) + (b - bb)
    return s, e


def fast_two_sum(a, b):
    """s + e == a + b exactly, with s = fl(a + b), provided |a| >= |b| (or a == 0)."""
    s = a + b
    e = b - (s - a)
    return s, e


def split(a):
    """a == hi + lo exactly, each with at most 26 significant bits."""
    c = SPLITTER * a
    hi = c - (c - a)
    return hi, a - hi


def two_prod(a, b):
    """p + e == a * b exactly, with p = fl(a * b) (no overflow or underflow assumed)."""
    p = a * b
    ah, al = split(a)
    bh, bl = split(b)
    e = ((ah * bh - p) + ah * bl + al * bh) + al * bl
    return p, e


class DDResidual:
    """r = b - A x in double-double for a fixed sparse A (CSR), x given as (xh, xl)."""

    def __init__(self, A):
        A = A.tocsr()
        A.sort_indices()
        self.A = A
        self.n = A.shape[0]
        nnz = np.diff(A.indptr)
        self.order = np.argsort(-nnz, kind='stable')            # rows by decreasing number of entries
        self.counts = nnz[self.order]
        self.neg_counts = -self.counts                         # ascending, for searchsorted
        self.start = A.indptr[:-1][self.order]
        self.kmax = int(nnz.max()) if self.n else 0
        self.absA = abs(A)

    def __call__(self, b, xh, xl):
        S = np.asarray(b, dtype=np.float64)[self.order].copy()
        E = np.zeros(self.n)
        data, idx = self.A.data, self.A.indices
        for k in range(self.kmax):
            m = int(np.searchsorted(self.neg_counts, -k, side='left'))   # rows with more than k entries: a prefix
            if m == 0:
                break
            pos = self.start[:m] + k
            a = data[pos]
            j = idx[pos]
            p, e = two_prod(a, xh[j])
            q = a * xl[j]
            s, t = two_sum(S[:m], -p)
            S[:m] = s
            E[:m] += t - e - q
        rh, rl = two_sum(S, E)
        out_h = np.empty(self.n)
        out_l = np.empty(self.n)
        out_h[self.order] = rh
        out_l[self.order] = rl
        return out_h, out_l

    def componentwise_backward_error(self, r, b, x):
        """max_i |r_i| / (|A| |x| + |b|)_i, with the denominator in double (it is a scale)."""
        den = self.absA @ np.abs(x) + np.abs(b)
        with np.errstate(divide='ignore', invalid='ignore'):
            w = np.where(den > 0, np.abs(r) / den, np.where(np.abs(r) > 0, np.inf, 0.0))
        return float(np.max(w))


def refine(A, b, solve, x0, out_index, sign=-1.0, max_it=40, rtol=1e-14, n_stable=2, residual=None):
    """Mixed-precision iterative refinement of A x = b from x0 with residuals in double-double.

    solve(r) returns an approximate solution of A d = r in double. The output is G = sign * x[out_index]. Converged
    when |G_k - G_{k-1}| <= rtol |G_k| holds for n_stable consecutive iterations. Returns a dict with the final G (as
    hi, lo), the history, and whether it converged.
    """
    res = residual or DDResidual(A)
    b = np.asarray(b, dtype=np.float64)
    xh = np.array(x0, dtype=np.float64)
    xl = np.zeros_like(xh)
    hist, G_prev, stable, converged = [], None, 0, False
    for k in range(max_it + 1):
        rh, rl = res(b, xh, xl)
        r = rh + rl
        Gh, Gl = sign * xh[out_index], sign * xl[out_index]
        G = float(Gh + Gl)
        omega = res.componentwise_backward_error(r, b, xh)
        hist.append(dict(it=k, G=G, G_lo=float(Gl), omega_dd=omega, res_inf=float(np.max(np.abs(r)))))
        if G_prev is not None:
            stable = stable + 1 if abs(G - G_prev) <= rtol * abs(G) else 0
            if stable >= n_stable:
                converged = True
                break
        if k == max_it:
            break
        d = solve(r)
        s, t = two_sum(xh, d)
        xh, xl = fast_two_sum(s, t + xl)
        G_prev = G
    return dict(G=hist[-1]['G'], G_lo=hist[-1]['G_lo'], converged=converged, iterations=len(hist) - 1, history=hist,
                x_hi=xh, x_lo=xl)
