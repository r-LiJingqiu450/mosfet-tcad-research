"""Independent 1-D FV equilibrium solver on the DEVSIM column y-lines (1COL data).

Purpose (C-1 deck repair, conditions D1/D2 of R04 v2.1 5.2): provide the lateral
psi Dirichlet data psi_1(y; V_G, 0) on exactly the 2-D x-lines (here DEVSIM's y-lines),
computed by the same FV discretisation a laterally uniform 2-D solution satisfies:

  oxide node i : eps_ox * sum_j (psi_j - psi_i)/h_ij = 0
  Si node i    : eps_si * sum_j (psi_j - psi_i)/h_ij + q (p_i - n_i + N) w_i = 0
  interface    : oxide and Si flux terms summed, Si charge with the Si half-cell w_0
  gate (top)   : psi = V_G + psi_b0 ;  bottom : psi = psi_b0
  n = n_i exp(psi/V_T), p = n_i exp(-psi/V_T) (phi_n = phi_p = 0, frozen holes).
w_i are the 1-D dual widths (trapezoid weights, half cells at ends) = D4.
Units cm. Newton with an exact tridiagonal Jacobian; stops when max|dpsi| <= tol.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.linalg import solve_banded

from phase3.c1 import deck as D


def solve_1col(y_ox, y_si, vg, tol=1e-14, max_it=200, psi0=None, phin=0.0):
    """y_ox: oxide y-lines ascending ending at 0; y_si: Si y-lines from 0 ascending.
    phin: spatially uniform electron quasi-Fermi potential (V); n = n_i exp((psi-phin)/V_T).
    phin = 0 reproduces the original arithmetic exactly (psi - 0.0 == psi)."""
    y_ox = np.asarray(y_ox, float); y_si = np.asarray(y_si, float)
    assert abs(y_ox[-1]) < 1e-20 and abs(y_si[0]) < 1e-20
    y = np.concatenate([y_ox[:-1], y_si])          # shared interface node
    k = len(y_ox) - 1                                # index of interface node
    N = len(y); h = np.diff(y)
    eps = np.where(np.arange(N - 1) < k, 3.9, 11.7) * D.EPS0_CM   # per segment
    w = np.zeros(N)                                  # Si dual widths (charge)
    hs = np.diff(y_si)
    ws = np.zeros(len(y_si)); ws[:-1] += hs / 2; ws[1:] += hs / 2
    w[k:] = ws
    VT, ni, q, Nd = D.VT, D.NI, D.Q, -D.NA
    psi = np.full(N, D.PSI_B0) if psi0 is None else np.array(psi0, float)
    psi[0] = vg + D.PSI_B0; psi[-1] = D.PSI_B0
    hist = []
    for it in range(max_it):
        n = ni * np.exp((psi - phin) / VT); p = ni * np.exp(-psi / VT)
        F = np.zeros(N); lo = np.zeros(N); di = np.zeros(N); up = np.zeros(N)
        g = eps / h                                       # conductances per segment
        # flux terms
        F[:-1] += g * (psi[1:] - psi[:-1]); F[1:] += g * (psi[:-1] - psi[1:])
        di[:-1] -= g; di[1:] -= g
        up[:-1] += g; lo[1:] += g
        # charge (Si nodes only)
        si = np.arange(N) >= k
        F[si] += q * (p[si] - n[si] + Nd) * w[si]
        di[si] += q * (-p[si] / VT - n[si] / VT) * w[si]
        # Dirichlet rows
        for i, val in ((0, vg + D.PSI_B0), (N - 1, D.PSI_B0)):
            F[i] = psi[i] - val; di[i] = 1.; lo[i] = 0.; up[i] = 0.
        ab = np.zeros((3, N)); ab[0, 1:] = up[:-1]; ab[1] = di; ab[2, :-1] = lo[1:]
        # fix Dirichlet coupling: rows 0 and N-1 have no off-diagonals
        ab[0, 1] = 0.; ab[2, N - 2] = 0.
        dpsi = solve_banded((1, 1), ab, -F)
        # damping for robustness far from the solution (|step| <= 1 V_T); undamped near it
        m = np.max(np.abs(dpsi))
        if m > VT:
            dpsi *= VT / m
        psi += dpsi
        hist.append(float(np.max(np.abs(dpsi))))
        if hist[-1] <= tol and m <= VT:
            break
    return dict(y=y, psi=psi, k=k, converged=hist[-1] <= tol, iterations=len(hist),
                last_update_V=hist[-1], history=hist[-6:],
                n=ni * np.exp((psi - phin) / VT), w_si=w, phin=phin)
