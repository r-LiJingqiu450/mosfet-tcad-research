"""C-11: parameter-free prediction of the far-field decay length lam and the connection coefficient B of the zero-bias
subthreshold end correction (deck W geometry, F0 physics), from 1-D quantities only (no 2-D solve).

Model H-lin (linearised eigen-projection):
 1. Far field: exact 1-D Boltzmann profile b(y) = (psi_inf(y) - psi_b0)/V_T from the first integral
    (db/dy)^2 = (2/L_D^2)[(e^{-b} + b - 1) + (n0/p0)(e^{b} - b - 1)], with b(0) = b_s from the 1-D gate equation.
 2. Linearised end perturbation u = psi - psi_inf obeys div(eps grad u) = q (p_inf + n_inf) u / V_T in Si and
    div(eps grad u) = 0 in the oxide.  Separable modes u = Phi(y) e^{-k x}:
        -(eps Phi')' + q (p_inf + n_inf)/V_T Phi = k^2 eps Phi,   Phi(-t_ox) = 0 (gate),  Phi(Y) = 0 (deep bulk).
    lam = 1/k0 for the lowest mode; Phi0 normalised to Phi0(0) = 1 (Si/oxide interface).
 3. Junction face x = 0 carries u(0, y) = U(y):  U = psi_bi - psi_inf(y) on the Dirichlet face 0 <= y <= h_w;
    U = 0 below h_w (deck W has Neumann there; irrelevant for h_w >> W_d); linear in the oxide from U(0) to 0 at the
    gate.  B = <U, Phi0>_eps / <Phi0, Phi0>_eps (Sturm-Liouville projection with weight eps).
 4. dL_pred = lam * Ein(B/V_T)   (C-10 structure; post hoc in C-10, challenged here with predicted lam and B).
Known omission (stated before testing): the linearisation cannot represent the nonlinear lateral depletion bulge of the
junction below the channel depletion edge, nor electron screening at the face (n_face/N_A); both are expected to make the
true B larger than B_lin if they matter.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eigh_tridiagonal
from scipy.special import exp1

Q = 1.602176634e-19
VT = 0.025851999786435316
EPS0 = 8.8541878128e-14
ESI, EOX = 11.7 * EPS0, 3.9 * EPS0
GAMMA = 0.5772156649015329


def ein(z):
    return GAMMA + math.log(z) + float(exp1(z))


def far_profile(na, b_s, ni=1e10, n_ld=40.0, dy_ld=0.01):
    """1-D far-field b(y) (y in units of L_D) on a uniform grid from 0 to Y = W/L_D + n_ld."""
    p0 = 0.5 * (na + math.sqrt(na * na + 4 * ni * ni)); n0 = ni * ni / p0
    r = n0 / p0
    F = lambda b: (math.exp(-b) + b - 1) + r * (math.expm1(b) - b)
    wd = math.sqrt(2 * max(b_s - 1, 1e-6))
    Y = wd + n_ld
    ys = np.arange(0.0, Y + dy_ld / 2, dy_ld)
    ys = ys[ys <= Y]
    sol = solve_ivp(lambda y, b: [-math.sqrt(max(2 * F(max(b[0], 0.0)), 0.0))], (0.0, Y), [b_s], t_eval=ys,
                    rtol=1e-11, atol=1e-13, max_step=0.05)
    b = np.maximum(sol.y[0], 0.0)
    return ys, b, p0, n0


def predict(na, tox_nm, b_s, psibi_hat, hw_nm, ni=1e10, n_ox=40, dy_ld=0.01):
    """Return dict(lam_nm, B_V, dL_pred_nm, ...); b_s and psibi_hat in units of V_T relative to the bulk."""
    LD = math.sqrt(ESI * VT / (Q * na))                       # cm
    ys, b, p0, n0 = far_profile(na, b_s, ni, dy_ld=dy_ld)
    tox = tox_nm * 1e-7
    # grid: oxide (uniform, n_ox cells) + silicon (uniform dy_ld L_D); unknowns at interior nodes
    yo = np.linspace(-tox, 0.0, n_ox + 1)                    # cm
    ysi = ys * LD
    y = np.concatenate([yo, ysi[1:]])
    eps = np.where(y[:-1] + np.diff(y) / 2 < 0, EOX, ESI)    # edge permittivities
    h = np.diff(y)
    # node reaction term q (p + n)/V_T (Si only) and node weights (control volumes)
    bb = np.concatenate([np.full(n_ox, np.nan), b])          # b at Si nodes incl. interface (index n_ox)
    react = np.zeros(len(y))
    si = np.arange(n_ox, len(y))
    react[si] = Q * (p0 * np.exp(-bb[si]) + n0 * np.exp(bb[si])) / VT
    vol = np.zeros(len(y)); vol[:-1] += h / 2; vol[1:] += h / 2
    epsn = np.zeros(len(y)); epsn[:-1] += eps * h / 2; epsn[1:] += eps * h / 2      # int eps over control volume
    # interior nodes 1..N-2 (Dirichlet at both ends); symmetric generalized problem A phi = k^2 M phi, M = diag(epsn)
    c = eps / h
    idx = np.arange(1, len(y) - 1)
    diag = c[idx - 1] + c[idx] + react[idx] * vol[idx]
    off = -c[idx[:-1]]
    m = epsn[idx]
    s = 1 / np.sqrt(m)
    w, v = eigh_tridiagonal(diag * s * s, off * s[:-1] * s[1:], select='i', select_range=(0, 0))
    k2 = w[0]; phi = np.zeros(len(y)); phi[idx] = v[:, 0] * s
    phi /= phi[n_ox]                                          # Phi0(0) = 1 at the interface
    lam = 1 / math.sqrt(k2)
    # boundary data U(y) on x = 0
    psi_inf = np.zeros(len(y)); psi_inf[si] = bb[si] * VT    # relative to the bulk
    U = np.zeros(len(y))
    face = si[(y[si] <= hw_nm * 1e-7 * (1 + 1e-12))]
    U[face] = psibi_hat * VT - psi_inf[face]
    U[:n_ox + 1] = (psibi_hat * VT - b_s * VT) * (y[:n_ox + 1] + tox) / tox
    B = float(np.sum(epsn * U * phi) / np.sum(epsn * phi * phi))
    return dict(lam_nm=lam * 1e7, B_V=B, B_over_VT=B / VT, dL_pred_nm=lam * ein(B / VT) * 1e7 if B > 0 else float('nan'),
                Wd_ld=math.sqrt(2 * max(b_s - 1, 1e-6)), LD_nm=LD * 1e7)
