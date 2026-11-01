"""C-11 independent check: a separate 2-D finite-volume Poisson solver (scipy, no DEVSIM, no shared code with the
decks) for the deck-W geometry at equilibrium, with the zero-bias conductance computed from the resulting electron
density by an independently written network and by the Rayleigh slice/rows bounds.

Geometry (as deck W, half device by symmetry):
- Si 0 <= x <= L/2, 0 <= y <= Y = 1.2 um. Oxide -t_ox <= y <= 0 above the channel.
- Junction face x = 0, 0 <= y <= h_w: Dirichlet psi_bi = V_T ln(N_D/n_i). Below h_w: Neumann.
- Gate y = -t_ox: psi_g = V_G + psi_b0 (V_FB = 0 convention). Bottom: psi_b0. Symmetry plane x = L/2: Neumann.
Physics: Boltzmann, p = n_i e^{-psi/V_T}, n = n_i e^{psi/V_T}, N_A uniform, holes and electrons in equilibrium (F0).
Discretisation: tensor-product mesh, two-point fluxes, control volumes. The conductance network uses Scharfetter-Gummel
equilibrium edge weights c_ij = q mu (h_perp/h_par) n_i e^{psi_m/V_T} with the exact SG average
n_ij = (n_i - n_j)/ln(n_i/n_j). G is taken from the near-zero (source) side.
Usage: python phase3/c11/fd2d_indep.py --na 1e18 --tox 2 --vg 0.783500 --hw 0.5 --L 1 --out <json>
"""
from __future__ import annotations

import argparse
import json
import math

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

Q = 1.602176634e-19
VT = 0.025851999786435316
EPS0 = 8.8541878128e-14
ESI, EOX = 11.7 * EPS0, 3.9 * EPS0
NI, ND = 1e10, 1e20


def graded(a, h0, hmax, total):
    """Nodes from 0 to total: geometric growth (ratio 1.08) from h0 up to hmax, then uniform."""
    xs = [0.0]; h = h0
    while xs[-1] < total - 1e-15:
        xs.append(min(total, xs[-1] + h)); h = min(hmax, h * a)
    return np.array(xs)


def solve(na, tox_nm, vg, hw_um, L_um, ny_ox=8):
    psib = -VT * math.log(na / NI); psibi = VT * math.log(ND / NI); psig = vg + psib
    X, Y, tox = L_um * 0.5e-4, 1.2e-4, tox_nm * 1e-7
    x = graded(1.08, 0.25e-7, 8e-7, X)
    ysi = graded(1.06, 0.2e-7, 8e-7, Y)
    yox = np.linspace(-tox, 0.0, ny_ox + 1)
    y = np.concatenate([yox[:-1], ysi])
    nx, ny = len(x), len(y); N = nx * ny
    idx = lambda i, j: i * ny + j
    jsi0 = ny_ox                                   # index of y = 0 (interface)
    hx, hy = np.diff(x), np.diff(y)
    # control-volume widths
    wx = np.zeros(nx); wx[:-1] += hx / 2; wx[1:] += hx / 2
    wy = np.zeros(ny); wy[:-1] += hy / 2; wy[1:] += hy / 2
    wy_si = np.zeros(ny); wy_si[jsi0:-1] += hy[jsi0:] / 2; wy_si[jsi0 + 1:] += hy[jsi0:] / 2
    eps_y = np.where(np.arange(ny - 1) < jsi0, EOX, ESI)              # vertical edges (between j, j+1)
    # lateral edge permittivity: weighted by the portion of the CV height in oxide / silicon
    wy_ox = wy - wy_si
    eps_x_node = (EOX * wy_ox + ESI * wy_si)                          # eps * height of CV row j
    # Dirichlet sets
    dir_val = np.full(N, np.nan)
    for j in range(ny):
        if j >= jsi0 and y[j] <= hw_um * 1e-4 * (1 + 1e-12):
            dir_val[idx(0, j)] = psibi
    for i in range(nx):
        dir_val[idx(i, 0)] = psig
        dir_val[idx(i, ny - 1)] = psib
    isdir = ~np.isnan(dir_val)
    # edge list with conductance-like coefficients (eps * area / length)
    I, J, C = [], [], []
    for i in range(nx - 1):
        for j in range(ny):
            I.append(idx(i, j)); J.append(idx(i + 1, j)); C.append(eps_x_node[j] / hx[i])
    for i in range(nx):
        for j in range(ny - 1):
            I.append(idx(i, j)); J.append(idx(i, j + 1)); C.append(eps_y[j] * wx[i] / hy[j])
    I, J, C = map(np.array, (I, J, C))
    K = sp.coo_matrix((np.concatenate([C, C]), (np.concatenate([I, J]), np.concatenate([J, I]))), shape=(N, N)).tocsr()
    Lap = (sp.diags(np.asarray(K.sum(axis=1)).ravel()) - K).tocsr()      # div(eps grad) -> -Lap psi
    vol_si = np.outer(wx, wy_si).ravel()
    # initial guess: 1-D far field from a 1-D Newton on the centre column, then a lateral ramp near the face
    psi = np.where(np.isnan(dir_val), psib, dir_val)
    for j in range(ny):
        psi[idx(np.arange(nx), j)] = psig if j == 0 else psib
    free = ~isdir
    def residual(p):
        rho = Q * (NI * np.exp(-p / VT) - NI * np.exp(p / VT) - na)
        return -(Lap @ p) + rho * vol_si
    def jac(p):
        drho = -Q * NI * (np.exp(-p / VT) + np.exp(p / VT)) / VT
        return (-Lap + sp.diags(drho * vol_si)).tocsr()
    # homotopy on the junction-face value (as in the decks): from the column value to psi_bi
    face = np.where(isdir & (np.arange(N) // ny == 0) & ((np.arange(N) % ny) >= jsi0))[0]
    base = None
    hist = []
    for t in np.concatenate([[0.0], np.linspace(0.05, 1.0, 20)]):
        if base is None:
            # t = 0: replace the face Dirichlet by Neumann-like (column value): solve with face values = column values
            dv = dir_val.copy(); dv[face] = np.nan
            isd = ~np.isnan(dv); fr = ~isd; p = psi.copy(); p[isd] = dv[isd]
            for it in range(200):
                F = residual(p); A = jac(p)
                d = spla.spsolve(A[fr][:, fr].tocsc(), -F[fr])
                d = np.clip(d, -2 * VT, 2 * VT); p[fr] += d
                if np.max(np.abs(d)) < 1e-13:
                    break
            base = p.copy(); psi = p.copy(); hist.append(dict(t=0.0, it=it)); continue
        dv = dir_val.copy(); dv[face] = (1 - t) * base[face] + t * psibi
        isd = ~np.isnan(dv); fr = ~isd; psi[isd] = dv[isd]
        for it in range(200):
            F = residual(psi); A = jac(psi)
            d = spla.spsolve(A[fr][:, fr].tocsc(), -F[fr])
            d = np.clip(d, -2 * VT, 2 * VT); psi[fr] += d
            if np.max(np.abs(d)) < 1e-13:
                break
        hist.append(dict(t=float(t), it=it, last=float(np.max(np.abs(d)))))
    P = psi.reshape(nx, ny)
    n = NI * np.exp(P / VT)                                   # electrons (only Si rows used)
    nsi = n[:, jsi0:]; ysi_ = y[jsi0:]; wys = wy_si[jsi0:]
    # slice and rows bounds (continuum-style integrals on the mesh), relative to the centre column
    sigma = nsi @ wys                                          # sheet density per column (mu q cancels)
    sig_inf = sigma[-1]
    integ = 1 - sig_inf / sigma
    dL_slice = float(np.sum(0.5 * (integ[1:] + integ[:-1]) * hx))
    # rows: per row j, R_j = int dx / n ; G_rows = sum_j wy_j / R_j ; compare with sig_inf/(X - dL)
    Rj = np.sum(0.5 * (1 / nsi[1:, :] + 1 / nsi[:-1, :]) * hx[:, None], axis=0)
    touch = ysi_ <= hw_um * 1e-4 * (1 + 1e-12)                 # only rows that touch the source face carry current
    G_rows = float(np.sum((wys / Rj)[touch]))
    dL_rows = X - sig_inf / G_rows
    # network (SG equilibrium weights), Si nodes only; source = face nodes (x = 0, y <= h_w), drain-side: symmetry plane
    # half device: impose phi = 0 on the face and phi = 1/2 on the symmetry column (antisymmetric response), G_half
    ns_ = nsi.ravel(); nxs, nys = nsi.shape
    sid = lambda i, j: i * nys + j
    def sgavg(a, b):
        with np.errstate(divide='ignore', invalid='ignore'):
            r = (a - b) / np.log(a / b)
        return np.where(np.isclose(a, b, rtol=1e-12), 0.5 * (a + b), r)
    eI, eJ, eC = [], [], []
    for i in range(nxs - 1):
        a, b = nsi[i, :], nsi[i + 1, :]
        eI += [sid(i, j) for j in range(nys)]; eJ += [sid(i + 1, j) for j in range(nys)]; eC += list(sgavg(a, b) * wys / hx[i])
    for i in range(nxs):
        a, b = nsi[i, :-1], nsi[i, 1:]
        eI += [sid(i, j) for j in range(nys - 1)]; eJ += [sid(i, j + 1) for j in range(nys - 1)]; eC += list(sgavg(a, b) * wx[i] / np.diff(ysi_))
    eI, eJ, eC = map(np.array, (eI, eJ, eC))
    W = sp.coo_matrix((np.concatenate([eC, eC]), (np.concatenate([eI, eJ]), np.concatenate([eJ, eI]))), shape=(len(ns_),) * 2).tocsr()
    Lw = (sp.diags(np.asarray(W.sum(axis=1)).ravel()) - W).tocsr()
    src = np.array([sid(0, j) for j in range(nys) if ysi_[j] <= hw_um * 1e-4 * (1 + 1e-12)])
    drn = np.array([sid(nxs - 1, j) for j in range(nys)])
    fixed = np.zeros(len(ns_), bool); fixed[src] = True; fixed[drn] = True
    phi = np.zeros(len(ns_)); phi[drn] = 1.0
    fr = np.where(~fixed)[0]
    A = Lw[fr][:, fr].tocsc(); bvec = -Lw[fr][:, np.where(fixed)[0]] @ phi[fixed]
    s = 1 / np.sqrt(A.diagonal()); y_ = spla.splu((sp.diags(s) @ A @ sp.diags(s)).tocsc()).solve(s * bvec)
    phi[fr] = s * y_
    G_src = -float(np.sum((Lw @ phi)[src]))                   # source-side current (near-zero side)
    dL_net = X - sig_inf / G_src                              # half device: G_half = sig_inf/(X - dL)
    surf = P[:, jsi0]
    return dict(args=dict(na=na, tox_nm=tox_nm, vg=vg, hw_um=hw_um, L_um=L_um), nx=nx, ny=ny,
                dL_net_nm=dL_net * 1e7, dL_slice_nm=dL_slice * 1e7, dL_rows_nm=dL_rows * 1e7,
                psi_s_centre=float(surf[-1]), psi_s_face=float(surf[0]), homotopy=hist)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--na', type=float, required=True); ap.add_argument('--tox', type=float, required=True)
    ap.add_argument('--vg', type=float, required=True); ap.add_argument('--hw', type=float, default=0.5)
    ap.add_argument('--L', type=float, default=1.0); ap.add_argument('--out', required=True)
    a = ap.parse_args()
    r = solve(a.na, a.tox, a.vg, a.hw, a.L)
    json.dump(r, open(a.out, 'w'), indent=1)
    print(f"nx {r['nx']} ny {r['ny']} | dL net {r['dL_net_nm']:.3f} slice {r['dL_slice_nm']:.3f} rows {r['dL_rows_nm']:.3f} nm | psi_s centre {r['psi_s_centre']:.6f}")


if __name__ == '__main__':
    main()
