#!/usr/bin/env python3
"""WP2 path P4: independent finite-volume reference for the resolved-well deck R (pre-registration section 5).

Independence rules: imports numpy/scipy only (nothing from phase3/, semiconductor_sim/, devsim or the other WP2
scripts); constants re-declared below; own tensor mesh, own box-integration assembly, own damped Newton solve of
the equilibrium Poisson-Boltzmann problem, own Scharfetter-Gummel zero-bias network and terminal solves.

Device (pre-registration section 4): Si x in [-l_w, L + l_w], y in [0, 1.2 um]; SiO2 y in [-t_ox, 0] over the full
width; gate (psi = V_G + psi_b0) on the oxide top for 0 <= x <= L; ohmic source/drain on the far Si edges for
0 <= y <= h_w (psi = V_T asinh(N_net/2n_i)); body y = 1.2 um psi = psi_b0 (all x); other boundaries homogeneous
Neumann. Donors N_D = 1e20 f_x(x) f_y(y) (erfc laterally and vertically, sigma = 15 nm, metallurgical junction at
x = 0 and L, f_y centred at h_w = 0.5 um). Boltzmann electrons and holes at equilibrium (phi_n = phi_p = 0).

Zero-bias network: Si edge conductance c = q mu_n (a_Si / h) n_i exp(psi_a/V_T) B((psi_a - psi_b)/V_T); terminals
source, drain (+ body for B1). Outputs G_DD (edge energy and -(I_S + I_B)), G_SD (-I_S and reciprocity), G_DB.

Mesh rules (pre-registration): spacing <= 0.5 nm at the interface and at x = 0, L; <= 10 nm everywhere.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import os
import platform
import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.special import erfc, erfcinv

VERSION = 'wp2-fv-rdeck-2026-11-04'
# ---- constants (re-declared; CODATA 2018 exact SI values)
Q = 1.602176634e-19
KB = 1.380649e-23
TEMP = 300.0
VT = KB * TEMP / Q
EPS0 = 8.8541878128e-14            # F/cm
EPS_SI, EPS_OX = 11.7 * EPS0, 3.9 * EPS0
NI = 1.0e10
MU_N = 1350.0
ND_PEAK = 1.0e20
HW, LW, SIG, SIGY, YMAX = 0.5e-4, 0.3e-4, 15e-7, 15e-7, 1.2e-4   # cm


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec='milliseconds')


# ------------------------------------------------------------------------------------------------- mesh
def march(a, b, hfun):
    """Nodes a = x_0 < ... < x_N = b with local spacing ~ hfun(x), rescaled to end exactly at b."""
    xs = [a]
    while xs[-1] < b:
        xs.append(xs[-1] + hfun(xs[-1]))
    h = np.diff(xs)
    h = h * (b - a) / h.sum()                         # sum(h) >= b - a: cells only shrink (caps are kept)
    return a + np.concatenate([[0.0], np.cumsum(h)])


def axis(hard, hfun):
    out = [np.array([hard[0]])]
    for a, b in zip(hard[:-1], hard[1:]):
        out.append(march(a, b, hfun)[1:])
    x = np.concatenate(out)
    x[-1] = hard[-1]
    return x


def make_mesh(L, tox, h0x=0.25e-7, h0y=0.05e-7, gx=1.08, gy=1.05, hmax=10e-7, hj=2.5e-7, y_hard=()):
    """Graded tensor mesh. Lateral: fine (h0x) at x = 0 and L, geometric growth gx, cap hmax.
    Depth: fine (h0y) at y = 0, growth gy, cap hmax; cap hj in the well-bottom junction band
    [h_w - 60 nm, h_w + 120 nm]. Oxide: uniform, >= 8 cells, cell <= 0.25 nm or t_ox/8."""
    def hx(x):
        d = min(abs(x), abs(x - L))
        return min(hmax, h0x + (gx - 1.0) * d)

    def hy(y):
        h = min(hmax, h0y + (gy - 1.0) * y)
        if HW - 60e-7 <= y <= HW + 120e-7:
            h = min(h, hj)
        return h
    xs = axis([-LW, 0.0, L, L + LW], hx)
    ysi = axis(sorted({0.0, HW, YMAX, *[float(v) for v in y_hard]}), hy)
    nox = max(8, int(math.ceil(tox / 0.25e-7))) if tox <= 2e-7 else 8
    yox = np.linspace(-tox, 0.0, nox + 1)
    return xs, yox, ysi


# ------------------------------------------------------------------------------------------------- device
def psi_b0(na):
    p0 = 0.5 * (na + math.sqrt(na * na + 4 * NI * NI))
    return -VT * math.log(p0 / NI)


def donors(x, y, L, na):
    s2 = math.sqrt(2.0) * SIG
    xc = -s2 * float(erfcinv(2.0 * na / ND_PEAK))
    fx = 0.5 * erfc((x - xc) / s2) + 0.5 * erfc(((L - x) - xc) / s2)
    fy = 0.5 * erfc((y - HW) / (math.sqrt(2.0) * SIGY))
    return ND_PEAK * fx * fy


class Grid:
    def __init__(self, xs, yox, ysi):
        self.x = xs
        self.y = np.concatenate([yox[:-1], ysi])
        self.kif = len(yox) - 1
        self.nx, self.ny = len(xs), len(self.y)
        self.N = self.nx * self.ny
        X, Y = np.meshgrid(self.x, self.y, indexing='ij')         # node (i, j) -> i * ny + j
        self.X, self.Y = X.ravel(), Y.ravel()
        hx, hy = np.diff(self.x), np.diff(self.y)
        si_band = self.y[1:] > 0.0                                 # band [y_j, y_j+1] is Si
        idx = np.arange(self.N).reshape(self.nx, self.ny)
        # half-widths of dual cells
        wx = np.zeros(self.nx); wx[:-1] += hx / 2; wx[1:] += hx / 2
        up = np.zeros(self.ny); dn = np.zeros(self.ny)               # half heights above / below node j
        up[1:] = hy / 2; dn[:-1] = hy / 2
        up_si = np.zeros(self.ny); dn_si = np.zeros(self.ny)
        up_si[1:] = np.where(si_band, hy / 2, 0.0); dn_si[:-1] = np.where(si_band, hy / 2, 0.0)
        up_eps = np.zeros(self.ny); dn_eps = np.zeros(self.ny)
        up_eps[1:] = np.where(si_band, EPS_SI, EPS_OX) * hy / 2
        dn_eps[:-1] = np.where(si_band, EPS_SI, EPS_OX) * hy / 2
        # horizontal edges (i, j)-(i+1, j)
        a_h = idx[:-1, :].ravel(); b_h = idx[1:, :].ravel()
        len_h = np.repeat(hx, self.ny)
        face_eps_h = np.tile(up_eps + dn_eps, self.nx - 1)
        face_si_h = np.tile(up_si + dn_si, self.nx - 1)
        # vertical edges (i, j)-(i, j+1)
        a_v = idx[:, :-1].ravel(); b_v = idx[:, 1:].ravel()
        len_v = np.tile(hy, self.nx)
        face_v = np.repeat(wx, self.ny - 1)
        band_si_v = np.tile(si_band, self.nx)
        self.ea = np.concatenate([a_h, a_v]); self.eb = np.concatenate([b_h, b_v])
        self.elen = np.concatenate([len_h, len_v])
        self.kpsi = np.concatenate([face_eps_h, np.where(band_si_v, EPS_SI, EPS_OX) * face_v]) / self.elen
        self.ksi = np.concatenate([face_si_h, np.where(band_si_v, face_v, 0.0)]) / self.elen
        self.vol_si = (np.outer(wx, up_si + dn_si)).ravel()          # Si part of the dual cell
        self.is_si = self.Y >= 0.0


def assemble_couplings(N, a, b, k):
    W = sp.coo_matrix((np.concatenate([k, k]), (np.concatenate([a, b]), np.concatenate([b, a]))), shape=(N, N)).tocsr()
    return (sp.diags(np.asarray(W.sum(axis=1)).ravel()) - W).tocsr()


def bern(v):
    v = np.asarray(v, float)
    out = np.empty_like(v)
    s = np.abs(v) < 1e-6
    out[s] = 1 - v[s] / 2 + v[s] ** 2 / 12
    out[~s] = v[~s] / np.expm1(v[~s])
    return out


# ------------------------------------------------------------------------------------------------- equilibrium
def equilibrium(g, L, tox, na, vg, tol=1e-12, max_it=200, contact_extent=HW):
    pb = psi_b0(na)
    netd = np.where(g.is_si, donors(g.X, np.maximum(g.Y, 0.0), L, na) - na, 0.0)
    tol_x = 1e-6 * min(np.diff(g.x).min(), 1e-9)
    left, right = np.abs(g.X - g.x[0]) < tol_x, np.abs(g.X - g.x[-1]) < tol_x
    contact = (left | right) & g.is_si & (g.Y <= contact_extent * (1 + 1e-12))
    gate = (np.abs(g.Y - g.y[0]) < 1e-18) & (g.X >= -tol_x) & (g.X <= L + tol_x)
    body = np.abs(g.Y - g.y[-1]) < 1e-18
    fixed = contact | gate | body
    psi = np.where(g.is_si, VT * np.arcsinh(netd / (2 * NI)), pb)
    psi[gate] = vg + pb
    psi[body] = pb
    psi[contact] = VT * np.arcsinh(netd[contact] / (2 * NI))
    Kp = assemble_couplings(g.N, g.ea, g.eb, g.kpsi)
    fr = np.where(~fixed)[0]
    Kff = Kp[fr][:, fr].tocsc()
    hist = []
    for it in range(max_it):
        n = np.where(g.is_si, NI * np.exp(psi / VT), 0.0)
        p = np.where(g.is_si, NI * np.exp(-psi / VT), 0.0)
        F = -(Kp @ psi) + Q * g.vol_si * (p - n + netd)                 # div(eps grad psi) + rho = 0
        dF = -Q * g.vol_si * (p + n) / VT
        J = (-Kff + sp.diags(dF[fr])).tocsc()
        dx = spla.spsolve(J, -F[fr])
        step = np.sign(dx) * VT * np.log1p(np.abs(dx) / VT)            # logarithmic damping
        psi[fr] += step
        hist.append(dict(it=it, max_update_V=float(np.max(np.abs(step))), max_raw_V=float(np.max(np.abs(dx)))))
        if np.max(np.abs(step)) <= tol:
            break
    for _ in range(2):                                                  # two extra iterations, recorded
        n = np.where(g.is_si, NI * np.exp(psi / VT), 0.0); p = np.where(g.is_si, NI * np.exp(-psi / VT), 0.0)
        F = -(Kp @ psi) + Q * g.vol_si * (p - n + netd)
        dF = -Q * g.vol_si * (p + n) / VT
        dx = spla.spsolve((-Kff + sp.diags(dF[fr])).tocsc(), -F[fr])
        psi[fr] += dx
        hist.append(dict(it='extra', max_update_V=float(np.max(np.abs(dx)))))
    n = np.where(g.is_si, NI * np.exp(psi / VT), 0.0); p = np.where(g.is_si, NI * np.exp(-psi / VT), 0.0)
    F = -(Kp @ psi) + Q * g.vol_si * (p - n + netd)
    scale = np.abs(Kp) @ np.abs(psi) + Q * g.vol_si * (p + n + np.abs(netd))
    res = float(np.max(np.abs(F[fr]) / np.maximum(scale[fr], 1e-300)))
    converged = hist[-3]['max_update_V'] <= tol
    return psi, dict(iterations=len(hist), converged=bool(converged), last=hist[-4:], rel_residual=res,
                     psi_b0=pb, masks=dict(contact=contact, gate=gate, body=body))


# ------------------------------------------------------------------------------------------------- network
def network_couplings(g, psi):
    a, b = g.ea, g.eb
    m = g.ksi > 0
    v = (psi[a[m]] - psi[b[m]]) / VT
    c = Q * MU_N * g.ksi[m] * NI * np.exp(psi[a[m]] / VT) * bern(v)
    return a[m], b[m], c


def terminal_solve(N, a, b, c, terminals, excite):
    Lap = assemble_couplings(N, a, b, c)
    fixed = np.zeros(N, bool)
    for idx in terminals.values():
        fixed[idx] = True
    active = np.zeros(N, bool); active[a] = True; active[b] = True     # Si nodes only
    phi = np.zeros(N); phi[terminals[excite]] = 1.0
    fr = np.where(~fixed & active)[0]
    A = Lap[fr][:, fr].tocsc(); rhs = -Lap[fr][:, np.where(fixed)[0]] @ phi[fixed]
    s = 1.0 / np.sqrt(A.diagonal())
    As = (sp.diags(s) @ A @ sp.diags(s)).tocsc()
    phi[fr] = s * spla.splu(As).solve(s * rhs)
    res = float(np.max(np.abs(A @ phi[fr] - rhs)) / np.max(np.abs(rhs)))
    I = Lap @ phi
    return dict(I={t: float(np.sum(I[idx])) for t, idx in terminals.items()},
                energy=float(np.sum(c * (phi[a] - phi[b]) ** 2)), residual=res)


def conductances(N, a, b, c, terms):
    rd = terminal_solve(N, a, b, c, terms, 'drain')
    rs = terminal_solve(N, a, b, c, terms, 'source')
    I_S, I_B = rd['I']['source'], rd['I'].get('body', 0.0)
    out = dict(G_DD=rd['energy'], G_DD_second=-(I_S + I_B), G_SD=-I_S, G_SD_second=-rs['I']['drain'],
               G_DD_direct_drain_sum=rd['I']['drain'], residual=max(rd['residual'], rs['residual']))
    out['G_DD_rel_forms'] = out['G_DD_second'] / out['G_DD'] - 1
    out['G_SD_rel_forms'] = out['G_SD_second'] / out['G_SD'] - 1
    if 'body' in terms:
        out['G_DB'] = -I_B
        out['kcl_closure_rel'] = abs(out['G_DD'] - out['G_SD'] - out['G_DB']) / out['G_DD']
    return out


def run(L_um, tox_nm, na, vg, contact_extent_um=None, **mesh_kw):
    L, tox = L_um * 1e-4, tox_nm * 1e-7
    xs, yox, ysi = make_mesh(L, tox, **mesh_kw)
    g = Grid(xs, yox, ysi)
    t0 = time.time()
    ce = HW if contact_extent_um is None else contact_extent_um * 1e-4
    psi, eq = equilibrium(g, L, tox, na, vg, contact_extent=ce)
    masks = eq.pop('masks')
    a, b, c = network_couplings(g, psi)
    left = masks['contact'] & (g.X < L / 2); right = masks['contact'] & (g.X > L / 2)
    body = masks['body'] & ~masks['contact']
    res = {}
    for name, terms in (('B0', dict(source=np.where(left)[0], drain=np.where(right)[0])),
                        ('B1', dict(source=np.where(left)[0], drain=np.where(right)[0], body=np.where(body)[0]))):
        res[name] = conductances(g.N, a, b, c, terms)
    # centre-column sheet conductance and surface potential (diagnostics)
    ic = int(np.argmin(np.abs(g.x - L / 2)))
    col = np.arange(ic * g.ny, (ic + 1) * g.ny)
    n = np.where(g.is_si, NI * np.exp(psi / VT), 0.0)
    wsi = np.zeros(g.ny); hy = np.diff(g.y); si_band = g.y[1:] > 0
    wsi[:-1] += np.where(si_band, hy / 2, 0); wsi[1:] += np.where(si_band, hy / 2, 0)
    Ns_centre = float(np.sum(n[col] * wsi))
    psi_s = float(psi[col][g.kif])
    hxs, hys = np.diff(xs), np.diff(ysi)
    mesh = dict(nx=int(g.nx), ny=int(g.ny), nodes=int(g.N), n_ox_cells=int(len(yox) - 1),
                dx_min_nm=float(hxs.min() * 1e7), dx_max_nm=float(hxs.max() * 1e7),
                dx_at_x0_nm=float(min(hxs[np.argmin(np.abs(xs[:-1]))], hxs[np.argmin(np.abs(xs[1:]))]) * 1e7),
                dx_at_L_nm=float(min(hxs[np.argmin(np.abs(xs[:-1] - L))], hxs[np.argmin(np.abs(xs[1:] - L))]) * 1e7),
                dy_si_min_nm=float(hys.min() * 1e7), dy_si_max_nm=float(hys.max() * 1e7), dy_at_interface_nm=float(hys[0] * 1e7),
                h_w_is_mesh_line=bool(np.any(np.abs(ysi - HW) < 1e-15)), contact_extent_um=ce * 1e4,
                kwargs={k: (list(v) if isinstance(v, tuple) else v) for k, v in mesh_kw.items()})
    assert mesh['dx_at_x0_nm'] <= 0.5 and mesh['dx_at_L_nm'] <= 0.5 and mesh['dy_at_interface_nm'] <= 0.5
    assert mesh['dx_max_nm'] <= 10.0 + 1e-9 and mesh['dy_si_max_nm'] <= 10.0 + 1e-9
    return dict(mesh=mesh, equilibrium=eq, results=res, psi_s_centre_V=psi_s,
                psi_hat_centre=(psi_s - eq['psi_b0']) / VT, Ns_centre_cm2=Ns_centre,
                G_centre_sheet_over_L=Q * MU_N * Ns_centre / L, t_s=time.time() - t0)


MESH_PRESETS = {
    'base': {},
    'refined': dict(h0x=0.125e-7, h0y=0.025e-7, gx=1.04, gy=1.025, hmax=5e-7, hj=1.25e-7),   # production (amendment A3)
}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--id', required=True)
    ap.add_argument('--na', type=float, required=True)
    ap.add_argument('--tox', type=float, required=True)
    ap.add_argument('--L', type=float, required=True)
    ap.add_argument('--vg', type=float, required=True)
    ap.add_argument('--mesh', choices=sorted(MESH_PRESETS), default='refined')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rec = dict(version=VERSION, path='P4', start_utc=now(), args=vars(a), argv=sys.argv,
               environment=dict(python=sys.version.split()[0], platform=platform.platform(),
                                numpy=np.__version__, scipy=__import__('scipy').__version__))
    rec.update(run(a.L, a.tox, a.na, a.vg, **MESH_PRESETS[a.mesh]))
    rec['mesh_preset'] = a.mesh
    rec['end_utc'] = now()
    with open(os.path.join(a.out, f'P4_{a.id}_{a.mesh}.json'), 'w', encoding='utf-8') as f:
        json.dump(rec, f, indent=1)
    r = rec['results']
    print(f"{a.id} P4[{a.mesh}] nodes={rec['mesh']['nodes']} conv={rec['equilibrium']['converged']} res={rec['equilibrium']['rel_residual']:.1e} "
          f"psi_hat_c={rec['psi_hat_centre']:.3f} B1: G_DD={r['B1']['G_DD']:.10e} G_SD={r['B1']['G_SD']:.10e} G_DB={r['B1']['G_DB']:.4e} "
          f"| B0: G0={r['B0']['G_SD']:.10e} t={rec['t_s']:.0f}s", flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
