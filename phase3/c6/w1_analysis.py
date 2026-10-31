"""Analysis of the C-6 W1 probe (phase3/c6/wprobe.py outputs).

For every state: dL_net (2-D truth), dL_slice (Ward's quasi-1-D current with the exact 2-D charge,
an upper bound on G), dL_spread = dL_slice - dL_net (>= 0, current spreading), dL_rows (lower bound).
Mechanism split of dL_slice (local, continuum, exact 1-D electrostatics with Boltzmann carriers):
  N_loc(x) = N_1D(psi_s(x)): the inversion sheet density of a 1-D slice whose surface potential is the
  2-D surface potential at x ('fringe' = surface-potential relaxation with 1-D vertical field).
  dL_fringe = int_0^{L/2} (1 - N_inf/N_loc(x)) dx;  dL_2Dfield = dL_slice - dL_fringe
  (the remainder: 2-D modification of the vertical field at given psi_s = charge sharing; Ward's alpha term).
Ward's assumed profile ([W88] p. 121): psi_s(x) = psi_s,inf + (psi_bi - psi_s,inf) exp(-pi x / (2 t_ox)).
Ward a-priori dL_W from phase3/c5/ward64_quadrature.py (evaluated here for every V_G present).
"""
from __future__ import annotations

import glob
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
from scipy.integrate import quad

from phase3.c5 import ward64_quadrature as W

Q = 1.602176634e-19
VT = 0.025851999786435316
EPS = 11.7 * 8.8541878128e-14
NI, NA = 1e10, 1e16
P0 = .5 * (NA + math.sqrt(NA * NA + 4 * NI * NI)); N0 = NI * NI / P0
PSI_B0 = -VT * math.log(P0 / NI)
TOX = 20e-7


def n_sheet_1d(psi_s):
    """Inversion sheet density (cm^-2) of a semi-infinite 1-D slice, frozen holes, phi_n = 0."""
    def E(psi):
        b = (psi - PSI_B0) / VT
        f = P0 * (math.exp(-b) + b - 1) + N0 * (math.expm1(b) - b)
        return math.sqrt(max(2 * Q * VT * f / EPS, 0.0))
    # N = int n dy = int_{psi_b}^{psi_s} n(psi) / E(psi) dpsi ; split to resolve the exponential near psi_s
    f = lambda p: NI * math.exp(p / VT) / E(p) if p > PSI_B0 + 1e-9 else 0.0
    pts = np.linspace(PSI_B0 + 1e-6, psi_s, 40)
    return sum(quad(f, a, b, limit=200)[0] for a, b in zip(pts[:-1], pts[1:]))


def load(d):
    out = []
    for p in sorted(glob.glob(os.path.join(d, 'w1_*.json'))):
        r = json.load(open(p))
        if 'end_utc' in r:
            out.append(r)
    return out


def per_state(r, cache):
    dc = r['decomposition']; a = r['args']
    x = np.array(dc['surface']['x_cm']); ps = np.array(dc['surface']['psi_V'])
    L = r['geometry']['n_lateral'] and x[-1]
    half = x <= L / 2 + 1e-15
    xm, pm = x[half], ps[half]
    psi_center = float(ps[np.argmin(np.abs(x - L / 2))])
    psi_cont = W.ws_from_vg(a['vg']) * W.LL * VT       # continuum 1-D surface potential (uniform channel)
    # reference far-field value: the discrete 2-D channel-centre value (end layers decayed for L >= 2 um;
    # flagged for L = 1 um), so that discretization offsets of the reference do not integrate over L
    psi_inf = psi_center
    N_inf = n_sheet_1d(psi_inf)
    # Ward-hybrid: Ward's Q_cs (p. 139) with v_s(x) replaced by the 2-D surface potential (alpha_0l kept)
    a0_of_yt = W.make_conv('parameter')
    yt = xm * 1e4 / W.EPSL_UM
    vs = pm / (VT * W.LL); wsl = psi_inf / (VT * W.LL)
    qm = W.qabs(wsl, W.A_INF)
    qcs = np.array([W.qabs(v, a0_of_yt(max(t, 1e-9))) for v, t in zip(vs, yt)])
    dL_hyb = float(np.trapezoid(1 - qm / qcs, xm)) * 1e7
    qcs_v = np.array([W.qabs(v, W.A_INF) for v in vs])                 # Ward profile channel only, 2-D psi
    dL_hyb_v = float(np.trapezoid(1 - qm / qcs_v, xm)) * 1e7
    Nloc = np.array([n_sheet_1d(p) for p in pm])
    integrand = 1 - N_inf / Nloc
    dL_fringe = float(np.trapezoid(integrand, xm)) * 1e7
    # Ward's assumed profile vs 2-D surface potential (source half)
    pw = psi_inf + (W_PSIBI - psi_inf) * np.exp(-math.pi * xm / (2 * TOX))
    int2d = float(np.trapezoid(pm - psi_inf, xm)) * 1e7
    intw = float(np.trapezoid(pw - psi_inf, xm)) * 1e7
    # 1/e decay length of (psi_s - psi_inf) from its value at x = 0
    dev = (pm - psi_inf) / (pm[0] - psi_inf)
    k = np.where(dev < math.exp(-1))[0]
    l_e = float(xm[k[0]]) * 1e7 if len(k) else float('nan')
    ws_lnlam = (psi_inf / VT) / math.log(NA / NI)
    return dict(L=a['L'], level=a['level'], VG=a['vg'], hw=a['hw'],
                dL_net=dc['dL_net_nm'], dL_slice=dc['dL_slice_nm'], dL_spread=dc['dL_spread_nm'], dL_rows=dc['dL_rows_nm'],
                dL_fringe_local=dL_fringe, dL_2Dfield=dc['dL_slice_nm'] - dL_fringe,
                dL_Ward_hybrid=dL_hyb, dL_Ward_hybrid_psi_only=dL_hyb_v, reference_ok=a['L'] >= 2,
                psi_center_minus_continuum=psi_center - psi_cont,
                G2D_over_GCA=dc['G_net'] / dc['G_GCA'], Ns_GCA=dc['Ns_GCA_cm2'],
                dL_net_times_Ns=dc['dL_net_nm'] * dc['Ns_GCA_cm2'],
                psi_s_inf=psi_inf, psi_s_center_minus_inf=psi_center - psi_inf, ws_lnlam=ws_lnlam, psi_s_edge=float(pm[0]),
                surf_excess_int_2D_nmV=int2d, surf_excess_int_Ward_nmV=intw, surf_decay_1e_nm=l_e,
                ward_profile_decay_nm=2 * TOX / math.pi * 1e7,
                cross_net_vs_devsim=r['cross_check']['net_vs_devsim_rel'] if r['devsim'] else float('nan'),
                bounds_ok=dc['bounds_ok'],
                devsim_E=(r['devsim']['E_dwr'] + r['devsim']['E_round']) if r['devsim'] else float('nan'), wall_s=r['wall_s'])


W_PSIBI = VT * math.log(1e20 / NI)


def main(d='results/phase3/C6/w1', out='results/phase3/C6/analysis_w1'):
    os.makedirs(out, exist_ok=True)
    rows = [per_state(r, {}) for r in load(d)]
    rows.sort(key=lambda z: (z['hw'], z['VG'], z['L'], z['level']))
    # Ward predictions at the V_G present
    ward = {}
    for vg in sorted({z['VG'] for z in rows}):
        ws = W.ws_from_vg(vg)
        ward[vg] = dict(ws_lnlam=ws * W.LL, dL_W=W.dL_end_nm(ws) if ws > 1 else None,
                        dL_W_fringe=W.dL_end_nm(ws, mech='v') if ws > 1 else None,
                        dL_W_alpha=W.dL_end_nm(ws, mech='alpha') if ws > 1 else None)
    for z in rows:
        w = ward[z['VG']]
        z['dL_Ward'] = w['dL_W']; z['rho_net_over_Ward'] = z['dL_net'] / w['dL_W'] if w['dL_W'] else None
    json.dump(dict(rows=rows, ward=ward), open(os.path.join(out, 'w1_summary.json'), 'w'), indent=1)
    import csv
    with open(os.path.join(out, 'w1_table.csv'), 'w', newline='') as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
    hdr = 'hw  VG    L  lev  G2D/GCA   dL_net  dL_slice spread  dL_rows  fringe 2Dfield  Ward  hybrid hyb_psi rho   l_e(nm) net-vs-dev'
    print(hdr)
    for z in rows:
        print(f"{z['hw']:.1f} {z['VG']:4.2f} {z['L']:4g} {z['level']:3d} {z['G2D_over_GCA']:8.5f} {z['dL_net']:8.3f} {z['dL_slice']:8.3f} "
              f"{z['dL_spread']:6.3f} {z['dL_rows']:8.3f} {z['dL_fringe_local']:7.2f} {z['dL_2Dfield']:7.2f} "
              f"{(z['dL_Ward'] or float('nan')):6.1f} {z['dL_Ward_hybrid']:6.1f} {z['dL_Ward_hybrid_psi_only']:6.1f} {(z['rho_net_over_Ward'] or float('nan')):5.2f} {z['surf_decay_1e_nm']:7.1f} {z['cross_net_vs_devsim']:.1e}")
    return rows, ward


if __name__ == '__main__':
    main(*sys.argv[1:])
