"""C-10: far-field (tail) analysis of the zero-bias subthreshold end correction on stored 2-D states.

For every stored endprobe/probe9 state that carries the surface-potential profile psi_s(x) (y = 0), this computes
  - the far-field surface potential psi_inf (channel centre) and u(x) = psi_s(x) - psi_inf;
  - an exponential tail fit u(x) ~ B exp(-x/lam) on 1e-4 V < u < 0.5 V_T, x < L/2 - 4 lam  (lam_fit, B_fit);
  - the asymptotic end correction (C-10 identity, see the C-10 theory note)
        dL_asym = lam_fit * [Ein(B_fit/V_T) + r],   Ein(z) = gamma + ln z + E1(z),   r = C_ox V_T / Q_dep,inf
    and the surface-only integral dL_surf = int_0^{L/2} (1 - exp(-u/V_T)) dx;
  - the 1-D far-field quantities (exact Boltzmann 1-D charge, depletion width, two-layer eigen-length l0) from
    phase3.c7.model_m7.Device, and the similarity variables of the C-10 note:
        t_hat = (eps_si/eps_ox) t_ox / L_D,  h_hat = h_w / L_D,  psi_hat = (psi - psi_b0)/V_T,  L_D = sqrt(eps V_T/(q N_A)).
Writes results/phase3/C10/tail/tail_<tag>.json and prints a table.
Usage: python phase3/c10/tail_analysis.py --glob 'results/phase3/C8/subth/e7_*.json' --tag c8subth
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np  # noqa: E402
from scipy.special import exp1  # noqa: E402

from phase3.c7.model_m7 import Device, VT, Q, ESI, EOX  # noqa: E402

GAMMA = 0.5772156649015329


def ein(z):
    return GAMMA + math.log(z) + float(exp1(z))


def far(r):
    a = r['args']; dev = Device(a['na'], a['tox'])
    psig, pinf1, Wd = dev.far_field(a['vg'])
    l0 = dev.l0(Wd)
    b = pinf1 - dev.psib
    qdep = math.sqrt(2 * Q * ESI * a['na'] * max(b - VT, 1e-9))
    LD = math.sqrt(ESI * VT / (Q * a['na']))
    alpha = dev.ainf * math.sqrt(dev.ll) * dev.lam ** ((1 - pinf1 / (VT * dev.ll)) / 2)
    hw = r.get('geometry', {}).get('h_w_cm', a.get('hw', 0.5) * 1e-4)
    return dict(dev=dev, psi_inf_1d=pinf1, Wd_nm=Wd * 1e7, l0_nm=l0 * 1e7, qdep=qdep, cox=dev.cox, LD_nm=LD * 1e7,
                alpha=alpha, t_hat=(ESI / EOX) * dev.tox / LD, h_hat=hw / LD, psi_hat=b / VT,
                psibi_hat=(dev.psibi - dev.psib) / VT, r=dev.cox * VT / qdep)


def tail(r, fr):
    s = r['result']['surface']; x = np.array(s['x_cm']); p = np.array(s['psi_V'])
    L = r['args']['L'] * 1e-4
    m = (x >= -1e-15) & (x <= L / 2 + 1e-12)
    x, p = x[m], p[m]
    pinf = float(p[-1])                       # value at the channel centre (x = L/2 is a node of the symmetric map)
    u = p - pinf
    lam_guess = fr['l0_nm'] * 1e-7
    sel = (u > 1e-4) & (u < 0.5 * VT) & (x < L / 2 - 4 * lam_guess)
    out = dict(psi_inf_2d=pinf, psi_s0=float(p[0]), n_fit=int(sel.sum()))
    if sel.sum() >= 4:
        A = np.polyfit(x[sel], np.log(u[sel]), 1)
        lam = -1.0 / A[0]; B = math.exp(A[1])
        resid = float(np.max(np.abs(np.polyval(A, x[sel]) - np.log(u[sel]))))
        out.update(lam_fit_nm=lam * 1e7, B_fit=B, B_over_VT=B / VT, fit_max_log_resid=resid,
                   dL_asym_nm=lam * (ein(B / VT) + fr['r']) * 1e7, dL_asym_nor_nm=lam * ein(B / VT) * 1e7)
    xs = np.linspace(0, x[-1], 20001)
    out['dL_surf_nm'] = float(np.trapezoid(1 - np.exp(-np.interp(xs, x, u) / VT), xs)) * 1e7
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--glob', required=True, action='append'); ap.add_argument('--tag', required=True)
    ap.add_argument('--alpha-min', type=float, default=0.0)
    a = ap.parse_args()
    paths = sorted(set(p for g in a.glob for p in glob.glob(g)))
    rows = []
    for p in paths:
        r = json.load(open(p))
        if 'end_utc' not in r or 'surface' not in r.get('result', {}):
            continue
        fr = far(r)
        if fr['alpha'] < a.alpha_min:
            continue
        c = r['result']['const_mu']
        t = tail(r, fr)
        row = dict(file=os.path.relpath(p), deck=r['args']['deck'], na=r['args']['na'], tox=r['args']['tox'], L=r['args']['L'],
                   vg=r['args']['vg'], level=r['args'].get('level', 0), hw_um=r['args'].get('hw', 0.5), halo=r['args'].get('halo_n', 0.0),
                   dL_net_nm=c['dL_net_nm'], dL_slice_nm=c['dL_slice_nm'], dL_rows_nm=c['dL_rows_nm'], Ns_cm2=r['result']['Ns_cm2'],
                   **{k: v for k, v in fr.items() if k != 'dev'}, **t)
        row['psibi_minus_psis_VT'] = (fr['dev'].psibi - t['psi_inf_2d']) / VT
        row['kappa'] = row['dL_net_nm'] / row['l0_nm']; row['dL_over_Wd'] = row['dL_net_nm'] / row['Wd_nm']
        row['dL_over_LD'] = row['dL_net_nm'] / row['LD_nm']
        if 'B_fit' in t:
            row['beta'] = t['B_fit'] / (fr['dev'].psibi - t['psi_inf_2d'])
            row['asym_rel_err'] = t['dL_asym_nm'] / row['dL_net_nm'] - 1
        rows.append(row)
    os.makedirs('results/phase3/C10/tail', exist_ok=True)
    json.dump(dict(globs=a.glob, rows=rows), open(f'results/phase3/C10/tail/tail_{a.tag}.json', 'w'), indent=1)
    hdr = 'deck NA     tox  L   VG     alpha    psih  dpsi/VT  dL_net  dL_asym  err%   lam_fit l0    B/VT   beta  dL/Wd dL/LD  t_hat h_hat'
    print(hdr)
    for w in sorted(rows, key=lambda z: (z['na'], z['tox'], z['deck'], z['hw_um'], z['L'], z['vg'])):
        print(f"{w['deck']}  {w['na']:.0e} {w['tox']:4g} {w['L']:3g} {w['vg']:6.4f} {w['alpha']:8.2f} {w['psi_hat']:6.2f} {w['psibi_minus_psis_VT']:6.2f}"
              f"  {w['dL_net_nm']:7.2f} {w.get('dL_asym_nm', float('nan')):7.2f} {100 * w.get('asym_rel_err', float('nan')):6.2f}"
              f"  {w.get('lam_fit_nm', float('nan')):6.2f} {w['l0_nm']:6.2f} {w.get('B_over_VT', float('nan')):6.1f} {w.get('beta', float('nan')):5.2f}"
              f" {w['dL_over_Wd']:5.2f} {w['dL_over_LD']:6.2f} {w['t_hat']:5.2f} {w['h_hat']:6.1f}")


if __name__ == '__main__':
    main()
