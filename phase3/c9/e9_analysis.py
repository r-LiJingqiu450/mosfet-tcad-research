"""C-9 analysis of direct 2-D solves (no affine extrapolation): H9-1..H9-4 (pre-registered), Terada-Muta on
direct data, split-C-V from terminal charges vs the C-8 electron proxy, two-length differential mobility.
Output: results/phase3/C9/analysis/e9_result.json (the figure is made by e9_figures.py).

Post-hoc option (C-9 erratum E9-1): --halo-runs DIR replaces the halo states (halo_n > 0) by those in DIR (probe9 v2,
corrected Lombardi impurity term) and --tag TAG writes e9_result_TAG.json; without options the original analysis
(e9_result.json) is reproduced unchanged."""
from __future__ import annotations

import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
from scipy.interpolate import CubicSpline

Q, MU = 1.602176634e-19, 1350.0
OUT = 'results/phase3/C9/analysis'


def load(halo_dir=None):
    D = {}
    paths = glob.glob('results/phase3/C9/runs/e9_*.json')
    if halo_dir:
        paths = [p for p in paths if json.load(open(p))['args']['halo_n'] == 0] + glob.glob(os.path.join(halo_dir, 'e9_*.json'))
    for p in paths:
        r = json.load(open(p))
        if 'end_utc' not in r:
            continue
        a = r['args']; key = (a['deck'], a['halo_n'], a['L'])
        c, m, t = r['result']['const_mu'], r['result']['lombardi_mu'], r['terminal_charges']
        D.setdefault(key, {})[a['vg']] = dict(G=c['G_net'], Gmu=m['G_net'], GGCA=r['result']['G_GCA_1COL'],
                                              GGCAmu=m['G_GCA_1COL_mu'], dL=c['dL_net_nm'], dLmu=m['dL_net_nm'],
                                              Ns=r['result']['Ns_cm2'], QSD=t['Q_SD'], Qe=t['Q_electrons'])
    return D


def series(D, key, f):
    vs = sorted(D[key]); return np.array(vs), np.array([D[key][v][f] for v in vs])


def max_gm_vt(v, G, vmax):
    m = v <= vmax + 1e-12
    cs = CubicSpline(v[m], G[m]); f = np.linspace(v[m][0], v[m][-1], 6001); g = cs(f); gm = cs(f, 1)
    k = int(np.argmax(gm))
    return float(f[k] - g[k] / gm[k]), float(f[k]), bool(0 < k < len(f) - 1)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--halo-runs', default=None); ap.add_argument('--tag', default=None)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    D = load(args.halo_runs); res = {}
    # H9-1 affine law (direct)
    h1 = {}
    for deck in ('W', 'R'):
        v, a = series(D, (deck, 0.0, 0.5), 'dL'); _, b = series(D, (deck, 0.0, 2.0), 'dL'); _, c1 = series(D, (deck, 0.0, 1.0), 'dL')
        m = v >= 1.24
        h1[deck] = dict(VG=v.tolist(), dL_05=a.tolist(), dL_1=c1.tolist(), dL_2=b.tolist(), max_abs_diff_05_2=float(np.max(np.abs(a[m] - b[m]))))
    h1['verdict'] = 'SUPPORTED' if all(h1[d]['max_abs_diff_05_2'] <= 1.0 for d in ('W', 'R')) else 'REJECTED (see max_abs_diff)'
    res['H9_1'] = h1
    # H9-2 split-C-V: terminal charge vs electron-only charge (both direct 2-D, same L), V_lo = 1.08
    h2 = {}
    for deck in ('W', 'R'):
        for L in (0.5, 1.0, 2.0):
            key = (deck, 0.0, L); v, G = series(D, key, 'G'); _, QSD = series(D, key, 'QSD'); _, Qe = series(D, key, 'Qe')
            _, Ns = series(D, key, 'Ns'); _, GG = series(D, key, 'GGCA')
            Lcm = L * 1e-4; sig = GG * Lcm
            Qt = -(QSD - QSD[0]); Qel = -(Qe - Qe[0]); Qlong = Q * (Ns - Ns[0]) * Lcm
            mu_t = Lcm ** 2 * G[1:] / Qt[1:]; mu_e = Lcm ** 2 * G[1:] / Qel[1:]; mu_long = sig[1:] / (Q * (Ns[1:] - Ns[0]))
            h2[f'{deck}_L{L:g}'] = dict(VG=v[1:].tolist(), err_terminal=(mu_t / mu_long - 1).tolist(), err_electrons=(mu_e / mu_long - 1).tolist(),
                                        contact_share=((Qt - Qel)[1:] / Qt[1:]).tolist())
    k = h2['W_L1']['VG'].index(1.32)
    ev = {d: (h2[f'{d}_L1']['err_terminal'][k], h2[f'{d}_L1']['err_electrons'][k]) for d in ('W', 'R')}
    h2['at_1.32V_L1'] = ev
    h2['verdict'] = {d: ('EQUIVALENT' if abs(t - e) <= 0.2 * abs(e) else 'NOT EQUIVALENT') for d, (t, e) in ev.items()}
    res['H9_2'] = h2
    # two-length differential mobility (terminal charge), deck R, L = 1 and 2
    v, G1 = series(D, ('R', 0.0, 1.0), 'G'); _, G2 = series(D, ('R', 0.0, 2.0), 'G')
    _, Q1 = series(D, ('R', 0.0, 1.0), 'QSD'); _, Q2 = series(D, ('R', 0.0, 2.0), 'QSD')
    dR = (1 / G2 - 1 / G1) / 1e-4; dQ = (-(Q2 - Q2[0]) + (Q1 - Q1[0])) / 1e-4
    _, Ns = series(D, ('R', 0.0, 1.0), 'Ns')
    res['two_length_R'] = dict(VG=v[1:].tolist(), mu_diff=(1 / dR[1:] / dQ[1:]).tolist(), mu_true=MU,
                               rel_err=((1 / dR[1:] / dQ[1:]) / MU - 1).tolist())
    # H9-3 halo
    h3 = {}
    for hn in (5e17, 2e18):
        for L in (0.5, 1.0, 2.0):
            if ('R', hn, L) in D:
                v, dl = series(D, ('R', hn, L), 'dL'); _, dlm = series(D, ('R', hn, L), 'dLmu')
                h3[f'{hn:.0e}_L{L:g}'] = dict(VG=v.tolist(), dL=dl.tolist(), dL_mu=dlm.tolist())
    k = h3['5e+17_L1']['VG'].index(1.32)
    h3['verdict_5e17'] = 'SUPPORTED (dL<0)' if h3['5e+17_L1']['dL'][k] < 0 else 'REJECTED (dL>=0)'
    res['H9_3'] = h3
    # H9-4 max-gm threshold with Lombardi mobility (direct), deck W and R; reference: 1COL-mobility sheet conductance
    h4 = {}
    for deck, hn in (('W', 0.0), ('R', 0.0), ('R', 2e18), ('R', 5e17)):
        for vmax in (1.80, 2.00):
            vtinf = None; rows = []
            for L in (0.5, 1.0, 2.0):
                if (deck, hn, L) not in D:
                    continue
                v, Gm = series(D, (deck, hn, L), 'Gmu'); _, GGm = series(D, (deck, hn, L), 'GGCAmu')
                if vtinf is None:
                    vtinf = max_gm_vt(v, GGm, vmax)
                vt = max_gm_vt(v, Gm, vmax)
                rows.append(dict(L=L, VT=vt[0], VG_gm_max=vt[1], interior=vt[2], dVT_mV=(vt[0] - vtinf[0]) * 1e3))
            h4[f'{deck}_halo{hn:.0e}_vmax{vmax}'] = dict(VT_inf=vtinf, rows=rows)
    w = {r['L']: r['dVT_mV'] for r in h4['W_halo0e+00_vmax2.0']['rows']}
    ok = (-15 <= w[1.0] <= -3) and (1.7 <= w[0.5] / w[1.0] <= 2.3)
    h4['verdict'] = ('SUPPORTED' if ok else 'REJECTED') + f" (dVT(1um)={w[1.0]:.2f} mV, ratio={w[0.5] / w[1.0]:.2f})"
    res['H9_4'] = h4
    # Terada-Muta on direct data (regression of 1/G on L over L = 0.5, 1, 2), const mu and Lombardi
    tm = {}
    for deck, hn in (('R', 0.0), ('R', 2e18), ('W', 0.0)):
        for f in ('G', 'Gmu'):
            Ls = np.array([0.5, 1.0, 2.0]) * 1e-4
            vv = sorted(D[(deck, hn, 1.0)])
            lines = {}
            for vg in (1.32, 1.40, 1.56, 2.00):
                y = np.array([1 / D[(deck, hn, L)][vg][f] for L in (0.5, 1.0, 2.0)])
                A = np.polyfit(Ls, y, 1); lines[vg] = (A, float(np.max(np.abs(np.polyval(A, Ls) - y) / y)))
            out = []
            for v1, v2 in ((1.32, 1.56), (1.40, 2.00), (1.56, 2.00)):
                (a1, b1), (a2, b2) = lines[v1][0], lines[v2][0]
                Lx = (b2 - b1) / (a1 - a2); Rx = a1 * Lx + b1
                out.append(dict(V1=v1, V2=v2, dL_TM_nm=Lx / 2 * 1e7, Rsd_TM_ohm_um=Rx * 1e4))
            tm[f'{deck}_halo{hn:.0e}_{f}'] = dict(pairs=out, max_line_nonlinearity={str(k): v[1] for k, v in lines.items()})
    res['terada_muta_direct'] = tm
    if args.halo_runs:
        res['source'] = dict(no_halo='results/phase3/C9/runs', halo=args.halo_runs, note='post-hoc erratum E9-1 (Lombardi impurity term)')
    json.dump(res, open(os.path.join(OUT, f"e9_result{'_' + args.tag if args.tag else ''}.json"), 'w'), indent=1)
    print('H9-1', h1['verdict'], {d: round(h1[d]['max_abs_diff_05_2'], 3) for d in ('W', 'R')})
    print('H9-2', h2['verdict'], {d: tuple(round(x, 4) for x in t) for d, t in ev.items()})
    print('H9-3', h3['verdict_5e17'], {k2: [round(x, 1) for x in v2['dL'][:4]] for k2, v2 in h3.items() if k2.startswith(('5e', '2e'))})
    print('H9-4', h4['verdict'])
    for k2, v2 in h4.items():
        if k2 != 'verdict':
            print('  ', k2, 'VTinf', round(v2['VT_inf'][0], 4), [(r['L'], round(r['dVT_mV'], 2), r['interior']) for r in v2['rows']])
    for k2, v2 in tm.items():
        print('TM', k2, [(p['V1'], p['V2'], round(p['dL_TM_nm'], 1), round(p['Rsd_TM_ohm_um'], 1)) for p in v2['pairs']],
              'nonlin', round(max(v2['max_line_nonlinearity'].values()), 5))
    print('two-length R:', [(v_, round(e, 4)) for v_, e in zip(res['two_length_R']['VG'], res['two_length_R']['rel_err'])][:6])
    for kk in ('W_L1', 'R_L1', 'R_L2'):
        e = h2[kk]; print('split', kk, [(v_, round(t, 3), round(x, 3), round(s, 3)) for v_, t, x, s in zip(e['VG'], e['err_terminal'], e['err_electrons'], e['contact_share'])][:5])


if __name__ == '__main__':
    main()
