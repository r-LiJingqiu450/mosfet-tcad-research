"""C-7 analysis: realism test (registered rules S1/S0) and M7 pre-registered test (rules R1-R4).
Outputs results/phase3/C7/analysis/{c7_summary.json, realism_table.csv, m7_test_table.csv, C7_figure.png}."""
from __future__ import annotations

import csv
import glob
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np

OUT = 'results/phase3/C7/analysis'
WDEP = 304.0   # nm, F0 (N_A = 1e16) planning value used in the registration


def load(d):
    rows = []
    for p in sorted(glob.glob(os.path.join(d, 'e7_*.json'))):
        r = json.load(open(p))
        if 'end_utc' not in r:
            continue
        a, c, m = r['args'], r['result']['const_mu'], r['result']['lombardi_mu']
        rows.append(dict(file=os.path.basename(p), deck=a['deck'], VG=a['vg'], tox=a['tox'], NA=a['na'], ov=a['ov'],
                         sigma=a['sigma'], dpsi=a['psibi_shift'], L=a['L'], level=a['level'],
                         dL_net=c['dL_net_nm'], dL_slice=c['dL_slice_nm'], dL_spread=c['dL_spread_nm'], dL_rows=c['dL_rows_nm'],
                         dL_top20=c.get('dL_slice_top20nm_nm'), dL_top50=c.get('dL_slice_top50nm_nm'),
                         dL_mu=m['dL_net_nm'], dL_mu_slice=m['dL_slice_nm'], bounds=c['bounds_ok'] and m['bounds_ok'],
                         Ns=r['result']['Ns_cm2'], newton_it=r['equilibrium']['final']['iterations'],
                         final_update=r['equilibrium']['final']['last'][-1], wall_s=r['wall_s']))
    return rows


def realism(rows):
    R = {(z['VG'], z['ov'], z['sigma']): z for z in rows if z['deck'] == 'R' and z['L'] == 2 and z['level'] == 0}
    comb = [(R[(vg, 30.0, 30.0)]) for vg in (1.0, 2.0)]
    s1 = dict(dL_comb_1V=comb[0]['dL_net'], dL_comb_1V_mu=comb[0]['dL_mu'], ratio_1V_2V=comb[0]['dL_net'] / comb[1]['dL_net'],
              ratio_1V_2V_mu=comb[0]['dL_mu'] / comb[1]['dL_mu'])
    s1['S1'] = bool(s1['dL_comb_1V'] >= 0.25 * WDEP and s1['dL_comb_1V_mu'] >= 0.25 * WDEP
                    and s1['ratio_1V_2V'] >= 2 and s1['ratio_1V_2V_mu'] >= 2)
    real = [R[(1.0, ov, sg)] for ov in (10.0, 30.0) for sg in (10.0, 30.0) if (1.0, ov, sg) in R]
    s1['S0_any'] = any(z['dL_net'] < 0.1 * WDEP or z['dL_mu'] < 0.1 * WDEP for z in real)
    s1['S0_cases_available'] = [(z['ov'], z['sigma']) for z in real]
    s1['verdict'] = 'SURVIVES (S1)' if s1['S1'] else ('DOES NOT SURVIVE (S0)' if s1['S0_any'] else 'INTERMEDIATE')
    return s1


def m7(rows):
    P = json.load(open('results/phase3/C7/m7/m7_predictions.json'))
    out = []
    for p in P['predictions']:
        z = [r for r in rows if r['deck'] == 'W' and abs(r['VG'] - p['VG']) < 1e-9 and r['tox'] == p['tox_nm'] and r['NA'] == p['N_A']]
        if not z:
            out.append(dict(**p, dL_2D=None)); continue
        z = z[0]
        out.append(dict(state=p['state'], label=p['level_label'], tox=p['tox_nm'], NA=p['N_A'], VG=p['VG'],
                        Ns_pred=p['Ns_1D'], Ns_2Dref=z['Ns'], dL_2D=z['dL_net'], dL_2D_slice=z['dL_slice'], dL_mu=z['dL_mu'],
                        dL_M7=p['dL_M7_nm'], dL_Ward=p['dL_Ward_nm'], r_M7=z['dL_net'] / p['dL_M7_nm'],
                        r_Ward=z['dL_net'] / p['dL_Ward_nm'], alpha_inf=p['alpha_inf']))
    done = [o for o in out if o.get('dL_2D') is not None]
    v = {}
    mod = [o for o in done if o['label'] in ('moderate', 'strong')]
    if len(mod) == 6:
        e = [abs(o['r_M7'] - 1) for o in mod]
        v['R1'] = 'SUPPORTED' if max(e) <= 0.15 else ('REFUTED' if max(e) > 0.25 else 'PARTIAL')
        v['R1_max_abs_dev'] = max(e)
    thr = [o for o in done if o['label'] in ('threshold', 'unseen bias')]
    if len(thr) == 4:
        inb = all((1.1 <= o['r_M7'] <= 1.6) if o['label'] == 'threshold' else (1.0 <= o['r_M7'] <= 1.4) for o in thr)
        refuted = any(not (0.9 <= o['r_M7'] <= 1.8) for o in thr)
        v['R2'] = 'REFUTED' if refuted else ('AS PREDICTED' if inb else 'WITHIN WIDE BAND ONLY')
    if len(done) == 10:
        sm = sum(abs(math.log(o['r_M7'])) for o in done); sw = sum(abs(math.log(o['r_Ward'])) for o in done)
        v['R3'] = dict(sum_abs_log_M7=sm, sum_abs_log_Ward=sw, M7_better=sm < sw)
        g = {(o['state'], o['label']): o for o in done}
        t = g[('S2_tox40', 'threshold')]['dL_2D'] / g[('S1_tox10', 'threshold')]['dL_2D']
        n = g[('S3_NA4e16', 'threshold')]['dL_2D'] / 369.873
        v['R4'] = dict(tox_ratio_2D=t, tox_ok=0.8 <= t <= 1.25, NA_ratio_2D=n, NA_ok=0.4 <= n <= 0.7)
    return out, v


def figure(rows, mt):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    C = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100']; INK, MUTED, GRID = '#1f1f1e', '#6b6a65', '#e6e5e0'
    plt.rcParams.update({'font.size': 9, 'axes.grid': True, 'grid.color': GRID, 'axes.spines.top': False,
                         'axes.spines.right': False, 'axes.edgecolor': MUTED, 'lines.linewidth': 2})
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.8))
    x = ax[0]
    R = {(z['VG'], z['ov'], z['sigma']): z for z in rows if z['deck'] == 'R' and z['L'] == 2 and z['level'] == 0}
    for (ov, sg), col, mk, lab in (((0., 0.), C[0], 'o', 'abrupt, no overlap'), ((30., 0.), C[1], 's', 'overlap 30 nm'),
                                   ((0., 30.), C[2], 'D', 'graded σ = 30 nm'), ((30., 30.), C[3], '^', 'overlap 30 + graded 30')):
        vg = [v for v in (1.0, 1.25, 2.0) if (v, ov, sg) in R]
        x.plot(vg, [R[(v, ov, sg)]['dL_net'] for v in vg], marker=mk, color=col, ms=7, label=lab)
        x.plot(vg, [R[(v, ov, sg)]['dL_mu'] for v in vg], ls=':', marker=mk, mfc='none', color=col, ms=6, lw=1.2)
    x.axhline(0, color=MUTED, lw=.8)
    x.set_xlabel('V_G (V)'); x.set_ylabel('ΔL per end (nm)')
    x.set_title('(a) Realism: solid = constant μ, dotted = Lombardi μ', loc='left', color=INK, fontsize=9)
    x.legend(frameon=False, fontsize=7.5)
    x = ax[1]
    d = [o for o in mt if o.get('dL_2D') is not None]
    lab = {'S0_original': ('t_ox 20, N_A 1e16 (new V_G)', 'o'), 'S1_tox10': ('t_ox 10 nm', 's'), 'S2_tox40': ('t_ox 40 nm', 'D'), 'S3_NA4e16': ('N_A 4e16', '^')}
    for k, (l, mk) in lab.items():
        e = [o for o in d if o['state'] == k]
        if e:
            x.plot([o['dL_2D'] for o in e], [o['dL_M7'] for o in e], mk, color=C[0], ms=7, label=f'M7: {l}')
            x.plot([o['dL_2D'] for o in e], [o['dL_Ward'] for o in e], mk, color=C[2], ms=6, mfc='none')
    lim = [10, 450]
    x.plot(lim, lim, color=MUTED, lw=1); x.fill_between(lim, [0.85 * v for v in lim], [1.15 * v for v in lim], color=GRID, alpha=.6)
    x.set_xscale('log'); x.set_yscale('log'); x.set_xlim(lim); x.set_ylim(lim)
    x.set_xlabel('2-D ΔL (nm), unseen states'); x.set_ylabel('predicted ΔL (nm)')
    x.set_title('(b) Pre-registered test: filled = M7, open = Ward §6.4', loc='left', color=INK, fontsize=9)
    x.legend(frameon=False, fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, 'C7_figure.png'), dpi=170)


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = load('results/phase3/C7/realism') + load('results/phase3/C7/m7')
    s = realism(rows)
    mt, v = m7(rows)
    json.dump(dict(realism=s, m7_test=mt, m7_verdicts=v, rows=rows), open(os.path.join(OUT, 'c7_summary.json'), 'w'), indent=1)
    for name, data in (('realism_table.csv', [r for r in rows if r['file'] and ('R_' in r['file'] or 'dpsi0.05' in r['file'] or 'dpsi-0.05' in r['file'])]),
                       ('m7_test_table.csv', [o for o in mt if o.get('dL_2D') is not None])):
        if data:
            with open(os.path.join(OUT, name), 'w', newline='') as fh:
                w = csv.DictWriter(fh, fieldnames=list(data[0].keys())); w.writeheader(); w.writerows(data)
    figure(rows, mt)
    print(json.dumps(s, indent=1)); print(json.dumps(v, indent=1))
    for o in mt:
        if o.get('dL_2D') is not None:
            print(f"{o['state']:12s} {o['label']:11s} VG {o['VG']:.4f} 2D {o['dL_2D']:7.1f} M7 {o['dL_M7']:7.1f} r {o['r_M7']:.3f} | Ward {o['dL_Ward']:6.1f} r {o['r_Ward']:.3f} | mu {o['dL_mu']:.1f}")


if __name__ == '__main__':
    main()
