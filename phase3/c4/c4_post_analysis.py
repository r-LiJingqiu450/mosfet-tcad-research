"""C-4 analysis of the l_res diagnostic (fb4) and of TC-1. Applies the rules registered in
phase3/specs/C4_T1v2_Preregistration_2026-10-30.md (section 4) and
phase3/reports/C4_Small_Bias_Perturbation_Analysis_2026-10-30.md (section 4). Post-hoc
quantities are labelled as such in the output.
Usage: python phase3/c4/c4_post_analysis.py --fb4 DIR --tc1 DIR --out OUT
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from phase3 import meshmap as MM
from phase3.c3.t1_analysis import classify


def model_key(cfg):
    return cfg['model'] + (f"_f{cfg['lres']:g}" if cfg['model'] == 'c' else '')


def fb4_analysis(d):
    recs = {}
    for p in sorted(glob.glob(os.path.join(d, 'fb4_*.json'))):
        r = json.load(open(p)); c = r['config']
        recs.setdefault(model_key(c), {})[c['level']] = r
    dm = MM.f0_depth_map()
    h = [1.2e-4 / (len(dm.nodes(k)) - 1) for k in (3, 2, 1, 0)]
    out, rows = {}, []
    for m, byl in sorted(recs.items()):
        levs = sorted(byl)
        for k in levs:
            s = byl[k]['result']
            rows.append(dict(model=m, level=k, eps_h=s['eps_h'], eps_M0=s['eps_M0'], I_2D=s['I_A_per_cm'],
                             I_GCA_h=s['I_GCA_h'], quad_est_rel=s['I_GCA_h_quad_est'] / s['I_GCA_h'],
                             delta_cont=s['delta_cont'], delta_KCL=s['delta_KCL'], delta_terminal=s['delta_terminal'],
                             final_dpsi=s['final_update']['dpsi'], exploratory=byl[k]['config'].get('exploratory', False)))
        o = dict(levels=levs, eps_h=[byl[k]['result']['eps_h'] for k in levs])
        if levs == [0, 1, 2, 3]:
            f = [byl[k]['result']['eps_h'] for k in (3, 2, 1, 0)]
            c = classify(f, h)
            o.update(cls=c['cls'], p123=c['p123'], p234=c['p234'], U_disc=c['U'],
                     eps_extrap=c.get('richardson'), eps_finest=f[0])
            fm = [byl[k]['result']['eps_M0'] for k in (3, 2, 1, 0)]
            cm = classify(fm, h)
            o.update(M0_cls=cm['cls'], M0_eps_extrap=cm.get('richardson'), M0_U=cm['U'])
        out[m] = o
    # pairwise sensitivity S(X,Y) = I_X/I_Y - 1 at each level (sensitivity between two specified models)
    pairs = {}
    ms = sorted(recs)
    for i, x in enumerate(ms):
        for y in ms[i + 1:]:
            pairs[f'{x}|{y}'] = {k: recs[x][k]['result']['I_A_per_cm'] / recs[y][k]['result']['I_A_per_cm'] - 1
                                 for k in recs[x] if k in recs[y]}
    return out, rows, pairs


def tc1_analysis(d):
    res = {}
    for p in sorted(glob.glob(os.path.join(d, 'tc1_*.json'))):
        r = json.load(open(p)); c = r['config']
        for e, v in r['ends'].items():
            res.setdefault((c['L'], e), {})[c['level']] = v['fit']
    table = []
    for (L, e), byl in sorted(res.items()):
        for k, f in sorted(byl.items()):
            table.append(dict(L_um=L, ends=e, level=k, **f))
    # registered decision rule
    Ls = sorted({L for L, _ in res})
    sym_ok = all(res[(L, e)][k]['c1_over_s1'] <= 3 for L in Ls for e in ('aa', 'bb') for k in res[(L, e)])
    asym_ok = all(res[(L, e)][k]['c1_over_s1'] > 10 for L in Ls for e in ('ab', 'ba') for k in res[(L, e)])
    conv = {}
    for L in Ls:
        for e in ('ab', 'ba'):
            byl = res[(L, e)]
            if 1 in byl and 2 in byl:
                conv[f'{L:g}:{e}'] = abs(byl[2]['c1'] - byl[1]['c1']) / abs(byl[2]['c1'])
    conv_ok = all(v <= 0.10 for v in conv.values()) and len(conv) == 2 * len(Ls)
    sym_bad = [(L, e) for L in Ls for e in ('aa', 'bb') if sum(res[(L, e)][k]['c1_over_s1'] > 10 for k in res[(L, e)]) >= 2]
    asym_none = all(res[(L, e)][k]['c1_over_s1'] <= 3 for L in Ls for e in ('ab', 'ba') for k in res[(L, e)])
    if sym_ok and asym_ok and conv_ok:
        verdict = 'SUPPORTED'
    elif sym_bad or asym_none:
        verdict = 'FALSIFIED'
    else:
        verdict = 'INCONCLUSIVE'
    # post-hoc (not part of the rule): size of the symmetric-deck linear term relative to the mixed decks,
    # end-additivity of c2, c1(ab)+c1(ba), and the L-scaling of c1 between L = 1 and 3 um
    post = {}
    for L in Ls:
        for k in (0, 1, 2):
            try:
                g = {e: res[(L, e)][k] for e in ('aa', 'bb', 'ab', 'ba')}
            except KeyError:
                continue
            post[f'L{L:g}_lev{k}'] = dict(
                max_sym_over_asym_c1=max(abs(g['aa']['c1']), abs(g['bb']['c1'])) / min(abs(g['ab']['c1']), abs(g['ba']['c1'])),
                c1_ab_plus_ba_rel=(g['ab']['c1'] + g['ba']['c1']) / abs(g['ab']['c1']),
                c2_additivity=((g['ab']['c2'] + g['ba']['c2']) - (g['aa']['c2'] + g['bb']['c2'])) / abs(g['aa']['c2']),
                c2=dict((e, g[e]['c2']) for e in g), c1=dict((e, g[e]['c1']) for e in g))
    if len(Ls) >= 2:
        for k in (0, 1, 2):
            try:
                post[f'c1_Lslope_lev{k}_ab'] = math.log(abs(res[(Ls[1], 'ab')][k]['c1'] / res[(Ls[0], 'ab')][k]['c1'])) / math.log(Ls[1] / Ls[0])
                post[f'c2_Lslope_lev{k}_aa'] = math.log(abs(res[(Ls[1], 'aa')][k]['c2'] / res[(Ls[0], 'aa')][k]['c2'])) / math.log(Ls[1] / Ls[0])
            except (KeyError, ValueError, ZeroDivisionError):
                pass
    return dict(verdict=verdict, rule=dict(symmetric_all_le_3sigma=sym_ok, asymmetric_all_gt_10sigma=asym_ok,
                                           asym_c1_convergence_lev1_to_2=conv, convergence_ok=conv_ok,
                                           symmetric_resolved_at_2plus_levels=sym_bad, asymmetric_none_resolved=asym_none),
                post_hoc=post), table


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--fb4'); ap.add_argument('--tc1'); ap.add_argument('--out', required=True)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    summ = {}
    if a.fb4:
        o, rows, pairs = fb4_analysis(a.fb4)
        summ['lres_diagnostic'] = dict(models=o, pairwise_S=pairs)
        with open(os.path.join(a.out, 'lres_diagnostic.csv'), 'w', newline='') as fh:
            fh.write('# l_res diagnostic: L = 1 um, V_G = 1.25 V, V_D = 50 mV; eps_h = I_GCA,h/I_X,h - 1\n')
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
        try:
            import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
            fig, ax = plt.subplots(1, 2, figsize=(11, 4))
            for m, v in o.items():
                ax[0].plot(v['levels'], v['eps_h'], marker='o', label=m)
            ax[0].set(xlabel='mesh level', ylabel='eps_h at V_D = 50 mV', title='Same-mesh GCA discrepancy, L = 1 um, V_G = 1.25 V')
            ax[0].axhline(0, c='k', lw=.6); ax[0].legend(fontsize=7)
            cs = sorted([(float(m.split('_f')[1]), v) for m, v in o.items() if m.startswith('c_f') and v.get('eps_extrap') is not None])
            if cs:
                ax[1].errorbar([x for x, _ in cs], [v['eps_extrap'] for _, v in cs], yerr=[2 * (v['U_disc'] or 0) for _, v in cs],
                               marker='s', color='C2', ls='', capsize=3, label='R-I(c), extrapolated (bars 2 U_disc)')
            na = sorted([(float(m.split('_f')[1]), v) for m, v in o.items() if m.startswith('c_f') and v.get('eps_extrap') is None and 'eps_finest' in v])
            if na:
                ax[1].errorbar([x for x, _ in na], [v['eps_finest'] for _, v in na], yerr=[2 * (v['U_disc'] or 0) for _, v in na],
                               marker='o', mfc='none', color='C2', ls='', capsize=3,
                               label='R-I(c), finest level, non-asymptotic (bars 2 U, exploratory f=3)')
            for m, col in (('a', 'C0'), ('b', 'C1')):
                if m in o and o[m].get('eps_extrap') is not None:
                    ax[1].axhline(o[m]['eps_extrap'], color=col, ls='--', label=f'R-I({m}) extrapolated')
            ax[1].set(xscale='log', xlabel='l_res / W_dep', ylabel='eps', title='Dependence on reservoir-strip length (L = 1 um, 1.25 V, 50 mV)')
            ax[1].axhline(0, c='k', lw=.6); ax[1].legend(fontsize=7)
            fig.tight_layout(); fig.savefig(os.path.join(a.out, 'lres_diagnostic.png'), dpi=150); fig.savefig(os.path.join(a.out, 'lres_diagnostic.svg'))
        except Exception as exc:
            summ['figure_error'] = repr(exc)
    if a.tc1:
        t, table = tc1_analysis(a.tc1)
        summ['TC1'] = t
        with open(os.path.join(a.out, 'tc1_fits.csv'), 'w', newline='') as fh:
            fh.write('# TC-1 fits r(d) = c1 d + c2 d^2 + c3 d^3 (V^-1, V^-2, V^-3); sigma scaled by max(1, chi2/dof)\n')
            w = csv.DictWriter(fh, fieldnames=list(table[0].keys())); w.writeheader(); w.writerows(table)
    json.dump(summ, open(os.path.join(a.out, 'c4_post_summary.json'), 'w'), indent=1, default=str)
    print(json.dumps(summ, default=str, indent=1)[:6000])
    return 0


if __name__ == '__main__':
    sys.exit(main())
