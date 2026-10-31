"""C-10 evaluation of the pre-registered tests P10-1..P10-6 (phase3/specs/C10_Preregistration_2026-10-31.md).

Reads the new states (results/phase3/C10/runs) and the existing C-8 subthreshold states (results/phase3/C8/subth),
uses phase3/c10/tail_analysis.py for the far-field (tail) quantities, and writes
results/phase3/C10/analysis/c10_result.json.
"""
from __future__ import annotations

import glob
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from phase3.c10 import tail_analysis as T  # noqa: E402
from phase3.c7.model_m7 import VT, Q, ESI  # noqa: E402

OUT = 'results/phase3/C10/analysis'
DESIGN = json.load(open('results/phase3/C10/design.json'))
SHIFT = DESIGN['shift_V']


def load(globs):
    rows = []
    for g in globs:
        for p in sorted(glob.glob(g)):
            r = json.load(open(p))
            if 'end_utc' not in r or 'surface' not in r.get('result', {}):
                continue
            fr = T.far(r); t = T.tail(r, fr); a = r['args']; c = r['result']['const_mu']
            # E10-1: the drain-side current (G_net, dL_net_nm) loses accuracy in deep subthreshold; the source-side
            # current is cancellation-free (verified against the swapped-BC solve, phase3/c10/probe10.py).
            G_src = -c['src_current']; Lcm = a['L'] * 1e-4
            dL_src = 0.5 * Lcm * (1 - r['result']['G_GCA_1COL'] / G_src) * 1e7
            row = dict(file=os.path.relpath(p), na=a['na'], tox=a['tox'], L=a['L'], vg=a['vg'], level=a['level'], hw=a['hw'],
                       hw_real_um=r['geometry']['h_w_cm'] * 1e4, shift=a.get('psibi_shift', 0.0), dL=dL_src,
                       dL_drain_prereg=c['dL_net_nm'], bounds_ok_drain=c['bounds_ok'],
                       bounds_ok_src=bool(c['G_rows'] * (1 - 1e-12) <= G_src <= c['G_slice'] * (1 + 1e-12)),
                       dL_slice=c['dL_slice_nm'], dL_rows=c['dL_rows_nm'], Ns=r['result']['Ns_cm2'],
                       converged=r['equilibrium']['final']['converged'], **{k: v for k, v in fr.items() if k != 'dev'}, **t)
            row['psibi_hat_eff'] = row['psibi_hat'] + row['shift'] / VT
            row['dpsi_VT'] = (fr['dev'].psibi + row['shift'] - t['psi_inf_2d']) / VT      # psi_bi - psi_s,inf (2-D centre)
            row['x_s_FS93_nm'] = math.sqrt(2 * ESI * row['dpsi_VT'] * VT / (Q * row['na'])) * 1e7
            row['ratio_FS93'] = row['dL'] / row['x_s_FS93_nm']
            row['dL_over_LD'] = row['dL'] / row['LD_nm']
            if 'B_fit' in t:
                row['identity_rel_err'] = t['dL_asym_nm'] / row['dL'] - 1
                row['kappa_fit'] = row['dL'] / t['lam_fit_nm']
                row['ln_B_VT'] = math.log(t['B_over_VT'])
                row['beta'] = t['B_fit'] / (row['dpsi_VT'] * VT)
            rows.append(row)
    return rows


def pick(rows, **kw):
    out = []
    for r in rows:
        ok = True
        for k, v in kw.items():
            if isinstance(v, float):
                ok &= abs(r[k] - v) <= 1e-6 * max(1.0, abs(v))
            else:
                ok &= r[k] == v
        if ok:
            out.append(r)
    return out


def one(rows, **kw):
    m = pick(rows, **kw)
    assert len(m) == 1, (kw, [x['file'] for x in m])
    return m[0]


def main():
    os.makedirs(OUT, exist_ok=True)
    new = load(['results/phase3/C10/runs/e7_*.json'])
    old = load(['results/phase3/C8/subth/e7_*.json'])
    allr = new + old
    res = dict(n_new=len(new), failed_or_unconverged=[r['file'] for r in new if not r['converged']],
               E10_1=dict(note='dL from the source-side current (cancellation-free); dL_drain_prereg is the pre-registered '
                               'quantity computed from the drain-side current, numerically invalid where bounds_ok_drain is False',
                          drain_side_bound_violations=[(r['file'], r['dL_drain_prereg'], r['dL']) for r in allr if not r['bounds_ok_drain']],
                          src_side_bound_violations=[r['file'] for r in allr if not r['bounds_ok_src']]))
    vg1, vg2 = 0.7835, round(DESIGN['vg_deep_device1'], 6)
    hw2s, hw1s = DESIGN['hw_shallow_um'], DESIGN['hw_device1_similarity_um']
    # P10-1 similarity image
    p1 = []
    for vg in (vg1, vg2):
        d = one(new, na=1e18, vg=vg, hw=hw2s, level=0, tox=2.0, shift=round(SHIFT, 6))
        s = one(new, na=1e16, vg=vg, hw=hw1s, level=0)
        ratio = (d['dL'] / d['LD_nm']) / (s['dL'] / s['LD_nm'])
        p1.append(dict(vg=vg, psi_hat_2prime=d['psi_hat'], psi_hat_1=s['psi_hat'], h_hat_2prime=d['hw_real_um'] * 1e3 / d['LD_nm'],
                       h_hat_1=s['hw_real_um'] * 1e3 / s['LD_nm'], dL_2prime=d['dL'], dL_1=s['dL'], dLLD_2prime=d['dL'] / d['LD_nm'],
                       dLLD_1=s['dL'] / s['LD_nm'], ratio=ratio, alpha_1=s['alpha'], alpha_2prime=d['alpha']))
    res['P10_1'] = dict(rows=p1, verdict='SUPPORTED' if all(abs(r['ratio'] - 1) <= 0.02 for r in p1) else 'REJECTED')
    # post-hoc (not pre-registered): exact image with n_i scaled to 1e12 so that n_face/N_A = 1e4 as in device 1
    # (phase3/c10/probe10ni.py); only dL (source-side current) and L_D are used, both independent of n_i bookkeeping
    ph = []
    for vg in (vg1, vg2):
        f = glob.glob(f'results/phase3/C10/ni_image/e7_W_L1_lev0_VG{vg:g}_*.json')
        if not f:
            continue
        r = json.load(open(f[0])); c = r['result']['const_mu']
        dl = 0.5 * 1e-4 * (1 - r['result']['G_GCA_1COL'] / (-c['src_current'])) * 1e7
        s = one(new, na=1e16, vg=vg, hw=hw1s, level=0)
        LD2 = math.sqrt(ESI * VT / (Q * 1e18)) * 1e7
        ph.append(dict(vg=vg, dL_image=dl, dLLD_image=dl / LD2, dLLD_1=s['dL'] / s['LD_nm'], ratio=(dl / LD2) / (s['dL'] / s['LD_nm']),
                       bounds_ok=bool(c['G_rows'] * (1 - 1e-12) <= -c['src_current'] <= c['G_slice'] * (1 + 1e-12))))
    res['P10_1_posthoc_ni_image'] = dict(rows=ph, note='post-hoc; n_i = 1e12, N_A = 1e18, psi_bi face value unchanged (N_D = 1e22): '
                                                       'every dimensionless parameter of device 1 matched; tolerance not pre-registered')
    # P10-2 factor decomposition at vg1
    dev1_hw05 = one(old, na=1e16, vg=vg1)
    dec = []
    for tag, hw, sh in (('A as is', 0.5, 0.0), ('B psi_bi matched', 0.5, round(SHIFT, 6)), ('C h_w matched', hw2s, 0.0),
                        ("D device 2' (both)", hw2s, round(SHIFT, 6))):
        r = one(new, na=1e18, vg=vg1, hw=hw, shift=sh, level=0, tox=2.0)
        dec.append(dict(tag=tag, dL=r['dL'], dLLD=r['dL'] / r['LD_nm'], rel_to_dev1_hw05=(r['dL'] / r['LD_nm']) / (dev1_hw05['dL'] / dev1_hw05['LD_nm']) - 1,
                        h_hat=r['hw_real_um'] * 1e3 / r['LD_nm'], psibi_hat=r['psibi_hat_eff'], beta=r.get('beta'), lam_fit=r.get('lam_fit_nm')))
    res['P10_2'] = dict(device1_hw05=dict(dL=dev1_hw05['dL'], dLLD=dev1_hw05['dL'] / dev1_hw05['LD_nm'], h_hat=dev1_hw05['hw_real_um'] * 1e3 / dev1_hw05['LD_nm']),
                        rows=dec)
    # P10-3 oxide scaling
    base305 = one(old, na=1e18, vg=1.08); baseA = one(new, na=1e18, vg=vg1, hw=0.5, shift=0.0, level=0, tox=2.0)
    p3 = []
    for s in DESIGN['tox_states']:
        r = one(new, na=1e18, tox=s['tox'], vg=round(s['vg'], 6), level=0)
        b = base305 if s['psi_hat'] > 25 else baseA
        key = f"{s['tox']:g}nm_psi{s['psi_hat']:.2f}"
        l0r = DESIGN['l0_ratio_pred'][key]
        ratio = r['dL'] / b['dL']
        p3.append(dict(key=key, dL=r['dL'], dL_base=b['dL'], ratio=ratio, l0_ratio=l0r, lamfit_ratio=r['lam_fit_nm'] / b['lam_fit_nm'],
                       in_Hl0=abs(ratio / l0r - 1) <= 0.10, in_HWd=abs(ratio - 1) <= 0.05, kappa=r['dL'] / r['l0_nm'], kappa_base=b['dL'] / b['l0_nm'],
                       dL_over_Wd=r['dL'] / r['Wd_nm'], dL_over_Wd_base=b['dL'] / b['Wd_nm']))
    six = [x for x in p3 if x['key'].startswith('6nm')]
    v3 = ('H-l0 SUPPORTED, H-Wd REJECTED' if all(x['in_Hl0'] and not x['in_HWd'] for x in six) else
          'H-Wd SUPPORTED, H-l0 REJECTED' if all(x['in_HWd'] and not x['in_Hl0'] for x in six) else 'INCONCLUSIVE')
    res['P10_3'] = dict(rows=p3, verdict=v3)
    # P10-4 bias range
    p4 = {}
    for na, tox, L in ((1e18, 2.0, 1.0), (1e16, 20.0, 2.0)):
        sw = sorted([r for r in allr if r['na'] == na and r['tox'] == tox and r['L'] == L and abs(r['hw'] - 0.5) < 1e-9
                     and abs(r['shift']) < 1e-12 and r['level'] == 0], key=lambda z: z['psi_hat'])
        a150 = min(sw, key=lambda z: abs(z['alpha'] - 150))
        rng = [z for z in sw if z['psi_hat'] <= a150['psi_hat'] + 1e-9]
        mx = max(z['dL'] for z in sw); d12 = min(sw, key=lambda z: abs(z['psi_hat'] - 12))
        hsat = all(abs(z['dL'] / a150['dL'] - 1) <= 0.05 for z in rng)
        hmax = d12['dL'] <= 0.92 * mx
        p4[f'{na:.0e}'] = dict(psi_hat=[z['psi_hat'] for z in sw], alpha=[z['alpha'] for z in sw], dL=[z['dL'] for z in sw],
                               dL_over_Wd=[z['dL'] / z['Wd_nm'] for z in sw], lam_fit=[z.get('lam_fit_nm') for z in sw],
                               x_s_FS93=[z['x_s_FS93_nm'] for z in sw], dL_alpha150=a150['dL'], dL_max=mx, dL_psi12=d12['dL'],
                               ratio_psi12_to_max=d12['dL'] / mx,
                               verdict='H-sat SUPPORTED' if hsat and not hmax else 'H-max SUPPORTED' if hmax and not hsat else 'INCONCLUSIVE')
    res['P10_4'] = p4
    # P10-5 identity on new states with alpha >= 100
    p5 = [dict(file=r['file'], alpha=r['alpha'], err=r.get('identity_rel_err')) for r in new if r['alpha'] >= 100]
    res['P10_5'] = dict(rows=p5, max_abs_err=max(abs(x['err']) for x in p5 if x['err'] is not None),
                        verdict='SUPPORTED' if all(x['err'] is not None and abs(x['err']) <= 0.05 for x in p5) else 'REJECTED')
    # P10-6 Fjeldly-Shur eq. (9)
    p6 = [dict(file=r['file'], ratio=r['ratio_FS93'], x_s=r['x_s_FS93_nm'], dL=r['dL']) for r in new]
    res['P10_6'] = dict(rows=p6, min_ratio=min(x['ratio'] for x in p6), max_ratio=max(x['ratio'] for x in p6),
                        verdict='H-FS SUPPORTED' if all(abs(x['ratio'] - 1) <= 0.2 for x in p6) else 'H-FS REJECTED')
    # mesh check
    mesh = []
    for na, vg, L in ((1e16, vg1, 2.0), (1e18, vg1, 1.0)):
        r1 = one(new, na=na, vg=vg, level=1, hw=0.5)
        r0 = one(allr, na=na, vg=vg, level=0, hw=0.5, shift=0.0, L=L)
        mesh.append(dict(na=na, vg=vg, dL_lev0=r0['dL'], dL_lev1=r1['dL'], rel=r1['dL'] / r0['dL'] - 1, hw_lev0=r0['hw_real_um'], hw_lev1=r1['hw_real_um']))
    res['mesh'] = mesh
    res['all_states'] = [{k: v for k, v in r.items() if k not in ('dev',)} for r in allr]
    json.dump(res, open(os.path.join(OUT, 'c10_result.json'), 'w'), indent=1)
    print('new states', len(new), 'unconverged', res['failed_or_unconverged'])
    print('P10-1', res['P10_1']['verdict'], [(round(x['vg'], 4), round(x['dLLD_2prime'], 4), round(x['dLLD_1'], 4), round(x['ratio'], 4)) for x in p1])
    print('P10-1 post-hoc n_i image', [(round(x['vg'], 4), round(x['dLLD_image'], 4), round(x['dLLD_1'], 4), round(x['ratio'], 4), x['bounds_ok']) for x in ph])
    print('E10-1 drain-side bound violations', len(res['E10_1']['drain_side_bound_violations']), 'source-side violations', len(res['E10_1']['src_side_bound_violations']))
    print('P10-2 device1(hw0.5) dL/LD', round(res['P10_2']['device1_hw05']['dLLD'], 4), [(x['tag'], round(x['dLLD'], 4), round(x['rel_to_dev1_hw05'], 4), round(x['beta'] or 0, 2)) for x in dec])
    print('P10-3', v3, [(x['key'], round(x['ratio'], 4), round(x['l0_ratio'], 4), round(x['lamfit_ratio'], 4), round(x['dL_over_Wd'], 3), round(x['dL_over_Wd_base'], 3)) for x in p3])
    for k, v in p4.items():
        print('P10-4', k, v['verdict'], 'psi12/max', round(v['ratio_psi12_to_max'], 4), [(round(p, 2), round(d, 2)) for p, d in zip(v['psi_hat'], v['dL'])])
    print('P10-5', res['P10_5']['verdict'], 'max |err|', round(res['P10_5']['max_abs_err'], 4))
    print('P10-6', res['P10_6']['verdict'], 'ratio range', round(res['P10_6']['min_ratio'], 3), round(res['P10_6']['max_ratio'], 3))
    print('mesh', [(f"{m['na']:.0e}", round(m['dL_lev0'], 3), round(m['dL_lev1'], 3), round(m['rel'], 5)) for m in mesh])


if __name__ == '__main__':
    main()
