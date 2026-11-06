#!/usr/bin/env python3
"""Evaluate the WP2 production runs against the frozen WP3 pre-registration (sections 7-8; amendments A1-A7).

Applies the registered criteria exactly; nothing here is tuned to the results:
  H1   reference quality (A-R, B1): cross-path <= 5e-3, refinement <= 5e-3, forms / KCL <= 1e-9 at every level,
       B0 companion inside its own-mesh bracket at every level; kill if fewer than 4 of 6 pass.
  H-B  body effect bands (level 2), supported if >= 5 of 6 in band for each Delta (never a kill).
  H2   held-out robustness (level 1): forms <= 1e-9, KCL <= 1e-9, B0 bracket; any failure kills the release.
  P-DD default DEVSIM drain contact current at +/-1 mV vs the P1 reference (same level), threshold 1e-2.
  P-SS coupled small-signal G_E vs G_DD for the dev1 states (B1, level 1), threshold 1e-2.
  E-CIRC (informed expectation, A6): Richardson circuit-current G within 1e-6 of the same-level network G_DD.
Also computes the section-7 reference values and tolerances. Writes benchmark/wp2/results/evaluation.json.
"""
from __future__ import annotations

import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from run_queue import AR, HO, OUT, LEDGER  # noqa: E402

PSI_HAT = {'AR-1': 12, 'AR-2': 16, 'AR-3': 21, 'AR-4': 12, 'AR-5': 16, 'AR-6': 21,
           'HO-1': 15, 'HO-2': 18, 'HO-3': 15, 'HO-4': 18, 'HO-5': 12}
BANDS = {  # pre-registration section 8, H-B: (Delta_DD band, Delta_SD band)
    'AR-1': ((1.2e-3, 1.2e-2), (-7.2e-3, -7.8e-4)), 'AR-2': ((2.6e-5, 2.6e-4), (-1.6e-4, -1.7e-5)),
    'AR-3': ((2.0e-7, 2.0e-6), (-1.2e-6, -1.3e-7)), 'AR-4': ((2.4e-4, 6.1e-3), (-1.1e-3, -1.0e-4)),
    'AR-5': ((5.0e-6, 1.3e-4), (-2.3e-5, -2.2e-6)), 'AR-6': ((3.8e-8, 9.7e-7), (-1.7e-7, -1.7e-8))}
PDD = dict(**{s: 'fail' for s in ('AR-1', 'AR-2', 'AR-3', 'AR-4', 'HO-1', 'HO-2', 'HO-3', 'HO-5')},
           **{'AR-6': 'pass', 'AR-5': 'none', 'HO-4': 'none'})
PSS = ('AR-1', 'AR-2', 'AR-3', 'HO-1', 'HO-2', 'HO-5')
TOL_CROSS, TOL_REF, TOL_FORMS = 5e-3, 5e-3, 1e-9


def load(name):
    p = os.path.join(OUT, name)
    return json.load(open(p)) if os.path.exists(p) else None


def p1(s, b, lev):
    return load(f'P1_{s}_{b}_lev{lev}.json')


def p2(s, b, lev):
    return load(f'P2_{s}_{b}_lev{lev}.json')


def p4(s, mesh='refined'):
    return load(f'P4_{s}_{mesh}.json')


def rel(a, b):
    return a / b - 1.0


def order(q0, q1, q2):
    d1, d2 = q1 - q0, q2 - q1
    if d1 == 0 or d2 == 0 or (d1 > 0) != (d2 > 0):
        return None
    return math.log(abs(d1) / abs(d2)) / math.log(math.sqrt(2.0))


def h1_and_reference():
    out = {}
    for sid, dev, L, vg in AR:
        r = dict(state=sid, missing=[])
        runs = {(b, lv): p1(sid, b, lv) for b in ('B0', 'B1') for lv in (0, 1, 2)}
        P4r, P4b = p4(sid, 'refined'), p4(sid, 'base')
        r['missing'] = [f'P1_{b}_lev{lv}' for (b, lv), v in runs.items() if v is None] + \
                       (['P4_refined'] if P4r is None else []) + (['P4_base'] if P4b is None else [])
        if r['missing']:
            out[sid] = r
            continue
        n = {k: v['network'] for k, v in runs.items()}
        r['psi_hat_P1_lev2'] = runs[('B1', 2)]['onecol']['psi_hat']
        checks, ref = {}, {}
        for q in ('G_DD', 'G_SD'):
            L0, L1, L2 = n[('B1', 0)][q], n[('B1', 1)][q], n[('B1', 2)][q]
            Q4 = P4r['results']['B1'][q]
            c_cross, c_ref = rel(Q4, L2), rel(L2, L1)
            checks[q] = dict(cross_P4_vs_L2=c_cross, pass_i=abs(c_cross) <= TOL_CROSS,
                             refine_L2_vs_L1=c_ref, pass_ii=abs(c_ref) <= TOL_REF)
            u_disc = max(abs(c_cross), abs(c_ref))
            u_round = abs(n[('B1', 2)][f'{q}_rel_forms'])
            u_solver = n[('B1', 2)]['residual_drain_excited']
            tol = 2.0 * (u_disc + u_round + u_solver)
            ref[q] = dict(value=L2, levels=[L0, L1, L2], P4_refined=Q4, P4_base=P4b['results']['B1'][q],
                          P4_self_convergence=rel(Q4, P4b['results']['B1'][q]), observed_order=order(L0, L1, L2),
                          u_disc=u_disc, u_round=u_round, u_solver=u_solver, tol_rel=tol)
        # G_DB and the B0 companion: same rule, reported (not part of H1)
        for label, body, q, p4q in (('G_DB', 'B1', 'G_DB', 'G_DB'), ('G0_B0', 'B0', 'G0', 'G_SD')):
            L0, L1, L2 = n[(body, 0)][q], n[(body, 1)][q], n[(body, 2)][q]
            Q4 = P4r['results'][body][p4q]
            u_disc = max(abs(rel(Q4, L2)), abs(rel(L2, L1)))
            fr = n[(body, 2)].get('G_DB_rel_forms' if label == 'G_DB' else 'G_SD_rel_forms') or 0.0
            ref[label] = dict(value=L2, levels=[L0, L1, L2], P4_refined=Q4, P4_base=P4b['results'][body][p4q],
                              cross_P4_vs_L2=rel(Q4, L2), refine_L2_vs_L1=rel(L2, L1), observed_order=order(L0, L1, L2),
                              u_disc=u_disc, u_round=abs(fr), u_solver=n[(body, 2)]['residual_drain_excited'],
                              tol_rel=2.0 * (u_disc + abs(fr) + n[(body, 2)]['residual_drain_excited']),
                              note='same section-7 rule, reported; not part of H1')
        iii = []
        for lv in (0, 1, 2):
            b1, b0 = n[('B1', lv)], n[('B0', lv)]
            iii.append(dict(level=lv, G_DD_forms=b1['G_DD_rel_forms'], G_SD_forms=b1['G_SD_rel_forms'],
                            kcl=b1['kcl_closure_rel'], B0_bracket_ok=b0['bounds_ok'],
                            ok=abs(b1['G_DD_rel_forms']) <= TOL_FORMS and abs(b1['G_SD_rel_forms']) <= TOL_FORMS
                            and b1['kcl_closure_rel'] <= TOL_FORMS and bool(b0['bounds_ok'])))
        r['checks'] = checks
        r['self_consistency'] = iii
        r['pass'] = all(c['pass_i'] and c['pass_ii'] for c in checks.values()) and all(x['ok'] for x in iii)
        G_GCA = runs[('B1', 2)]['network']['G_GCA_1COL']
        L_nm = float(L) * 1e3
        r['reference'] = ref
        r['dL_SD_nm'] = 0.5 * L_nm * (1 - G_GCA / ref['G_SD']['value'])
        r['G_GCA_1COL'] = G_GCA
        dG = ref['G_SD']['tol_rel']
        r['dL_tol_nm'] = 0.5 * (L_nm - 2 * r['dL_SD_nm']) * dG          # d ln G = 2 d(dL)/(L - 2 dL)
        out[sid] = r
    passed = [s for s, r in out.items() if r.get('pass')]
    return out, dict(passed=passed, n_pass=len(passed), kill=len(passed) < 4)


def hb():
    out = {}
    for sid, (bdd, bsd) in BANDS.items():
        b1, b0 = p1(sid, 'B1', 2), p1(sid, 'B0', 2)
        if b1 is None or b0 is None:
            out[sid] = dict(missing=True); continue
        G0 = b0['network']['G0']
        ddd, dsd = rel(b1['network']['G_DD'], G0), rel(b1['network']['G_SD'], G0)
        P = p4(sid)
        p4d = None
        if P:
            p4d = dict(Delta_DD=rel(P['results']['B1']['G_DD'], P['results']['B0']['G_SD']),
                       Delta_SD=rel(P['results']['B1']['G_SD'], P['results']['B0']['G_SD']))
        out[sid] = dict(Delta_DD=ddd, band_DD=bdd, in_DD=bdd[0] <= ddd <= bdd[1],
                        Delta_SD=dsd, band_SD=bsd, in_SD=bsd[0] <= dsd <= bsd[1], P4=p4d,
                        G_DB_over_G0=b1['network']['G_DB'] / G0)
    done = [r for r in out.values() if not r.get('missing')]
    n_dd, n_sd = sum(r['in_DD'] for r in done), sum(r['in_SD'] for r in done)
    return out, dict(in_band_DD=n_dd, in_band_SD=n_sd, n=len(done), supported=n_dd >= 5 and n_sd >= 5)


def h2():
    out = {}
    for sid, dev, L, vg in HO:
        b1, b0 = p1(sid, 'B1', 1), p1(sid, 'B0', 1)
        if b1 is None or b0 is None:
            out[sid] = dict(missing=True); continue
        n1 = b1['network']
        a_, b_, c_ = abs(n1['G_DD_rel_forms']), n1['kcl_closure_rel'], bool(b0['network']['bounds_ok'])
        out[sid] = dict(psi_hat=b1['onecol']['psi_hat'], G_DD=n1['G_DD'], G_SD=n1['G_SD'], G_DB=n1['G_DB'],
                        a_forms=a_, b_kcl=b_, c_B0_bracket=c_, G_SD_forms=n1['G_SD_rel_forms'],
                        **{'pass': a_ <= TOL_FORMS and b_ <= TOL_FORMS and c_})
    done = [r for r in out.values() if not r.get('missing')]
    return out, dict(n=len(done), n_pass=sum(r['pass'] for r in done),
                     kill=any(not r['pass'] for r in done), complete=len(done) == len(HO))


def p2_eval():
    rows = []
    for sid, dev, L, vg in AR + HO:
        for b, lv in (('B1', 1), ('B1', 2), ('B0', 1)):
            d = p2(sid, b, lv)
            if d is None:
                continue
            ref_run = p1(sid, b, lv)
            ref = ref_run['network']['G_DD'] if ref_run else None
            s = d['symmetric']
            row = dict(state=sid, body=b, level=lv, prediction=PDD[sid], G_ref_P1=ref,
                       G_ref_same_run=d['network_same_run']['G_DD'],
                       same_run_vs_P1=rel(d['network_same_run']['G_DD'], ref) if ref else None,
                       G_times_1mV=(ref or d['network_same_run']['G_DD']) * 1e-3)
            for k, lab in (('I_drain', 'default'), ('I_circuit', 'circuit'), ('I_weighted', 'weighted')):
                g1, g01 = s['1mV'][k]['G'], s['0.1mV'][k]['G']
                base = ref or d['network_same_run']['G_DD']
                row[f'{lab}_1mV'] = rel(g1, base)
                row[f'{lab}_0.1mV'] = rel(g01, base)
                row[f'{lab}_richardson'] = rel((100 * g01 - g1) / 99, base)
            row['source_minus_1mV'] = rel(s['1mV']['I_source_minus'], ref or d['network_same_run']['G_DD'])
            row['newton_iterations'] = [h['iterations'] for st in d['states'] for h in st['continuation']]
            row['indicator_check'] = d['indicator_check_max_rel']
            if PDD[sid] == 'fail':
                row['pdd_ok'] = abs(row['default_1mV']) > 1e-2
            elif PDD[sid] == 'pass':
                row['pdd_ok'] = abs(row['default_1mV']) <= 1e-2
            else:
                row['pdd_ok'] = None
            row['ecirc_ok'] = abs(row['circuit_richardson']) <= 1e-6
            rows.append(row)
    per_state = {}
    for r in rows:
        per_state.setdefault(r['state'], []).append(r)
    summ = {}
    for sid, rs in per_state.items():
        oks = [r['pdd_ok'] for r in rs if r['pdd_ok'] is not None]
        summ[sid] = dict(prediction=PDD[sid], runs=len(rs),
                         outcome=('no prediction' if PDD[sid] == 'none' else ('confirmed' if all(oks) else 'refuted')),
                         default_1mV=[r['default_1mV'] for r in rs])
    return rows, summ, dict(ecirc_all=all(r['ecirc_ok'] for r in rows) if rows else None,
                            ecirc_worst=max((abs(r['circuit_richardson']) for r in rows), default=None))


def pss_eval():
    rows = []
    for sid, dev, L, vg in AR + HO:
        for b in ('B1', 'B0'):
            d = p1(sid, b, 1)
            if d is None or 'smallsignal' not in d:
                continue
            ss = d['smallsignal']
            rows.append(dict(state=sid, body=b, device=dev, predicted=(sid in PSS and b == 'B1'),
                             G_E_rel=ss['rel_G_E_vs_G_DD'], G_M_rel=ss['rel_G_M_vs_G_DD'], kappa_G=ss['kappa_G'],
                             E_dwr=ss['E_dwr'], E_round=ss['E_round'], eta=ss['eta'],
                             ok=(abs(ss['rel_G_E_vs_G_DD']) > 1e-2) if (sid in PSS and b == 'B1') else None))
    pred = [r for r in rows if r['predicted']]
    return rows, dict(n_predicted=len(pred), n_confirmed=sum(bool(r['ok']) for r in pred))


def psi_hat_check():
    out = {}
    for sid, dev, L, vg in AR + HO:
        d = p1(sid, 'B1', 1)
        if d:
            out[sid] = dict(expected=PSI_HAT[sid], P1_lev1=d['onecol']['psi_hat'])
    return out


def ledger_summary():
    if not os.path.exists(LEDGER):
        return {}
    rows = [json.loads(x) for x in open(LEDGER)]
    return dict(runs=len(rows), failures=[r['job'] for r in rows if r['exit'] != 0],
                total_wall_h=sum(r['wall_s'] for r in rows) / 3600, max_rss_MB=max(r['maxrss_MB'] for r in rows))


def main():
    H1, H1s = h1_and_reference()
    HB, HBs = hb()
    H2, H2s = h2()
    P2rows, PDDs, Es = p2_eval()
    PSSr, PSSs = pss_eval()
    res = dict(H1=H1, H1_summary=H1s, HB=HB, HB_summary=HBs, H2=H2, H2_summary=H2s, P2_rows=P2rows,
               PDD_by_state=PDDs, ECIRC=Es, PSS_rows=PSSr, PSS_summary=PSSs, psi_hat=psi_hat_check(),
               ledger=ledger_summary())
    with open(os.path.join(HERE, 'results', 'evaluation.json'), 'w') as f:
        json.dump(res, f, indent=1, default=lambda o: bool(o) if isinstance(o, (bool,)) else str(o))
    print('== H1 (A-R reference quality) ==')
    for sid, r in H1.items():
        if r.get('missing'):
            print(f'  {sid}: missing {r["missing"]}'); continue
        c = r['checks']
        print(f"  {sid} psi_hat={r['psi_hat_P1_lev2']:.2f} PASS={r['pass']} | G_DD cross {c['G_DD']['cross_P4_vs_L2']:+.2e} refine "
              f"{c['G_DD']['refine_L2_vs_L1']:+.2e} | G_SD cross {c['G_SD']['cross_P4_vs_L2']:+.2e} refine {c['G_SD']['refine_L2_vs_L1']:+.2e} "
              f"| forms/KCL/bracket ok per level {[x['ok'] for x in r['self_consistency']]} | tol(G_SD) {r['reference']['G_SD']['tol_rel']:.2e}")
    print('  summary', H1s)
    print('== H-B (body effect, level 2) ==')
    for sid, r in HB.items():
        if not r.get('missing'):
            print(f"  {sid} Delta_DD={r['Delta_DD']:+.3e} in {r['band_DD']} -> {r['in_DD']} | Delta_SD={r['Delta_SD']:+.3e} in {r['band_SD']} -> {r['in_SD']}")
    print('  summary', HBs)
    print('== H2 (held-out robustness) ==')
    for sid, r in H2.items():
        if not r.get('missing'):
            print(f"  {sid} psi_hat={r['psi_hat']:.2f} forms={r['a_forms']:.1e} kcl={r['b_kcl']:.1e} B0bracket={r['c_B0_bracket']} PASS={r['pass']}")
    print('  summary', H2s)
    print('== P-DD / E-CIRC (P2 runs) ==')
    for r in P2rows:
        print(f"  {r['state']} {r['body']} lev{r['level']} GV1mV={r['G_times_1mV']:.1e} pred={r['prediction']:5s} default1mV={r['default_1mV']:+.3e} "
              f"ok={r['pdd_ok']} | circuit R={r['circuit_richardson']:+.1e} weighted R={r['weighted_richardson']:+.1e} | same-run vs P1 {r['same_run_vs_P1']:+.1e}")
    print('  by state', {k: v['outcome'] for k, v in PDDs.items()}, '| E-CIRC', Es)
    print('== P-SS (coupled small-signal, level 1) ==')
    for r in PSSr:
        print(f"  {r['state']} {r['body']} G_E {r['G_E_rel']:+.3e} G_M {r['G_M_rel']:+.3e} kappa {r['kappa_G']:.1e} predicted={r['predicted']} ok={r['ok']}")
    print('  summary', PSSs)
    print('psi_hat', {k: round(v['P1_lev1'], 3) for k, v in res['psi_hat'].items()})
    print('ledger', res['ledger'])
    return 0


if __name__ == '__main__':
    sys.exit(main())
