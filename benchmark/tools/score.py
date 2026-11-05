#!/usr/bin/env python3
"""Solver-agnostic scorer for the low-current MOSFET benchmark (WP1).

A submission gives, per case id, the zero-bias conductance G0 (S/cm per unit width) and/or the per-end
correction dL (nm), and declares whether it was computed on the identical mesh and discrete network as
the case reference (same_mesh). Scoring rules (see benchmark/README.md):

Two-terminal cases (schema 1.0 and B0 companions of schema 1.1):
1. Sign check: G0 <= 0 fails (a zero-bias conductance is positive).
2. Bracket check, ONLY when same_mesh is true: the Rayleigh bracket of the reference network,
   G_rows (1 - 1e-12) <= G <= G_slice (1 + 1e-12). It bounds that discrete network on that mesh and is not
   a bound on the continuum value, so it is never applied to other meshes or codes.
3. Tolerance check, ONLY when the case uncertainty is established:
   |G0/G0_ref - 1| <= tol.G0_rel and/or |dL - dL_ref| <= tol.dL_nm.
4. If no check applies, the case is "not_scorable" (with the reason).

A case passes when every applicable check passes. Stress cases (finite_bias_dd) have no reference value
and are always "not_scorable".

Three-terminal cases (schema 1.1, family zero_bias_3terminal, ohmic body B1): submissions give G_DD_S_per_cm,
G_SD_S_per_cm and/or G_DB_S_per_cm (and optionally dL_nm). Same-mesh submissions are checked within tau_num = 1e-6
relative (no bracket exists with an ohmic body) against the reference-path value on the declared mesh (mesh_id, one of
observables.same_mesh_references; default the case mesh); others against the established tolerance. G_DD and G_SD
(and dL) decide pass/fail; G_DB is compared and reported only (no pre-registered tolerance).
B0 companions (tier B-R) are two-terminal cases: G_DD or G_SD submissions are read as G0 (two_terminal_aliases).

Usage:  python benchmark/tools/score.py SUBMISSION.json [--cases benchmark/cases] [--out REPORT.json]
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BENCH = os.path.join(ROOT, 'benchmark')
BRACKET_RTOL = 1e-12


def load_cases(cases_dir=None):
    d = cases_dir or os.path.join(BENCH, 'cases')
    out = {}
    for p in sorted(glob.glob(os.path.join(d, '*.json'))):
        with open(p, encoding='utf-8') as f:
            c = json.load(f)
        out[c['id']] = c
    return out


def validate_submission(sub):
    try:
        import jsonschema
    except ImportError:            # scoring still works; validation is reported as skipped
        return 'skipped (jsonschema not installed)'
    with open(os.path.join(BENCH, 'schema', 'submission.schema.json'), encoding='utf-8') as f:
        jsonschema.validate(sub, json.load(f))
    return 'ok'


THREE_T = ('G_DD', 'G_SD', 'G_DB')


def score_case_3t(case, entry, res):
    """Three-terminal zero-bias case (schema 1.1, body B1): per-observable checks.

    same_mesh: |Q/Q_ref - 1| <= tau_num (pre-registered 1e-6; there is no Rayleigh bracket with an ohmic body), where Q_ref is
               the reference-path value on the declared mesh (entry mesh_id, a key of observables.same_mesh_references;
               default the case mesh). Scored observables: observables.scored (G_DD, G_SD); G_DB is compared and reported
               but does not decide pass/fail (no pre-registered tolerance).
    tolerance: |Q/Q_ref - 1| <= tolerance[Q_rel] and |dL - dL_ref| <= tolerance[dL_nm] where established."""
    obs, unc = case['observables'], case['uncertainty']
    scored = {s[:-len('_S_per_cm')] for s in obs.get('scored', []) if s.endswith('_S_per_cm')}
    vals = {q: entry.get(f'{q}_S_per_cm') for q in THREE_T}
    dL = entry.get('dL_nm')
    if all(v is None for v in vals.values()) and dL is None:
        res.update(status='not_scorable', reason='three-terminal case: submit G_DD_S_per_cm, G_SD_S_per_cm, G_DB_S_per_cm and/or dL_nm')
        return res
    for q, v in vals.items():
        if v is not None and not v > 0:
            res['checks'].append(dict(check=f'sign_{q}', passed=False, detail=f'{q} = {v!r} is not positive'))
    if entry.get('same_mesh'):
        tau = unc['tau_num_rel']
        mid = entry.get('mesh_id', case['mesh'].get('mesh_id'))
        smr = obs.get('same_mesh_references', {})
        if mid not in smr:
            res['checks'].append(dict(check='same_mesh_id', passed=False, detail=f'mesh_id {mid!r} is not one of {sorted(smr)}'))
        else:
            for q, v in vals.items():
                if v is not None and v > 0:
                    r = v / smr[mid][f'{q}_S_per_cm'] - 1.0
                    c = dict(check=f'same_mesh_tau_num_{q}', passed=abs(r) <= tau, rel_error=r, tol=tau, mesh_id=mid)
                    if q not in scored:
                        c['informational'] = True
                    res['checks'].append(c)
    tol = unc['tolerance']
    # A same-mesh submission on a coarser published mesh is scored by tau_num on that mesh; its distance from the (finer)
    # case reference is then discretization error of the declared mesh and is reported, not scored.
    coarse_same_mesh = bool(entry.get('same_mesh')) and entry.get('mesh_id', case['mesh'].get('mesh_id')) != case['mesh'].get('mesh_id')
    if tol != 'not_established':
        for q, v in vals.items():
            if v is not None and v > 0:
                r = v / obs[f'{q}_S_per_cm']['value'] - 1.0
                if f'{q}_rel' in tol:
                    c = dict(check=f'tolerance_{q}', passed=abs(r) <= tol[f'{q}_rel'], rel_error=r, tol=tol[f'{q}_rel'])
                    if coarse_same_mesh:
                        c.update(informational=True, detail='same-mesh submission on a coarser published mesh: reported, not scored')
                    res['checks'].append(c)
                else:
                    res['checks'].append(dict(check=f'reported_{q}', passed=True, rel_error=r, informational=True,
                                              detail='no pre-registered tolerance for this observable'))
        if dL is not None:
            e = dL - obs['dL_nm']['value']
            c = dict(check='tolerance_dL', passed=abs(e) <= tol['dL_nm'], error_nm=e, tol_nm=tol['dL_nm'])
            if coarse_same_mesh:
                c.update(informational=True, detail='same-mesh submission on a coarser published mesh: reported, not scored')
            res['checks'].append(c)
    if not [c for c in res['checks'] if not c.get('informational')]:
        res.update(status='not_scorable',
                   reason='no established tolerance for this case and the submission is not on a reference mesh')
        return res
    failed = [c['check'] for c in res['checks'] if not c['passed'] and not c.get('informational')]
    res['status'] = 'fail' if failed else 'pass'
    res['reason'] = ('failed: ' + ', '.join(failed)) if failed else 'all applicable checks passed'
    return res


def score_case(case, entry):
    res = dict(case=case['id'], tier=case['tier'], checks=[], status=None, reason='')
    if case['family'] == 'zero_bias_3terminal':
        return score_case_3t(case, entry, res)
    if case['family'] not in ('zero_bias_network', 'well_depth'):
        res.update(status='not_scorable', reason='stress case without a reference value')
        return res
    obs = case['observables']
    G_ref, dL_ref = obs['G0_S_per_cm']['value'], obs['dL_nm']['value']
    G_gca = obs['G_GCA_S_per_cm']['value']
    L_nm = case['device']['geometry']['L_um'] * 1e3
    G = entry.get('G0_S_per_cm')
    if G is None and case['observables'].get('two_terminal_aliases'):
        G = entry.get('G_SD_S_per_cm', entry.get('G_DD_S_per_cm'))   # B0 companion: G_DD = G_SD = G0
    dL = entry.get('dL_nm')
    if G is not None and not G > 0:
        res['checks'].append(dict(check='sign', passed=False, detail=f'G0 = {G!r} is not positive'))
    # G equivalent of a dL-only submission (same G_GCA reference by construction of dL)
    G_eq = G if G is not None else (G_gca / (1.0 - 2.0 * dL / L_nm) if dL is not None and 2.0 * dL < L_nm else None)
    if entry.get('same_mesh'):
        b = obs['bracket']
        lo, hi = b['G_rows_S_per_cm'] * (1 - BRACKET_RTOL), b['G_slice_S_per_cm'] * (1 + BRACKET_RTOL)
        ok = G_eq is not None and lo <= G_eq <= hi
        res['checks'].append(dict(check='same_mesh_bracket', passed=bool(ok), value=G_eq, lower=lo, upper=hi))
    tol = case['uncertainty']['tolerance']
    if tol != 'not_established':
        if G is not None:
            r = G / G_ref - 1.0
            res['checks'].append(dict(check='tolerance_G0', passed=abs(r) <= tol['G0_rel'], rel_error=r, tol=tol['G0_rel']))
        if dL is not None:
            e = dL - dL_ref
            res['checks'].append(dict(check='tolerance_dL', passed=abs(e) <= tol['dL_nm'], error_nm=e, tol_nm=tol['dL_nm']))
    if not res['checks']:
        res.update(status='not_scorable',
                   reason='no established tolerance for this case and the submission is not on the reference mesh (bracket not applicable)')
        return res
    failed = [c['check'] for c in res['checks'] if not c['passed']]
    res['status'] = 'fail' if failed else 'pass'
    res['reason'] = ('failed: ' + ', '.join(failed)) if failed else 'all applicable checks passed'
    return res


def score(sub, cases=None):
    cases = cases if cases is not None else load_cases()
    out = dict(submission=sub.get('submission', {}), validation=validate_submission(sub), results=[], summary={})
    for cid, entry in sorted(sub['results'].items()):
        if cid not in cases:
            out['results'].append(dict(case=cid, status='unknown_case', checks=[], reason='no such case id'))
            continue
        out['results'].append(score_case(cases[cid], entry))
    s = {}
    for r in out['results']:
        s[r['status']] = s.get(r['status'], 0) + 1
    out['summary'] = dict(sorted(s.items()))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('submission')
    ap.add_argument('--cases', default=None)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    with open(a.submission, encoding='utf-8') as f:
        sub = json.load(f)
    rep = score(sub, load_cases(a.cases))
    txt = json.dumps(rep, indent=1, sort_keys=True, allow_nan=False) + '\n'
    if a.out:
        with open(a.out, 'w', encoding='utf-8', newline='\n') as f:
            f.write(txt)
    for r in rep['results']:
        print(f"{r['status']:13s} {r['case']}  {r['reason']}")
    print('summary', json.dumps(rep['summary']))
    return 0


if __name__ == '__main__':
    sys.exit(main())
