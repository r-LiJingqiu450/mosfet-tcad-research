#!/usr/bin/env python3
"""Solver-agnostic scorer for the low-current MOSFET benchmark (WP1).

A submission gives, per case id, the zero-bias conductance G0 (S/cm per unit width) and/or the per-end
correction dL (nm), and declares whether it was computed on the identical mesh and discrete network as
the case reference (same_mesh). Scoring rules (see benchmark/README.md):

1. Sign check: G0 <= 0 fails (a zero-bias conductance is positive).
2. Bracket check, ONLY when same_mesh is true: the Rayleigh bracket of the reference network,
   G_rows (1 - 1e-12) <= G <= G_slice (1 + 1e-12). It bounds that discrete network on that mesh and is not
   a bound on the continuum value, so it is never applied to other meshes or codes.
3. Tolerance check, ONLY when the case uncertainty is established:
   |G0/G0_ref - 1| <= tol.G0_rel and/or |dL - dL_ref| <= tol.dL_nm.
4. If no check applies, the case is "not_scorable" (with the reason).

A case passes when every applicable check passes. Stress cases (finite_bias_dd) have no reference value
and are always "not_scorable".

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


def score_case(case, entry):
    res = dict(case=case['id'], tier=case['tier'], checks=[], status=None, reason='')
    if case['family'] not in ('zero_bias_network', 'well_depth'):
        res.update(status='not_scorable', reason='stress case without a reference value')
        return res
    obs = case['observables']
    G_ref, dL_ref = obs['G0_S_per_cm']['value'], obs['dL_nm']['value']
    G_gca = obs['G_GCA_S_per_cm']['value']
    L_nm = case['device']['geometry']['L_um'] * 1e3
    G = entry.get('G0_S_per_cm')
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
