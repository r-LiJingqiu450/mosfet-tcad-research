"""WP1 acceptance tests for the low-current benchmark package.

Expected counts are fixed in advance by the WP1 specification
(claude/WP1_Codex_Spec_Benchmark_Freeze_2026-11-04.md); they must not be tuned to match.
"""
from __future__ import annotations

import glob
import hashlib
import json
import os
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BENCH = os.path.join(ROOT, 'benchmark')
sys.path.insert(0, os.path.join(BENCH, 'tools'))

import score as S  # noqa: E402

CASES = S.load_cases()
AW = sorted(k for k, c in CASES.items() if c['tier'] == 'A-W')
BW = sorted(k for k, c in CASES.items() if c['tier'] == 'B-W')


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def sub(results, name='test'):
    return dict(submission=dict(name=name), results=results)


def statuses(rep):
    return {r['case']: r['status'] for r in rep['results']}


# ------------------------------------------------------------------ package integrity
def test_tier_counts():
    assert len(AW) == 3
    assert len(BW) == 28


def test_all_cases_validate_against_schema():
    jsonschema = pytest.importorskip('jsonschema')
    schema = json.load(open(os.path.join(BENCH, 'schema', 'case.schema.json')))
    for c in CASES.values():
        jsonschema.validate(c, schema)


def test_manifest_hashes_match():
    lines = open(os.path.join(BENCH, 'MANIFEST.sha256')).read().splitlines()
    assert len(lines) == len(CASES) + 1 + 2       # cases + INDEX.json + 2 schemas
    for ln in lines:
        h, p = ln.split('  ', 1)
        assert sha(os.path.join(ROOT, p)) == h, p


def test_sources_trace_to_base_commit_and_hash():
    for c in CASES.values():
        base = c['provenance']['base_commit']
        for s in c['provenance']['sources']:
            assert sha(os.path.join(ROOT, s['path'])) == s['sha256'], s['path']
            at_base = subprocess.run(['git', '-C', ROOT, 'rev-parse', f"{base}:{s['path']}"], capture_output=True, text=True, check=True).stdout.strip()
            here = subprocess.run(['git', '-C', ROOT, 'hash-object', s['path']], capture_output=True, text=True, check=True).stdout.strip()
            assert at_base == here, s['path']


def test_export_is_deterministic():
    files = sorted(glob.glob(os.path.join(BENCH, 'cases', '*.json'))) + [os.path.join(BENCH, 'INDEX.json'), os.path.join(BENCH, 'MANIFEST.sha256')]
    before = {p: sha(p) for p in files}
    r = subprocess.run([sys.executable, os.path.join(BENCH, 'tools', 'export_cases.py')], cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    after_files = sorted(glob.glob(os.path.join(BENCH, 'cases', '*.json'))) + [os.path.join(BENCH, 'INDEX.json'), os.path.join(BENCH, 'MANIFEST.sha256')]
    assert after_files == files
    assert {p: sha(p) for p in files} == before


# ------------------------------------------------------------------ retrospective discrimination (pre-registered counts)
def drain_value(c):
    for e in c['known_bad_outputs'] + c['cross_paths']:
        if e['name'] == 'drain_side_flux_sum':
            return e['G0_S_per_cm']
    raise KeyError(c['id'])


def test_bw_drain_side_flagged_16_of_28():
    rep = S.score(sub({k: dict(G0_S_per_cm=drain_value(CASES[k]), same_mesh=True) for k in BW}), CASES)
    st = statuses(rep)
    assert sum(v == 'fail' for v in st.values()) == 16
    # case-by-case agreement with the independent C10 flags (bounds_ok_drain)
    c10 = {os.path.normpath(r['file']): r['bounds_ok_drain'] for r in json.load(open(os.path.join(ROOT, 'results/phase3/C10/analysis/c10_result.json')))['all_states']}
    for k in BW:
        f = [s['path'] for s in CASES[k]['provenance']['sources'] if 'runs/' in s['path'] or 'subth/' in s['path']][0]
        assert (st[k] == 'fail') == (not c10[os.path.normpath(f)]), k


def test_bw_source_side_flagged_0_of_28():
    rep = S.score(sub({k: dict(G0_S_per_cm=CASES[k]['observables']['G0_S_per_cm']['value'], same_mesh=True) for k in BW}), CASES)
    assert all(v == 'pass' for v in statuses(rep).values())


@pytest.mark.parametrize('name', ['coupled_small_signal_G_E', 'coupled_small_signal_G_M'])
def test_aw_coupled_small_signal_flagged_3_of_3(name):
    res = {}
    for k in AW:
        v = [e for e in CASES[k]['known_bad_outputs'] if e['name'] == name][0]['G0_S_per_cm']
        res[k] = dict(G0_S_per_cm=v, same_mesh=True)
    st = statuses(S.score(sub(res), CASES))
    assert sum(v == 'fail' for v in st.values()) == 3


def test_aw_edge_energy_inside_own_bracket_3_of_3():
    res = {k: dict(G0_S_per_cm=[e for e in CASES[k]['cross_paths'] if e['name'] == 'edge_energy'][0]['G0_S_per_cm'], same_mesh=True) for k in AW}
    rep = S.score(sub(res), CASES)
    for r in rep['results']:
        br = [c for c in r['checks'] if c['check'] == 'same_mesh_bracket'][0]
        assert br['passed'], r['case']
        assert r['status'] == 'pass', r


def test_aw_independent_fv_passes_tolerance_consistency_only():
    """Consistency check, not independent validation: the FV spread is part of the uncertainty estimate."""
    res = {k: dict(dL_nm=[e for e in CASES[k]['cross_paths'] if e['name'] == 'independent_finite_volume'][0]['dL_nm'], same_mesh=False) for k in AW}
    rep = S.score(sub(res), CASES)
    assert all(r['status'] == 'pass' for r in rep['results'])
    assert all([c['check'] for c in r['checks']] == ['tolerance_dL'] for r in rep['results'])


def test_bracket_not_applied_to_other_meshes():
    """The FV dL lies outside the DEVSIM-mesh bracket at 2 of 3 A-W states; scoring it as same_mesh would wrongly fail."""
    outside = 0
    for k in AW:
        c = CASES[k]
        dl = [e for e in c['cross_paths'] if e['name'] == 'independent_finite_volume'][0]['dL_nm']
        b = c['observables']['bracket']
        outside += not (b['dL_rows_nm'] <= dl <= b['dL_slice_nm'])
    assert outside == 2


def test_stress_cases_not_scorable_and_unknown_ids_reported():
    fb = sorted(k for k, c in CASES.items() if c['family'] == 'finite_bias_dd')
    rep = S.score(sub({fb[0]: dict(G0_S_per_cm=1.0, same_mesh=False), 'no_such_case': dict(G0_S_per_cm=1.0, same_mesh=False)}), CASES)
    st = statuses(rep)
    assert st[fb[0]] == 'not_scorable' and st['no_such_case'] == 'unknown_case'


def test_nonpositive_conductance_fails():
    k = AW[0]
    rep = S.score(sub({k: dict(G0_S_per_cm=-1e-12, same_mesh=False)}), CASES)
    assert rep['results'][0]['status'] == 'fail'


def test_bw_without_mesh_identity_is_not_scorable():
    k = BW[0]
    rep = S.score(sub({k: dict(G0_S_per_cm=CASES[k]['observables']['G0_S_per_cm']['value'], same_mesh=False)}), CASES)
    assert rep['results'][0]['status'] == 'not_scorable'
