"""WP4 tests: DEVSIM's user-facing output paths scored as external submissions (examples/), not stored in the case files.

They check traceability (submitted values recomputed from the WP4 production records), reproducibility of the stored score
reports, agreement of the scorer with the registered evaluation (benchmark/wp4/results/evaluation_wp4.json, prediction X2),
and that WP4 wrote nothing into the A-R case files (WP4 pre-registration section 6). They are package checks, not validation.
"""
from __future__ import annotations

import json
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BENCH = os.path.join(ROOT, 'benchmark')
sys.path.insert(0, os.path.join(BENCH, 'tools'))

import score as S  # noqa: E402

CASES = S.load_cases()
EX = os.path.join(BENCH, 'examples')
PROD = os.path.join(BENCH, 'wp4', 'results', 'production')
EV = json.load(open(os.path.join(BENCH, 'wp4', 'results', 'evaluation_wp4.json')))
AR = sorted(k for k, c in CASES.items() if c['tier'] == 'A-R')
SID = {k: CASES[k]['device']['state_id'] for k in AR}
NAMES = ('three_terminal_devsim_contact_currents_extended_level1_mesh', 'three_terminal_devsim_ac_small_signal_double_level1_mesh',
         'three_terminal_devsim_ac_small_signal_extended_level1_mesh')


def example(name):
    return (json.load(open(os.path.join(EX, f'{name}.submission.json'))), json.load(open(os.path.join(EX, f'{name}.report.json'))))


def level1_mesh(k):
    smr = CASES[k]['observables']['same_mesh_references']
    return [m for m, v in smr.items() if v['level'] == 1][0]


def rich(g1, g01):
    return (100.0 * g01 - g1) / 99.0


@pytest.mark.parametrize('name', NAMES)
def test_wp4_example_report_reproduces(name):
    sub, rep = example(name)
    again = S.score(sub, CASES)
    assert again['summary'] == rep['summary']
    assert again['results'] == rep['results']


def test_extended_contact_current_values_trace_to_records():
    sub, _ = example(NAMES[0])
    assert sorted(sub['results']) == AR
    for k in AR:
        d = json.load(open(os.path.join(PROD, f"P2X_{SID[k]}_B1_lev1.json")))
        s = d['symmetric']
        I_S = {st['V_D']: st['I_source'] for st in d['states']}
        r = sub['results'][k]
        assert r['G_DD_S_per_cm'] == rich(s['1mV']['I_drain']['G'], s['0.1mV']['I_drain']['G'])
        assert r['G_SD_S_per_cm'] == rich(-(I_S[1e-3] - I_S[-1e-3]) / 2e-3, -(I_S[1e-4] - I_S[-1e-4]) / 2e-4)
        assert r['same_mesh'] is True and r['mesh_id'] == level1_mesh(k)
        assert d['precision']['parameters_after_run'] == dict(extended_solver=True, extended_model=True, extended_equation=True)


@pytest.mark.parametrize('prec', ['double', 'extended'])
def test_ac_values_trace_to_records(prec):
    sub, _ = example(f'three_terminal_devsim_ac_small_signal_{prec}_level1_mesh')
    assert sorted(sub['results']) == AR
    for k in AR:
        d = json.load(open(os.path.join(PROD, f"AC_{SID[k]}_B1_lev1_{prec}.json")))
        assert d['precision'] == prec and d['frequency_Hz'] == 1e-6
        assert sub['results'][k]['G_DD_S_per_cm'] == d['ac']['G_AC'] == -d['ac']['I_real']
        assert sub['results'][k]['mesh_id'] == level1_mesh(k)


def test_scorer_agrees_with_registered_x2():
    _, rep = example(NAMES[0])
    x2 = {r['state']: r['X2'] for r in EV['B1'] if not r.get('missing')}
    assert len(x2) == 6
    for res in rep['results']:
        g_dd = [c for c in res['checks'] if c['check'] == 'same_mesh_tau_num_G_DD']
        assert len(g_dd) == 1 and g_dd[0]['passed'] == x2[SID[res['case']]]


def test_wp4_paths_not_written_into_case_files():
    for k in AR:
        names = [e['name'] for e in CASES[k]['known_bad_outputs'] + CASES[k]['cross_paths']]
        assert not any('extended' in n or n.startswith('devsim_ac_') for n in names)
        assert 'WP4' not in json.dumps(CASES[k]['provenance'])
