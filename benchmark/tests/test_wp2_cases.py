"""WP2 acceptance tests: resolved-well three-terminal tier (A-R), held-out tier (HO) and B0 companions (B-R).

These tests check package integrity and that the scorer applies the frozen pre-registration rules
(benchmark/specs/WP3_Preregistration_AR_and_HeldOut_2026-11-04.md) consistently with the exported cases.
Expected tier counts are taken from benchmark/wp2/results/evaluation.json (the registered evaluation), not tuned here.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BENCH = os.path.join(ROOT, 'benchmark')
sys.path.insert(0, os.path.join(BENCH, 'tools'))

import score as S  # noqa: E402

CASES = S.load_cases()
EV = json.load(open(os.path.join(BENCH, 'wp2', 'results', 'evaluation.json')))
AR = sorted(k for k, c in CASES.items() if k.startswith('AR__'))
HO = sorted(k for k, c in CASES.items() if c['tier'] == 'HO')
BR = sorted(k for k, c in CASES.items() if c['tier'] == 'B-R')
THREE = AR + HO


def sub(results):
    return dict(submission=dict(name='test'), results=results)


def one(case_id, entry):
    return S.score(sub({case_id: entry}), CASES)['results'][0]


def refs(c, mesh_id=None):
    m = c['observables']['same_mesh_references'][mesh_id or c['mesh']['mesh_id']]
    return {q: m[f'{q}_S_per_cm'] for q in ('G_DD', 'G_SD', 'G_DB')}


# ------------------------------------------------------------------ counts follow the registered evaluation
def test_counts_follow_evaluation():
    h1 = EV['H1_summary']
    assert len(AR) == 6 and len(HO) == 5 and len(BR) == 6 * 3 + 5
    n_ar = sum(CASES[k]['tier'] == 'A-R' for k in AR)
    assert n_ar == (0 if h1['kill'] else h1['n_pass'])
    for k in AR:
        sid = CASES[k]['device']['state_id']
        assert (CASES[k]['tier'] == 'A-R') == (bool(EV['H1'][sid]['pass']) and not h1['kill']), k
        assert (CASES[k]['uncertainty']['status'] == 'established') == (CASES[k]['tier'] == 'A-R'), k
    assert EV['H2_summary']['complete'] and not EV['H2_summary']['kill']


def test_mesh_files_match_cases_and_records():
    for k in THREE + BR:
        c = CASES[k]
        p = os.path.join(ROOT, c['mesh']['node_lines_file'])
        assert hashlib.sha256(open(p, 'rb').read()).hexdigest() == c['mesh']['node_lines_file_sha256'], k
        m = json.load(open(p))
        assert (m['nx'], m['ny'], m['level']) == (c['mesh']['nx'], c['mesh']['ny'], c['mesh']['level']), k
        assert m['cross_check_vs_P1_records'] and all(m['cross_check_vs_P1_records'].values()), k
        assert len(m['x_cm']) == m['nx'] and len(m['y_oxide_cm']) + len(m['y_si_cm']) - 1 == m['ny']


# ------------------------------------------------------------------ same-mesh scoring with tau_num (section 7)
def test_reference_values_pass_on_every_reference_mesh():
    for k in THREE:
        c = CASES[k]
        for mid in c['observables']['same_mesh_references']:
            r = one(k, dict(**{f'{q}_S_per_cm': v for q, v in refs(c, mid).items()}, same_mesh=True, mesh_id=mid))
            assert r['status'] == 'pass', (k, mid, r)


@pytest.mark.parametrize('factor,expected', [(1 + 5e-7, 'pass'), (1 - 5e-7, 'pass'), (1 + 2e-6, 'fail'), (1 - 2e-6, 'fail')])
def test_tau_num_discriminates(factor, expected):
    for k in THREE:
        c = CASES[k]
        r = one(k, dict(G_DD_S_per_cm=refs(c)['G_DD'] * factor, G_SD_S_per_cm=refs(c)['G_SD'], same_mesh=True))
        assert r['status'] == expected, (k, r)


def test_g_db_is_reported_not_scored():
    for k in THREE:
        c = CASES[k]
        r = one(k, dict(**{f'{q}_S_per_cm': v for q, v in refs(c).items()}, same_mesh=True) | dict(G_DB_S_per_cm=refs(c)['G_DB'] * 1.5))
        assert r['status'] == 'pass', r
        assert [x for x in r['checks'] if x['check'] == 'same_mesh_tau_num_G_DB'][0]['informational']


def test_unknown_mesh_id_fails_and_other_mesh_without_tolerance_is_not_scorable():
    k = HO[0]
    c = CASES[k]
    assert one(k, dict(G_DD_S_per_cm=refs(c)['G_DD'], same_mesh=True, mesh_id='no_such_mesh'))['status'] == 'fail'
    assert one(k, dict(G_DD_S_per_cm=refs(c)['G_DD'], same_mesh=False))['status'] == 'not_scorable'


# ------------------------------------------------------------------ tolerance scoring (A-R only)
def test_tolerance_rule_is_section_7():
    for k in AR:
        c = CASES[k]
        if c['uncertainty']['status'] != 'established':
            continue
        tol = c['uncertainty']['tolerance']
        for q in ('G_DD', 'G_SD'):
            comp = c['uncertainty']['components'][q]
            assert abs(tol[f'{q}_rel'] - 2 * (comp['u_disc'] + comp['u_round'] + comp['u_solver'])) <= 1e-15 * tol[f'{q}_rel']
            ref = c['observables'][f'{q}_S_per_cm']['value']
            assert one(k, {f'{q}_S_per_cm': ref * (1 + 0.99 * tol[f'{q}_rel']), 'same_mesh': False})['status'] == 'pass'
            assert one(k, {f'{q}_S_per_cm': ref * (1 + 1.01 * tol[f'{q}_rel']), 'same_mesh': False})['status'] == 'fail'


# ------------------------------------------------------------------ the exported atlas entries agree with the scorer
def _entry_submission(c, e):
    res = {f'{q}_S_per_cm': e[f'{q}_S_per_cm'] for q in ('G_DD', 'G_SD') if e.get(f'{q}_S_per_cm') is not None}
    if not res:
        return None
    res['same_mesh'] = bool(e['same_mesh'])
    if e.get('mesh_id'):
        res['mesh_id'] = e['mesh_id']
    return res


def test_known_bad_outputs_fail_and_cross_paths_do_not():
    for k in THREE:
        c = CASES[k]
        for e in c['known_bad_outputs']:
            s = _entry_submission(c, e)
            if s is not None:
                assert one(k, s)['status'] == 'fail', (k, e['name'])
        for e in c['cross_paths']:
            s = _entry_submission(c, e)
            if s is not None:
                assert one(k, s)['status'] in ('pass', 'not_scorable'), (k, e['name'])


def test_b0_companions_bracket_scoring():
    for k in BR:
        c = CASES[k]
        g0 = c['observables']['G0_S_per_cm']['value']
        assert one(k, dict(G_SD_S_per_cm=g0, same_mesh=True))['status'] == 'pass'
        assert one(k, dict(G_DD_S_per_cm=g0, same_mesh=True))['status'] == 'pass'
        assert one(k, dict(G0_S_per_cm=g0, same_mesh=False))['status'] == 'not_scorable'
        for e in c['known_bad_outputs']:
            assert one(k, dict(G0_S_per_cm=e['G0_S_per_cm'], same_mesh=True))['status'] == 'fail', (k, e['name'])
        for e in c['cross_paths']:
            if e['same_mesh']:
                assert one(k, dict(G0_S_per_cm=e['G0_S_per_cm'], same_mesh=True))['status'] == 'pass', (k, e['name'])


def test_every_three_terminal_case_carries_the_atlas_paths():
    for k in THREE:
        c = CASES[k]
        lev = c['mesh']['level']
        names = {e['name'] for e in c['known_bad_outputs'] + c['cross_paths']}
        p2_levels = sorted(int(s['path'].rsplit('_lev', 1)[1][0]) for s in c['provenance']['sources'] if '/P2_' in s['path'])
        assert 1 in p2_levels, k                                   # every state has its level-1 P2 run
        for lv in p2_levels:
            for nm in ('devsim_drain_contact_current', 'circuit_node_current', 'weighting_function_current', 'devsim_source_contact_current'):
                assert f'{nm}_richardson_lev{lv}' in names, (k, nm, lv)
        if lev not in p2_levels:                                   # a missing P2 run is stated in the case notes (amendment A8)
            assert any(f'P2 (drift-diffusion) at level {lev}' in n for n in c['notes']), k
        assert f'drain_side_flux_sum_lev{lev}' in names
        assert 'coupled_small_signal_G_E_lev1' in names and 'coupled_small_signal_G_M_lev1' in names
