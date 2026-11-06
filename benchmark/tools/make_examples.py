#!/usr/bin/env python3
"""Write worked example submissions built from outputs already stored in the case files, and their score reports.

These are retrospective demonstrations: each example re-submits one known output path of the project's own
code (source-side reference, drain-side flux sum, coupled small-signal G_E, independent finite volume) so that
a reader can see what pass / fail / not_scorable look like. They are not independent validation.

Usage (repository root):  python benchmark/tools/make_examples.py
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import score as S  # noqa: E402

BENCH = os.path.dirname(HERE)
OUT = os.path.join(BENCH, 'examples')


def pick(c, name):
    for e in c['known_bad_outputs'] + c['cross_paths']:
        if e['name'] == name:
            return e
    return None


WP4_PROD = os.path.join(BENCH, 'wp4', 'results', 'production')
WP4_NOTE = ('values from benchmark/wp4/results/production/{files} (WP4 pre-registration '
            'benchmark/specs/WP4_Preregistration_External_Paths_2026-11-04.md; evaluation benchmark/wp4/results/evaluation_wp4.json)')
WP4_SAME_RUN_RTOL = 1e-9   # same mesh: the same-run network G_DD must reproduce the case's level-1 reference (extended equilibria shift it by < 1e-10)


def richardson(g1, g01):
    return (100.0 * g01 - g1) / 99.0


def wp4_examples(cases, tt, ref3):
    """DEVSIM default contact currents with extended precision (P2X), and DEVSIM built-in AC analysis (double, extended)."""
    ar = [k for k in tt if cases[k]['tier'] == 'A-R']
    files = {k: cases[k]['device']['state_id'] for k in ar}
    if not ar or not all(os.path.exists(os.path.join(WP4_PROD, f'P2X_{sid}_B1_lev1.json')) for sid in files.values()):
        return {}

    def same_run_check(k, d):
        mid, refs = ref3(cases[k], 1)
        g = d['network_same_run']['G_DD']
        if abs(g / refs['G_DD_S_per_cm'] - 1.0) > WP4_SAME_RUN_RTOL:
            raise SystemExit(f'{k}: WP4 record is not on the level-1 reference mesh/network ({g} vs {refs["G_DD_S_per_cm"]})')
        return mid

    p2x, acs = {}, {'double': {}, 'extended': {}}
    for k in ar:
        sid = files[k]
        d = json.load(open(os.path.join(WP4_PROD, f'P2X_{sid}_B1_lev1.json')))
        mid = same_run_check(k, d)
        s = d['symmetric']
        I_S = {st['V_D']: st['I_source'] for st in d['states']}
        p2x[k] = dict(G_DD_S_per_cm=richardson(s['1mV']['I_drain']['G'], s['0.1mV']['I_drain']['G']),
                      G_SD_S_per_cm=richardson(-(I_S[1e-3] - I_S[-1e-3]) / 2e-3, -(I_S[1e-4] - I_S[-1e-4]) / 2e-4),
                      same_mesh=True, mesh_id=mid)
        for prec in acs:
            path = os.path.join(WP4_PROD, f'AC_{sid}_B1_lev1_{prec}.json')
            if os.path.exists(path):
                a = json.load(open(path))
                acs[prec][k] = dict(G_DD_S_per_cm=a['ac']['G_AC'], same_mesh=True, mesh_id=same_run_check(k, a))
    out = {'three_terminal_devsim_contact_currents_extended_level1_mesh': dict(
        submission=dict(name='DEVSIM 2.11.2 default contact currents with extended precision (level-1 mesh)',
                        code='DEVSIM 2.11.2, extended_solver/extended_model/extended_equation = True (float128 models, equations and '
                             'Newton residuals; double-precision UMFPACK LU)',
                        current_evaluation='get_contact_current at the drain (G_DD) and at the source (G_SD); symmetric +-0.1/+-1 mV '
                                           'differences, Richardson',
                        notes=WP4_NOTE.format(files='P2X_<state>_B1_lev1.json') + '. G_DD is the registered quantity (X2); G_SD is post hoc.'),
        results=p2x)}
    for prec, how in (('double', 'double precision'), ('extended', 'extended precision (float128 assembly; one double-precision LU '
                                                                       'solve, no refinement)')):
        if len(acs[prec]) == len(ar):
            out[f'three_terminal_devsim_ac_small_signal_{prec}_level1_mesh'] = dict(
                submission=dict(name=f'DEVSIM 2.11.2 built-in AC analysis, {prec} precision (level-1 mesh)', code='DEVSIM 2.11.2',
                                current_evaluation=f"devsim.solve(type='ac', frequency=1e-6 Hz) with acreal = 1 V on the drain source; "
                                                   f"G_DD = -Re(V1.I); {how}",
                                notes=WP4_NOTE.format(files=f'AC_<state>_B1_lev1_{prec}.json')),
                results=acs[prec])
    return out


def main():
    cases = S.load_cases()
    zb = sorted(k for k, c in cases.items() if c['family'] in ('zero_bias_network', 'well_depth'))
    ex = {
        'source_side_reference_same_mesh': dict(
            submission=dict(name='source-side network current (reference path)', current_evaluation='cancellation-free source-side terminal sum'),
            results={k: dict(G0_S_per_cm=cases[k]['observables']['G0_S_per_cm']['value'], same_mesh=True) for k in zb}),
        'drain_side_flux_sum_same_mesh': dict(
            submission=dict(name='drain-side flux sum', current_evaluation='direct sum of edge currents at the drain contact'),
            results={k: dict(G0_S_per_cm=pick(cases[k], 'drain_side_flux_sum')['G0_S_per_cm'], same_mesh=True)
                     for k in zb if pick(cases[k], 'drain_side_flux_sum')}),
        'coupled_small_signal_G_E_same_mesh': dict(
            submission=dict(name='DEVSIM coupled-Jacobian small-signal G_E', current_evaluation='-x[V1.I] of the equilibrated SuperLU solve'),
            results={k: dict(G0_S_per_cm=pick(cases[k], 'coupled_small_signal_G_E')['G0_S_per_cm'], same_mesh=True)
                     for k in zb if pick(cases[k], 'coupled_small_signal_G_E')}),
        'independent_finite_volume_other_mesh': dict(
            submission=dict(name='independent finite-volume solver (phase3/c11/fd2d_indep.py)', current_evaluation='source-side network current on its own mesh'),
            results={k: dict(dL_nm=pick(cases[k], 'independent_finite_volume')['dL_nm'], same_mesh=False)
                     for k in zb if pick(cases[k], 'independent_finite_volume')}),
    }
    # schema 1.1 three-terminal cases (A-R / tier-C / HO, body B1)
    tt = sorted(k for k, c in cases.items() if c['family'] == 'zero_bias_3terminal')

    def ref3(c, level=None):
        smr = c['observables']['same_mesh_references']
        mid = c['mesh']['mesh_id'] if level is None else [m for m, v in smr.items() if v['level'] == level][0]
        return mid, {q: smr[mid][q] for q in ('G_DD_S_per_cm', 'G_SD_S_per_cm', 'G_DB_S_per_cm')}

    def at_level(c, name):
        return pick(c, f"{name}_lev{c['mesh']['level']}")

    ex.update({
        'three_terminal_reference_same_mesh': dict(
            submission=dict(name='P1 network reference values (B1)', current_evaluation='edge energy (G_DD), source side (G_SD), body side (G_DB)'),
            results={k: dict(**ref3(cases[k])[1], same_mesh=True) for k in tt}),
        'three_terminal_devsim_drain_contact_current_same_mesh': dict(
            submission=dict(name='DEVSIM default drain contact current', current_evaluation='get_contact_current, symmetric +-0.1/+-1 mV, Richardson'),
            results={k: dict(G_DD_S_per_cm=at_level(cases[k], 'devsim_drain_contact_current_richardson')['G_DD_S_per_cm'], same_mesh=True)
                     for k in tt if at_level(cases[k], 'devsim_drain_contact_current_richardson')}),
        'three_terminal_circuit_node_current_same_mesh': dict(
            submission=dict(name='DEVSIM circuit-node current', current_evaluation='-x[V1.I], symmetric +-0.1/+-1 mV, Richardson'),
            results={k: dict(G_DD_S_per_cm=at_level(cases[k], 'circuit_node_current_richardson')['G_DD_S_per_cm'], same_mesh=True)
                     for k in tt if at_level(cases[k], 'circuit_node_current_richardson')}),
        'three_terminal_coupled_small_signal_G_E_level1_mesh': dict(
            submission=dict(name='DEVSIM coupled-Jacobian small-signal G_E (level-1 mesh)', current_evaluation='-x[V1.I] of the equilibrated SuperLU solve'),
            results={k: dict(G_DD_S_per_cm=pick(cases[k], 'coupled_small_signal_G_E_lev1')['G_DD_S_per_cm'], same_mesh=True,
                             mesh_id=ref3(cases[k], 1)[0]) for k in tt}),
        'three_terminal_independent_finite_volume_other_mesh': dict(
            submission=dict(name='independent finite-volume solver (benchmark/wp2/fv_rdeck.py, refined mesh)', current_evaluation='network terminal currents'),
            results={k: dict(G_DD_S_per_cm=pick(cases[k], 'independent_finite_volume_refined')['G_DD_S_per_cm'],
                             G_SD_S_per_cm=pick(cases[k], 'independent_finite_volume_refined')['G_SD_S_per_cm'], same_mesh=False)
                     for k in tt if pick(cases[k], 'independent_finite_volume_refined')}),
    })
    # WP4 external paths (A-R cases, level-1 mesh), read from the WP4 production records. They are scored as submissions and
    # are not written into the case files (WP4 pre-registration section 6: WP1 and WP2 case files stay unchanged).
    ex.update(wp4_examples(cases, tt, ref3))
    os.makedirs(OUT, exist_ok=True)
    for name, sub in ex.items():
        with open(os.path.join(OUT, f'{name}.submission.json'), 'w', encoding='utf-8', newline='\n') as f:
            f.write(json.dumps(sub, indent=1, sort_keys=True, allow_nan=False) + '\n')
        rep = S.score(sub, cases)
        with open(os.path.join(OUT, f'{name}.report.json'), 'w', encoding='utf-8', newline='\n') as f:
            f.write(json.dumps(rep, indent=1, sort_keys=True, allow_nan=False) + '\n')
        print(f'{name:40s}', json.dumps(rep['summary']))
    return 0


if __name__ == '__main__':
    sys.exit(main())
