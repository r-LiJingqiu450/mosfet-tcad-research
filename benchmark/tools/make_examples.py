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
