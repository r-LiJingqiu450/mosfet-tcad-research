#!/usr/bin/env python3
"""POST-HOC DESCRIPTIVE ANALYSIS of the WP2 production outputs (not pre-registered; used in no registered decision).

1. DEVSIM post-solve contact currents (get_contact_current) at the n+ contacts: are the raw values at a given bias
   the same for different gate voltages on the same mesh? Identical values across states with channel conductances
   that differ by orders of magnitude mean the contact current carries no channel information at that level.
2. The signal they would have to resolve: G_DD * V_D, against the raw current magnitude.
3. Summary of the forms diagnostic (benchmark/wp2/forms_diagnostic.py) if its outputs exist.

Reads benchmark/wp2/results/production/P2_*.json and P1_*.json; writes benchmark/wp2/results/atlas_observations.json.
"""
from __future__ import annotations

import collections
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PROD = os.path.join(HERE, 'results', 'production')
DIAG = os.path.join(HERE, 'results', 'diagnostics')


def geom(a):
    dev = 'dev1' if a['na'] == 1e18 else 'dev2'
    return f"{dev}_L{a['L']:g}_lev{a['level']}"


def main():
    runs = []
    for p in sorted(glob.glob(os.path.join(PROD, 'P2_*.json'))):
        d = json.load(open(p))
        a = d['args']
        runs.append(dict(file=os.path.basename(p), state=a['id'], body=a['body'], mesh=geom(a), vg=a['vg'],
                         G=d['network_same_run']['G_DD'], states={s['V_D']: s for s in d['states']}))
    groups = collections.defaultdict(list)
    for r in runs:
        groups[(r['mesh'], r['body'])].append(r)
    ident = []
    for (mesh, body), rs in sorted(groups.items()):
        for vd in (1e-4, -1e-4, 1e-3, -1e-3):
            for key in ('I_drain', 'I_source'):
                vals = collections.Counter(repr(r['states'][vd][key]) for r in rs)
                modal, k = vals.most_common(1)[0]
                sharing = sorted(r['state'] for r in rs if repr(r['states'][vd][key]) == modal)
                Gs = [r['G'] for r in rs if r['state'] in sharing]
                ident.append(dict(mesh=mesh, body=body, V_D=vd, current=key, n_runs=len(rs), n_distinct=len(vals),
                                  modal_value_A_per_cm=float(modal), states_sharing_modal_value=sharing,
                                  G_DD_range_of_sharing_states=[min(Gs), max(Gs)],
                                  signal_range_A_per_cm=[min(Gs) * abs(vd), max(Gs) * abs(vd)]))
    # per run: absolute error of the DEVSIM post-solve contact currents, taking the circuit-node current (accurate to ~1e-9
    # relative, expectation E-CIRC) as the finite-bias current; DEVSIM's sign convention is I_drain = +I, I_source = -I
    per_run = []
    for r in runs:
        s = r['states']
        per_run.append(dict(run=r['file'][:-5], mesh=r['mesh'], body=r['body'], G_DD=r['G'], signal_1mV=r['G'] * 1e-3,
                            max_abs_err_drain=max(abs(s[v]['I_drain'] - s[v]['I_circuit']) for v in s),
                            max_abs_err_source=max(abs(s[v]['I_source'] + s[v]['I_circuit']) for v in s),
                            I_drain_raw={f'{v:+g}': s[v]['I_drain'] for v in sorted(s)},
                            I_source_raw={f'{v:+g}': s[v]['I_source'] for v in sorted(s)},
                            I_circuit_raw={f'{v:+g}': s[v]['I_circuit'] for v in sorted(s)},
                            drain_edge_abs_sum={f'{v:+g}': s[v]['drain_edge_abs_sum'] for v in sorted(s)}))
    diag = {}
    for p in sorted(glob.glob(os.path.join(DIAG, 'forms_diag_*.json'))):
        d = json.load(open(p))
        diag[os.path.basename(p)[:-5]] = dict(
            mesh=d['mesh'], energy_vs_source_rel=d['energy_vs_source_rel'],
            source_half_channel_cuts=d['source_half_channel_cuts'], channel_median_rel_vs_energy=d['channel_cuts']['median_rel_vs_energy'],
            drain_side_cuts_max_abs_rel_vs_source=d['drain_side_cuts_max_abs_rel_vs_source'], edge_conductance=d['edge_conductance'])
    p1_forms = []
    for p in sorted(glob.glob(os.path.join(PROD, 'P1_*_B1_lev*.json'))):
        d = json.load(open(p))
        n = d['network']
        p1_forms.append(dict(run=os.path.basename(p)[:-5], G_DD=n['G_DD'], G_DD_rel_forms=n['G_DD_rel_forms'],
                             G_SD_rel_forms=n['G_SD_rel_forms'], abs_forms_A_per_V_cm=n['G_DD'] * abs(n['G_DD_rel_forms'])))
    floor = {}
    for x in per_run:
        f = floor.setdefault(x['mesh'], dict(runs=0, drain=[], source=[], G_DD=[]))
        f['runs'] += 1
        f['drain'].append(x['max_abs_err_drain'])
        f['source'].append(x['max_abs_err_source'])
        f['G_DD'].append(x['G_DD'])
    floor = {m: dict(runs=f['runs'], max_abs_err_drain_range=[min(f['drain']), max(f['drain'])],
                     max_abs_err_source_range=[min(f['source']), max(f['source'])], G_DD_range=[min(f['G_DD']), max(f['G_DD'])])
             for m, f in sorted(floor.items())}
    out = dict(note='post-hoc descriptive analysis; not pre-registered; used in no registered decision',
               contact_current_identity=ident, noise_floor_by_mesh=floor, per_run=per_run, forms_diagnostic=diag, p1_G_DD_forms=p1_forms)
    with open(os.path.join(HERE, 'results', 'atlas_observations.json'), 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=1)
    print('DEVSIM raw contact currents shared across gate voltages (same mesh, body, bias):')
    for x in ident:
        if x['n_runs'] > 1:
            print(f"  {x['mesh']:13s} {x['body']} V_D={x['V_D']:+.0e} {x['current']:8s}: {len(x['states_sharing_modal_value'])}/{x['n_runs']} runs "
                  f"identical ({x['modal_value_A_per_cm']:+.6e} A/cm); their G_DD span {x['G_DD_range_of_sharing_states'][0]:.2e}.."
                  f"{x['G_DD_range_of_sharing_states'][1]:.2e} S/cm")
    print('DEVSIM contact-current absolute error (vs circuit-node current), max over the four biases, by mesh:')
    for m, f in floor.items():
        print(f"  {m:13s} {f['runs']} runs: drain {f['max_abs_err_drain_range'][0]:.1e}..{f['max_abs_err_drain_range'][1]:.1e} A/cm, "
              f"source {f['max_abs_err_source_range'][0]:.1e}..{f['max_abs_err_source_range'][1]:.1e} A/cm; G_DD {f['G_DD_range'][0]:.1e}..{f['G_DD_range'][1]:.1e} S/cm")
    print('P1 G_DD forms (B1):', ', '.join(f"{x['run'][3:]} {x['G_DD_rel_forms']:+.1e}" for x in p1_forms))
    for k, v in diag.items():
        print(f"  {k}: energy vs source {v['energy_vs_source_rel']:+.2e}; source-half cuts vs source <= "
              f"{v['source_half_channel_cuts']['max_abs_rel_vs_source']:.1e}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
