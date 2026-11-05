#!/usr/bin/env python3
"""Print the WP2 report tables (markdown) from benchmark/wp2/results/evaluation.json and atlas_observations.json,
so that no number in the report is transcribed by hand.

Usage (repository root):  python benchmark/wp2/report_tables.py > benchmark/wp2/results/report_tables.md
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, 'results')


def e(x, p=2):
    return 'n/a' if x is None else f'{x:+.{p}e}'


def main():
    ev = json.load(open(os.path.join(RES, 'evaluation.json')))
    at = json.load(open(os.path.join(RES, 'atlas_observations.json'))) if os.path.exists(os.path.join(RES, 'atlas_observations.json')) else None
    out = []
    w = out.append
    w('### H1: reference quality (A-R, body B1)\n')
    w('| State | ψ̂ (P1 L2) | G_DD ref (S/cm) | G_DD P4/L2−1 | G_DD L2/L1−1 | G_SD ref (S/cm) | G_SD P4/L2−1 | G_SD L2/L1−1 | max \\|forms\\| L0/L1/L2 | max KCL | B0 bracket L0/L1/L2 | H1 |')
    w('|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---|---|')
    for sid, r in ev['H1'].items():
        if r.get('missing'):
            w(f'| {sid} | missing: {", ".join(r["missing"])} |||||||||||')
            continue
        c, ref, sc = r['checks'], r['reference'], r['self_consistency']
        forms = '/'.join(f"{max(abs(x['G_DD_forms']), abs(x['G_SD_forms'])):.1e}" for x in sc)
        kcl = max(x['kcl'] for x in sc)
        br = '/'.join('ok' if x['B0_bracket_ok'] else 'FAIL' for x in sc)
        w(f"| {sid} | {r['psi_hat_P1_lev2']:.3f} | {ref['G_DD']['value']:.6e} | {e(c['G_DD']['cross_P4_vs_L2'])} | {e(c['G_DD']['refine_L2_vs_L1'])} | "
          f"{ref['G_SD']['value']:.6e} | {e(c['G_SD']['cross_P4_vs_L2'])} | {e(c['G_SD']['refine_L2_vs_L1'])} | {forms} | {kcl:.1e} | {br} | "
          f"{'**pass**' if r['pass'] else '**FAIL**'} |")
    s = ev['H1_summary']
    w(f"\nH1 summary: {s['n_pass']} of 6 pass ({', '.join(s['passed'])}); kill (fewer than 4) = {s['kill']}.\n")
    w('### Reference values and tolerances (section 7)\n')
    w('| State | G_DD tol (rel) | G_SD tol (rel) | u_disc G_DD / G_SD | u_round G_DD / G_SD | u_solver | observed order G_DD / G_SD | P4 refined vs base G_DD | ΔL_SD (nm) | ΔL tol (nm) | G_DB (S/cm) | G_DB/G_DD |')
    w('|---|---:|---:|---|---|---:|---|---:|---:|---:|---:|---:|')
    for sid, r in ev['H1'].items():
        if r.get('missing'):
            continue
        ref = r['reference']
        dd, sd, db = ref['G_DD'], ref['G_SD'], ref['G_DB']
        oo = '/'.join('n/a' if x['observed_order'] is None else f"{x['observed_order']:.2f}" for x in (dd, sd))
        w(f"| {sid} | {dd['tol_rel']:.2e} | {sd['tol_rel']:.2e} | {dd['u_disc']:.1e} / {sd['u_disc']:.1e} | {dd['u_round']:.1e} / {sd['u_round']:.1e} | "
          f"{dd['u_solver']:.1e} | {oo} | {e(dd['P4_self_convergence'])} | {r['dL_SD_nm']:.3f} | {r['dL_tol_nm']:.3f} | {db['value']:.3e} | {db['value'] / dd['value']:.2e} |")
    w('\n### H-B: body effect at level 2 (B1 vs B0)\n')
    w('| State | Δ_DD | band | in | Δ_SD | band | in | G_DB/G₀ | P4 Δ_DD | P4 Δ_SD |')
    w('|---|---:|---|---|---:|---|---|---:|---:|---:|')
    for sid, r in ev['HB'].items():
        if r.get('missing'):
            w(f'| {sid} | missing |||||||||')
            continue
        p4 = r.get('P4') or {}
        w(f"| {sid} | {e(r['Delta_DD'])} | [{r['band_DD'][0]:.1e}, {r['band_DD'][1]:.1e}] | {'yes' if r['in_DD'] else 'NO'} | {e(r['Delta_SD'])} | "
          f"[{r['band_SD'][0]:.1e}, {r['band_SD'][1]:.1e}] | {'yes' if r['in_SD'] else 'NO'} | {r['G_DB_over_G0']:.2e} | {e(p4.get('Delta_DD'))} | {e(p4.get('Delta_SD'))} |")
    s = ev['HB_summary']
    w(f"\nH-B summary: Δ_DD in band {s['in_band_DD']}/{s['n']}, Δ_SD in band {s['in_band_SD']}/{s['n']}; supported (≥ 5 of 6 each) = {s['supported']}.\n")
    w('### H2: held-out robustness (level 1)\n')
    w('| State | ψ̂ | G_DD (S/cm) | G_SD (S/cm) | G_DB (S/cm) | (a) G_DD forms | (b) KCL | (c) B0 bracket | H2 |')
    w('|---|---:|---:|---:|---:|---:|---:|---|---|')
    for sid, r in ev['H2'].items():
        if r.get('missing'):
            w(f'| {sid} | missing ||||||||')
            continue
        w(f"| {sid} | {r['psi_hat']:.3f} | {r['G_DD']:.6e} | {r['G_SD']:.6e} | {r['G_DB']:.3e} | {r['a_forms']:.1e} | {r['b_kcl']:.1e} | "
          f"{'ok' if r['c_B0_bracket'] else 'FAIL'} | {'**pass**' if r['pass'] else '**FAIL**'} |")
    s = ev['H2_summary']
    w(f"\nH2 summary: {s['n_pass']} of {s['n']} pass; complete = {s['complete']}; kill = {s['kill']}.\n")
    w('### P-DD and E-CIRC: drift–diffusion at V_D = ±0.1 and ±1 mV (relative error vs the P1 network G_DD on the same mesh)\n')
    w('| Run | G·1 mV (A/cm) | predicted | DEVSIM drain ±1 mV | ±0.1 mV | Richardson | P-DD | circuit ±1 mV | circuit Richardson | weighting Richardson | E-CIRC |')
    w('|---|---:|---|---:|---:|---:|---|---:|---:|---:|---|')
    for r in ev['P2_rows']:
        pr = {'fail': '> 1e-2', 'pass': '≤ 1e-2', 'none': 'none'}[r['prediction']]
        ok = {True: 'confirmed', False: '**refuted**', None: '—'}[r['pdd_ok']]
        w(f"| {r['state']} {r['body']} L{r['level']} | {r['G_times_1mV']:.1e} | {pr} | {e(r['default_1mV'])} | {e(r['default_0.1mV'])} | "
          f"{e(r['default_richardson'])} | {ok} | {e(r['circuit_1mV'])} | {e(r['circuit_richardson'], 1)} | {e(r['weighted_richardson'], 1)} | "
          f"{'yes' if r['ecirc_ok'] else '**no**'} |")
    w('\nP-DD by state: ' + '; '.join(f"{k} {v['outcome']}" for k, v in ev['PDD_by_state'].items()) +
      f".  E-CIRC: all = {ev['ECIRC']['ecirc_all']}, worst |rel| = {ev['ECIRC']['ecirc_worst']:.1e}.\n")
    w('### P-SS: DEVSIM coupled small-signal at level 1 (relative error vs the P1 network G_DD)\n')
    w('| State | body | G_E | G_M | κ_G | E_dwr | predicted > 1e-2 | outcome |')
    w('|---|---|---:|---:|---:|---:|---|---|')
    for r in ev['PSS_rows']:
        oc = '—' if r['ok'] is None else ('confirmed' if r['ok'] else '**refuted**')
        w(f"| {r['state']} | {r['body']} | {e(r['G_E_rel'])} | {e(r['G_M_rel'])} | {r['kappa_G']:.1e} | {r['E_dwr']:.2f} | {'yes' if r['predicted'] else 'no'} | {oc} |")
    s = ev['PSS_summary']
    w(f"\nP-SS summary: {s['n_confirmed']} of {s['n_predicted']} predicted cases confirmed.\n")
    w('### ψ̂ check (P1, B1, level 1)\n')
    w('| ' + ' | '.join(ev['psi_hat']) + ' |')
    w('|' + '---:|' * len(ev['psi_hat']))
    w('| ' + ' | '.join(f"{v['P1_lev1']:.3f} (target {v['expected']})" for v in ev['psi_hat'].values()) + ' |')
    L = ev['ledger']
    w(f"\nLedger: {L['runs']} runs, failures {L['failures']}, total wall {L['total_wall_h']:.2f} h, max RSS {L['max_rss_MB']} MB.\n")
    rows = [json.loads(x) for x in open(os.path.join(RES, 'ledger.jsonl'))]
    w('### Unsuccessful attempts (kept in the ledger)\n')
    w('| Job | start (UTC) | exit | wall (s) | peak RSS (MB) | cause |')
    w('|---|---|---:|---:|---:|---|')
    cause = {-9: 'killed by the memory cgroup (OOM)', -15: 'stopped by me (SIGTERM) to prevent a cascade of OOM kills (amendment A9)'}
    for r in rows:
        if r['exit'] != 0:
            w(f"| {r['job']} | {r['start_utc']} | {r['exit']} | {r['wall_s']} | {r['maxrss_MB']} | {cause.get(r['exit'], 'see log')} |")
    n_env = sum(1 for r in rows if r.get('env_extra'))
    w(f"\nRuns with `MALLOC_MMAP_THRESHOLD_=131072` (amendment A9): {n_env} of {len(rows)} ledger lines.\n")
    if at:
        w('### Post-hoc: absolute error of DEVSIM post-solve contact currents (max over ±0.1, ±1 mV; vs the circuit-node current)\n')
        w('| Mesh | runs | drain (A/cm) | source (A/cm) | G_DD of those runs (S/cm) |')
        w('|---|---:|---|---|---|')
        for m, f in at['noise_floor_by_mesh'].items():
            w(f"| {m} | {f['runs']} | {f['max_abs_err_drain_range'][0]:.1e} – {f['max_abs_err_drain_range'][1]:.1e} | "
              f"{f['max_abs_err_source_range'][0]:.1e} – {f['max_abs_err_source_range'][1]:.1e} | {f['G_DD_range'][0]:.1e} – {f['G_DD_range'][1]:.1e} |")
        w('')
        w('### Post-hoc: DEVSIM raw contact currents identical across gate voltages (same mesh, body and bias)\n')
        w('| Mesh | body | V_D | current | runs identical / runs | value (A/cm) | G_DD span of those runs (S/cm) |')
        w('|---|---|---:|---|---|---:|---|')
        for x in at['contact_current_identity']:
            if x['n_runs'] > 1 and len(x['states_sharing_modal_value']) > 1:
                w(f"| {x['mesh']} | {x['body']} | {x['V_D']:+.0e} | {x['current']} | {len(x['states_sharing_modal_value'])}/{x['n_runs']} | "
                  f"{x['modal_value_A_per_cm']:+.6e} | {x['G_DD_range_of_sharing_states'][0]:.2e} – {x['G_DD_range_of_sharing_states'][1]:.2e} |")
        if at['forms_diagnostic']:
            w('\n### Post-hoc: which G_DD form carries the round-off (body B0, cuts through the channel)\n')
            w('| Run | nodes | energy / (−I_S) − 1 | source-half channel cuts vs −I_S (max) | cuts vs energy (range) | max c_e / G |')
            w('|---|---:|---:|---:|---|---:|')
            for k, v in at['forms_diagnostic'].items():
                h = v['source_half_channel_cuts']
                w(f"| {k.replace('forms_diag_', '')} | {v['mesh']['nx'] * v['mesh']['ny']} | {v['energy_vs_source_rel']:+.2e} | {h['max_abs_rel_vs_source']:.1e} | "
                  f"[{h['min_rel_vs_energy']:+.2e}, {h['max_rel_vs_energy']:+.2e}] | {v['edge_conductance']['max_over_G']:.1e} |")
    print('\n'.join(out))
    return 0


if __name__ == '__main__':
    sys.exit(main())
