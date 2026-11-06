#!/usr/bin/env python3
"""Generate the per-task kit data of the WP6 agent study from the WP2 deck (unchanged code paths).

For each task state (the 6 A-R and 5 HO three-terminal cases, body B1, published level-1 mesh; plus the non-benchmark
smoke state SMOKE-sub on the level-0 mesh) this script, in a fresh DEVSIM process:
  1. builds the WP2 deck and converges the zero-bias state (benchmark/wp2/rdeck_common.build);
  2. computes the network G_DD with the WP2 P1 code (benchmark/wp2/rdeck_probe.analyse) and checks it against the
     case reference on that mesh (observables.same_mesh_references); the smoke state has no case, its value is kept;
  3. writes benchmark/wp6/tasks/<task>/ : device.json (neutral description), mesh.json (node lines and rules only),
     zero_bias_state.npz (psi_si, psi_ox, n_si: the converged arrays).
benchmark/wp6/tasks/tasks.json maps task ids to cases and holds the references; it is never copied to a workspace.

Usage (repository root):  python benchmark/wp6/harness/make_kit_states.py [--only T01]
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WP6 = os.path.dirname(HERE)
ROOT = os.path.dirname(os.path.dirname(WP6))
TASKS = os.path.join(WP6, 'tasks')
PY = sys.executable
ENV = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MALLOC_MMAP_THRESHOLD_='131072',
           DEVSIM_MATH_LIBS='/usr/lib/x86_64-linux-gnu/libblas.so.3:/usr/lib/x86_64-linux-gnu/liblapack.so.3')
# task order: A-R then HO by state id (the ids carry no information to the agent; workspaces use run ids)
ORDER = ['AR-1', 'AR-2', 'AR-3', 'AR-4', 'AR-5', 'AR-6', 'HO-1', 'HO-2', 'HO-3', 'HO-4', 'HO-5']
SMOKE = dict(task='SMOKE', state_id='SMOKE-sub', case=None, L_um=1.0, N_A_cm3=1e18, t_ox_nm=2.0, V_G_V=0.55, level=0,
             mesh_id='dev1_L1_lev0', reference_G_DD_S_per_cm=None)


def task_table():
    rows = []
    cases = {}
    for p in sorted(glob.glob(os.path.join(ROOT, 'benchmark', 'cases', '*.json'))):
        c = json.load(open(p))
        if c['family'] == 'zero_bias_3terminal':
            cases[c['device']['state_id']] = (os.path.relpath(p, ROOT), c)
    for k, sid in enumerate(ORDER, 1):
        path, c = cases[sid]
        smr = c['observables']['same_mesh_references']
        mid = [m for m, v in smr.items() if v['level'] == 1][0]
        d = c['device']
        rows.append(dict(task=f'T{k:02d}', state_id=sid, case=c['id'], case_file=path, tier=c['tier'],
                         L_um=d['geometry']['L_um'], N_A_cm3=d['doping']['N_A_cm3'], t_ox_nm=d['geometry']['t_ox_nm'],
                         V_G_V=d['bias']['V_G_V'], level=1, mesh_id=mid, reference_G_DD_S_per_cm=smr[mid]['G_DD_S_per_cm'],
                         current_at_1mV_A_per_cm=smr[mid]['G_DD_S_per_cm'] * 1e-3))
    return rows


def worker(t):
    """Runs in a fresh process: WP2 deck -> converged arrays + network G_DD."""
    import numpy as np
    sys.path.insert(0, os.path.join(ROOT, 'benchmark', 'wp2'))
    import rdeck_common as RC
    import rdeck_probe as RP
    from rdeck_common import D, ds
    _, geo, sol0, rec = RC.build(t['L_um'], t['level'], t['t_ox_nm'], t['N_A_cm3'], t['V_G_V'], 'B1')
    out, _ = RP.analyse(geo, sol0, 'B1')
    arr = dict(psi_si=np.array(ds.get_node_model_values(device=D.DEV, region=D.SI, name='Potential')),
               psi_ox=np.array(ds.get_node_model_values(device=D.DEV, region=D.OX, name='Potential')),
               n_si=np.array(ds.get_node_model_values(device=D.DEV, region=D.SI, name='Electrons')))
    d = os.path.join(TASKS, t['task'])
    os.makedirs(d, exist_ok=True)
    np.savez(os.path.join(d, 'zero_bias_state.npz'), **arr)
    sha = {k: hashlib.sha256(v.tobytes()).hexdigest() for k, v in arr.items()}
    print(json.dumps(dict(task=t['task'], G_DD=out['G_DD'], G_DD_second=out['G_DD_second'], psi_hat=rec['onecol']['psi_hat'],
                          equilibrium_iterations=rec['equilibrium']['iterations'], array_sha256=sha,
                          n_si=int(len(arr['psi_si'])), n_ox=int(len(arr['psi_ox'])))), flush=True)


def device_json(t, mesh):
    na = t['N_A_cm3']
    vt, ni = 1.380649e-23 * 300.0 / 1.602176634e-19, 1e10
    p0 = .5 * (na + math.sqrt(na * na + 4 * ni * ni))
    psi_b0 = -vt * math.log(p0 / ni)
    L = t['L_um']
    return dict(
        description='Two-dimensional n-channel MOSFET (quantities per unit width) with resolved n+ source and drain wells',
        geometry=dict(channel_length_L_um=L, oxide_thickness_nm=t['t_ox_nm'], silicon_x_range_um=[-0.3, L + 0.3],
                      gate_x_range_um=[0.0, L], silicon_depth_um=1.2, well_depth_h_w_um=0.5,
                      source_drain_contact_depth_um=mesh['contact_extent_cm'] * 1e4,
                      coordinates='x along the channel, source side first; y = depth into silicon, y = 0 at the Si/SiO2 '
                                  'interface, oxide at -t_ox <= y <= 0'),
        doping=dict(N_A_cm3=na, N_D_peak_cm3=1e20, h_w_um=0.5, sigma_nm=15.0, sigma_y_nm=15.0,
                    profile='uniform acceptors N_A plus donors N_D(x, y) = 1e20 f_x(x) f_y(y), with '
                            'f_x = 0.5 erfc((x - x_c)/(sqrt(2) sigma)) + 0.5 erfc(((L - x) - x_c)/(sqrt(2) sigma)), '
                            'x_c = -sqrt(2) sigma erfcinv(2 N_A / 1e20) (so N_D = N_A at x = 0 and x = L), '
                            'f_y = 0.5 erfc((y - h_w)/(sqrt(2) sigma_y))'),
        physics=dict(T_K=300.0, V_T_V=vt, n_i_cm3=ni, eps_si_rel=11.7, eps_ox_rel=3.9, mu_n_cm2_per_Vs=1350.0,
                     mobility='constant', statistics='Boltzmann, n = n_i exp((psi - phi_n)/V_T)',
                     holes='held at equilibrium, p = n_i exp(-psi/V_T) (not a solution variable)',
                     generation_recombination='none', psi_reference='intrinsic level (psi = 0 where n = p = n_i)',
                     psi_b0_V=psi_b0, oxide='SiO2 without charge'),
        contacts=dict(gate='oxide top y = -t_ox over 0 <= x <= L: psi = V_G + psi_b0 (flat-band voltage 0)',
                      source='silicon edge x = -0.3 um, 0 <= y <= contact depth: psi = V_T asinh(N_net/(2 n_i)), n = n_i exp(psi/V_T)',
                      drain='silicon edge x = L + 0.3 um, 0 <= y <= contact depth: psi = V_T asinh(N_net/(2 n_i)) + V_D, '
                            'n = n_i exp((psi - V_D)/V_T); driven by the circuit voltage source V1',
                      body='silicon bottom y = 1.2 um, full width: psi = psi_b0, n = n_i exp(psi_b0/V_T)',
                      other_boundaries='zero normal electric displacement and zero normal electron current',
                      oxide_corners='the two Si/SiO2 corner nodes at the far silicon edges are tied to the adjacent '
                                    'source/drain contact potential'),
        bias=dict(V_G_V=t['V_G_V'], V_S_V=0.0, V_B_V=0.0, V_D_V=0.0),
        mesh=dict(file='mesh.json', nx=mesh['nx'], ny=mesh['ny'],
                  note='fixed tensor-product mesh; each rectangle is split into two triangles'),
        units='cm, V, cm^-3; currents in A/cm and conductances in S/cm (per cm of device width)')


def mesh_json(mesh_id):
    m = json.load(open(os.path.join(ROOT, 'benchmark', 'meshes', f'mesh_{mesh_id}.json')))
    keep = ('x_cm', 'y_oxide_cm', 'y_si_cm', 'jw_si', 'kif', 'nx', 'ny', 'gate_columns', 'contact_extent_cm', 'x_junctions_cm',
            'units')
    out = {k: m[k] for k in keep}
    out['rules'] = {k: v for k, v in m['rules'].items() if k != 'triangles'}
    out['rules']['triangles'] = ('rectangle (i,j)-(i+1,j+1) is split into [(i,j),(i+1,j),(i+1,j+1)] and '
                                 '[(i,j),(i+1,j+1),(i,j+1)]')
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--worker', help=argparse.SUPPRESS)
    ap.add_argument('--only')
    a = ap.parse_args()
    if a.worker:
        worker(json.loads(a.worker))
        return 0
    rows = task_table() + [SMOKE]
    os.makedirs(TASKS, exist_ok=True)
    results = {}
    for t in rows:
        if a.only and t['task'] != a.only:
            continue
        mesh = mesh_json(t['mesh_id'])
        d = os.path.join(TASKS, t['task'])
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, 'mesh.json'), 'w') as f:
            json.dump(mesh, f, indent=1)
        with open(os.path.join(d, 'device.json'), 'w') as f:
            json.dump(device_json(t, mesh), f, indent=1)
        p = subprocess.run([PY, os.path.abspath(__file__), '--worker', json.dumps(t)], cwd=ROOT, env=ENV,
                           capture_output=True, text=True)
        line = [x for x in p.stdout.splitlines() if x.startswith('{')]
        if p.returncode != 0 or not line:
            print(p.stdout[-2000:], p.stderr[-2000:])
            raise SystemExit(f"{t['task']}: worker failed ({p.returncode})")
        r = json.loads(line[-1])
        ref = t['reference_G_DD_S_per_cm']
        r['rel_vs_reference'] = None if ref is None else r['G_DD'] / ref - 1.0
        if ref is None:
            t['reference_G_DD_S_per_cm'] = r['G_DD']
            t['reference_note'] = 'network G_DD computed here (smoke state; no case file)'
        results[t['task']] = r
        print(f"{t['task']} {t['state_id']:9s} G_DD={r['G_DD']:.10e} rel_vs_case={r['rel_vs_reference']} psi_hat={r['psi_hat']:.2f}",
              flush=True)
    if not a.only:
        for t in rows:
            t['kit_generation'] = results[t['task']]
        with open(os.path.join(TASKS, 'tasks.json'), 'w') as f:
            json.dump(dict(note='WP6 task table; private to the harness (never copied into a workspace)', tasks=rows), f, indent=1)
    return 0


if __name__ == '__main__':
    sys.exit(main())
