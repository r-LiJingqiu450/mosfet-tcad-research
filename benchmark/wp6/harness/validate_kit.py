#!/usr/bin/env python3
"""Validate the WP6 task kit: the agent-facing deck reproduces the benchmark device and its reference.

For each task (fresh DEVSIM process, kit files in a temporary directory exactly as in a workspace):
  1. kit build: mosfet_deck.build() (mesh from mesh.json, physics, zero-bias state from zero_bias_state.npz, Newton);
  2. differences between the kit's converged state and the stored arrays (max |d psi|, max |dn/n|);
  3. the WP2 network G_DD (benchmark/wp2 code, unchanged) evaluated on the kit's state, against the task reference;
  4. with --dd TASK: the kit's DEVSIM drift-diffusion at V_D = +-0.1 and +-1 mV, comparing the circuit-node current
     (Richardson) and DEVSIM's default drain contact current with the reference (a correct and a known-bad path).
Writes benchmark/wp6/results/kit_validation.json.

Usage (repository root):  python benchmark/wp6/harness/validate_kit.py [--dd T02]
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
WP6 = os.path.dirname(HERE)
ROOT = os.path.dirname(os.path.dirname(WP6))
PY = sys.executable
ENV = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MALLOC_MMAP_THRESHOLD_='131072',
           DEVSIM_MATH_LIBS='/usr/lib/x86_64-linux-gnu/libblas.so.3:/usr/lib/x86_64-linux-gnu/liblapack.so.3')


def worker(t, kitdir, dd):
    import numpy as np
    sys.path.insert(0, kitdir)
    import mosfet_deck as M
    info = M.build()
    ds = M.ds
    st = np.load(os.path.join(kitdir, 'zero_bias_state.npz'))
    psi = np.array(ds.get_node_model_values(device='mos', region='si', name='Potential'))
    n = np.array(ds.get_node_model_values(device='mos', region='si', name='Electrons'))
    out = dict(task=t['task'], newton_iterations=info['newton_iterations'],
               max_dpsi_V=float(np.max(np.abs(psi - st['psi_si']))), max_rel_dn=float(np.max(np.abs(n / st['n_si'] - 1))))
    # WP2 network on the kit state
    sys.path.insert(0, os.path.join(ROOT, 'benchmark', 'wp2'))
    import rdeck_common as RC
    from rdeck_common import WP
    _, geo = RC.make_geo(t['L_um'], t['level'], t['t_ox_nm'], t['N_A_cm3'], t['V_G_V'])
    xs, ys = map(np.array, RC.D.coords(RC.D.SI))
    n0i, n1i, c, _, _, _ = WP.network(geo)
    Mk = RC.terminal_masks(geo, xs, ys)
    T = dict(source=np.where(Mk['source'])[0], drain=np.where(Mk['drain'])[0], body=np.where(Mk['body'] & ~Mk['source'] & ~Mk['drain'])[0])
    rd = RC.terminal_solve(RC.laplacian(len(xs), n0i, n1i, c), n0i, n1i, c, T, 'drain')
    out['G_DD_network_on_kit'] = rd['energy']
    out['rel_vs_reference'] = rd['energy'] / t['reference_G_DD_S_per_cm'] - 1.0
    if dd:
        cur = {}
        for v in (1e-4, -1e-4, 1e-3, -1e-3):
            ds.circuit_alter(name='V1', value=v)
            M.solve_dc()
            cur[v] = dict(circuit=-float(ds.get_circuit_node_value(solution='dcop', node='V1.I')),
                          contact=float(ds.get_contact_current(device='mos', contact='drain', equation='ElectronContinuityEquation')))
            ds.circuit_alter(name='V1', value=0.0)
            M.solve_dc()
        ref = t['reference_G_DD_S_per_cm']
        for k in ('circuit', 'contact'):
            g1 = (cur[1e-3][k] - cur[-1e-3][k]) / 2e-3
            g01 = (cur[1e-4][k] - cur[-1e-4][k]) / 2e-4
            out[f'dd_{k}'] = dict(sym_1mV_rel=g1 / ref - 1, sym_0p1mV_rel=g01 / ref - 1, richardson_rel=((100 * g01 - g1) / 99) / ref - 1)
        out['dd_raw_currents'] = {f'{v:+g}': cur[v] for v in cur}
    print(json.dumps(out), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--worker', help=argparse.SUPPRESS)
    ap.add_argument('--kitdir', help=argparse.SUPPRESS)
    ap.add_argument('--dd', default='', help='comma-separated tasks for the drift-diffusion check')
    a = ap.parse_args()
    if a.worker:
        worker(json.loads(a.worker), a.kitdir, a.dd == 'yes')
        return 0
    tasks = json.load(open(os.path.join(WP6, 'tasks', 'tasks.json')))['tasks']
    dd = set(filter(None, a.dd.split(',')))
    res = []
    for t in tasks:
        tmp = tempfile.mkdtemp(prefix='wp6kit_')
        try:
            for f in ('mosfet_deck.py',):
                shutil.copy2(os.path.join(WP6, 'kit', f), tmp)
            for f in ('device.json', 'mesh.json', 'zero_bias_state.npz'):
                shutil.copy2(os.path.join(WP6, 'tasks', t['task'], f), tmp)
            p = subprocess.run([PY, os.path.abspath(__file__), '--worker', json.dumps(t), '--kitdir', tmp,
                                '--dd', 'yes' if t['task'] in dd else 'no'], cwd=ROOT, env=ENV, capture_output=True, text=True)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        line = [x for x in p.stdout.splitlines() if x.startswith('{')]
        if p.returncode != 0 or not line:
            print(p.stdout[-1500:], p.stderr[-1500:])
            raise SystemExit(f"{t['task']}: validation worker failed")
        r = json.loads(line[-1])
        res.append(r)
        msg = (f"{r['task']} newton={r['newton_iterations']} max|dpsi|={r['max_dpsi_V']:.1e} V max|dn/n|={r['max_rel_dn']:.1e} "
               f"network-on-kit vs reference {r['rel_vs_reference']:+.2e}")
        if 'dd_circuit' in r:
            msg += (f" | DD circuit Richardson {r['dd_circuit']['richardson_rel']:+.2e}, contact current +-1 mV "
                    f"{r['dd_contact']['sym_1mV_rel']:+.2e}")
        print(msg, flush=True)
    with open(os.path.join(WP6, 'results', 'kit_validation.json'), 'w') as f:
        json.dump(res, f, indent=1)
    return 0


if __name__ == '__main__':
    sys.exit(main())
