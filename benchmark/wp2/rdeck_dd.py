#!/usr/bin/env python3
"""WP2 path P2: full drift-diffusion at small drain bias on the resolved-well deck R (pre-registration sections 5, 8).

Same mesh and equilibrium as P1. V_D = +0.1, -0.1, +1, -1 mV, each reached from V_D = 0 by Newton (step halving if
needed; tau = 1e-12 V on psi and phi_n), standard ohmic drain (s_couple = 1). Recorded per bias point:
  * DEVSIM default contact currents (drain/source/body electron) and the circuit current -V1.I;
  * weighting-function current I_W = sum_e F_e (w_n0 - w_n1), F_e = EdgeCouple * ElectronCurrent, with w the
    zero-bias drain-excitation network potential (w = 1 drain, 0 source/body); sign calibrated against the DEVSIM
    drain contact current with the indicator weight.
Symmetric differences G(|V|) = (I(+V) - I(-V)) / (2|V|) are compared with the same-mesh network G_DD.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rdeck_common as RC  # noqa: E402
from rdeck_common import D, WP, ds  # noqa: E402
from rdeck_probe import env  # noqa: E402

VERSION = 'wp2-rdeck-dd-2026-11-04'
BIASES = (1e-4, -1e-4, 1e-3, -1e-3)
TAU = 1e-12


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec='milliseconds')


def edge_flux():
    ec = np.array(ds.get_edge_model_values(device=D.DEV, region=D.SI, name='EdgeCouple'))
    j = np.array(ds.get_edge_model_values(device=D.DEV, region=D.SI, name='ElectronCurrent'))
    return ec * j


def reach(v_target):
    """Newton at V1 = v_target from the current (V1 = 0) state; halve the step on failure."""
    hist, v = [], 0.0
    h = abs(v_target)
    while abs(v_target - v) > 1e-18:
        vn = v + math.copysign(min(h, abs(v_target - v)), v_target - v)
        saved = D._state()
        ds.circuit_alter(name='V1', value=float(vn))
        r = D.newton(TAU, max_it=60, extra=2)
        hist.append(dict(target_V=vn, step_V=h, converged=bool(r['converged']), iterations=int(r['iterations']),
                         last=r['history'][-3:]))
        if r['converged']:
            v = vn
        else:
            D._restore(saved); ds.circuit_alter(name='V1', value=float(v)); h /= 2
            if h < abs(v_target) / 64:
                raise D.NotConverged(f'V_D continuation failed at {vn}')
    return hist


def currents(n0i, n1i, w, ind, sgn):
    F = edge_flux()
    cc = {}
    for c in ('source', 'drain', 'body'):
        try:
            cc[c] = float(ds.get_contact_current(device=D.DEV, contact=c, equation='ElectronContinuityEquation'))
        except ds.error:
            cc[c] = None
    return dict(I_drain=cc['drain'], I_source=cc['source'], I_body=cc['body'],
                I_circuit=float(-ds.get_circuit_node_value(solution='dcop', node='V1.I')),
                I_weighted=float(sgn * np.sum(F * (w[n0i] - w[n1i]))),
                I_indicator=float(sgn * np.sum(F * (ind[n0i] - ind[n1i]))),
                drain_edge_abs_sum=float(np.sum(np.abs(F) * np.abs(ind[n0i] - ind[n1i]))))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--id', required=True)
    ap.add_argument('--na', type=float, required=True)
    ap.add_argument('--tox', type=float, required=True)
    ap.add_argument('--L', type=float, required=True)
    ap.add_argument('--vg', type=float, required=True)
    ap.add_argument('--level', type=int, required=True)
    ap.add_argument('--body', choices=['B0', 'B1'], required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, f'P2_{a.id}_{a.body}_lev{a.level}.json')
    rec = dict(version=VERSION, path='P2', start_utc=now(), args=vars(a), argv=sys.argv, environment=env(), tau_V=TAU)
    t0 = time.time()
    _, geo, sol0, b = RC.build(a.L, a.level, a.tox, a.na, a.vg, a.body)
    rec.update(b)
    xs, ys = map(np.array, D.coords(D.SI))
    n0i, n1i, c, asym, psi, n = WP.network(geo)
    M = RC.terminal_masks(geo, xs, ys)
    T = dict(source=np.where(M['source'])[0], drain=np.where(M['drain'])[0])
    if a.body == 'B1':
        T['body'] = np.where(M['body'] & ~M['source'] & ~M['drain'])[0]
    Lap = RC.laplacian(len(xs), n0i, n1i, c)
    rd = RC.terminal_solve(Lap, n0i, n1i, c, T, 'drain')
    w = rd['phi']
    ind = np.zeros(len(xs)); ind[T['drain']] = 1.0
    rec['network_same_run'] = dict(G_DD=rd['energy'], G_SD=-rd['I']['source'], residual=rd['residual'])
    # sign calibration at equilibrium is noise-only, so calibrate at the first bias point (indicator vs DEVSIM)
    states = []
    sgn = None
    for v in BIASES:
        hist = reach(v)
        if sgn is None:
            F = edge_flux()
            raw = float(np.sum(F * (ind[n0i] - ind[n1i])))
            Idev = float(ds.get_contact_current(device=D.DEV, contact='drain', equation='ElectronContinuityEquation'))
            sgn = 1.0 if abs(raw - Idev) <= abs(raw + Idev) else -1.0
        st = dict(V_D=v, continuation=hist, **currents(n0i, n1i, w, ind, sgn))
        st['indicator_vs_devsim_abs'] = abs(st['I_indicator'] - st['I_drain'])
        states.append(st)
        ds.circuit_alter(name='V1', value=0.0)
        r0 = D.newton(TAU, max_it=60, extra=1)
        st['return_to_zero'] = dict(converged=bool(r0['converged']), iterations=int(r0['iterations']))
    rec['sign_weighted'] = sgn
    rec['states'] = states
    Gref = rd['energy']
    sym = {}
    for mag in (1e-4, 1e-3):
        sp_, sm = [s for s in states if s['V_D'] == mag][0], [s for s in states if s['V_D'] == -mag][0]
        d = {}
        for k in ('I_drain', 'I_circuit', 'I_weighted'):
            G = (sp_[k] - sm[k]) / (2 * mag)
            d[k] = dict(G=G, rel_vs_network=G / Gref - 1, one_sided_plus=sp_[k] / mag, one_sided_rel=sp_[k] / mag / Gref - 1)
        d['I_source_minus'] = -(sp_['I_source'] - sm['I_source']) / (2 * mag)
        d['G_times_V_A_per_cm'] = Gref * mag
        sym[f'{mag * 1e3:g}mV'] = d
    rec['symmetric'] = sym
    rec['indicator_check_max_rel'] = max(s['indicator_vs_devsim_abs'] / max(s['drain_edge_abs_sum'], 1e-300) for s in states)
    rec['wall_s'] = time.time() - t0
    rec['end_utc'] = now()
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(rec, f, indent=1)
    s1 = sym['1mV']
    print(f"{a.id} {a.body} lev{a.level} G_net={Gref:.6e} | 1mV: default {s1['I_drain']['rel_vs_network']:+.3e} circuit "
          f"{s1['I_circuit']['rel_vs_network']:+.3e} weighted {s1['I_weighted']['rel_vs_network']:+.3e} | 0.1mV: default "
          f"{sym['0.1mV']['I_drain']['rel_vs_network']:+.3e} weighted {sym['0.1mV']['I_weighted']['rel_vs_network']:+.3e} "
          f"| indicator check {rec['indicator_check_max_rel']:.1e} | wall {rec['wall_s']:.0f}s", flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
