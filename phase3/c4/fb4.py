"""C-4 Task D l_res DIAGNOSTIC (pre-registered in phase3/specs/C4_T1v2_Preregistration_2026-10-30.md
section 4). NOT FB-1, NOT Package E, no H-C decision.

One process = one deck (model, L, level, l_res). V_D = 50 mV (continuation in 10 mV steps,
Newton tau = 1e-12 V with one extra iteration), V_G = 1.25 V by default.
  I_X,h    : 2-D drain-terminal current (A/cm) = -V1.I (drain contact + 'bd' contact if present)
  I_GCA,h  : (q mu/L) int_0^VD N_s,h^FV(phi) dphi from 1COL slices on the identical depth lines
             (Gauss-Legendre 8 points; 6-point value kept as the quadrature estimate)
  I_M0     : continuum M0 at this L
  eps_h    : I_GCA,h / I_X,h - 1  (same-mesh discrepancy, FB-1 definition)
  eps_M0   : I_M0 / I_X,h - 1
  conservation: delta_cont (3 channel sections vs drain-terminal electron current),
             delta_KCL (sum of all electron contact currents / drain-terminal current)
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from dataclasses import replace

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np

import devsim as ds

from phase3 import m0, runlog
from phase3.c1 import deck as D, measure as M
from phase3.c3.t1 import onecol_sol, g_gca, set_lateral_data, equilibrate
from phase3.c4 import deck4 as D4


def i_gca_h(geo, vg, vd, n):
    x, w = np.polynomial.legendre.leggauss(n)
    phis = 0.5 * vd * (x + 1)
    return D.Q * D.MU_N / geo.L_cm * 0.5 * vd * sum(
        ww * g_gca(geo, onecol_sol(geo, vg, phin=float(p)))[1] for p, ww in zip(phis, w))


def electron_contact_currents():
    out = {}
    for c in ds.get_contact_list(device=D.DEV):
        if 'ElectronContinuityEquation' in ds.get_contact_equation_list(device=D.DEV, contact=c):
            out[c] = ds.get_contact_current(device=D.DEV, contact=c, equation='ElectronContinuityEquation')
    return out


def state(geo, vg, model, vd):
    mode = D4.mode_of(model)
    sol0, eqrec = equilibrate(geo, vg, mode)
    its = []
    for v in np.linspace(vd / 5, vd, 5):
        if mode == 'a':
            set_lateral_data(geo, vg, sol0, onecol_sol(geo, vg, phin=float(v)))
        ds.circuit_alter(name='V1', value=float(v))
        r = D.newton(1e-12, extra=1, max_it=60)
        if not r['converged']:
            raise D.NotConverged(f'fb4 {model} {vg} {v}')
        its.append(r['iterations'])
    I = -ds.get_circuit_node_value(solution='dcop', node='V1.I')
    cc = electron_contact_currents()
    Id = cc['drain'] + cc.get('bd', 0.0)
    xs = geo.xs_ch; nch = len(xs) - 1; dx = xs[1] - xs[0]
    secs = [M.cross_section_dc(float(xs[k] + dx / 2))['I_A_per_cm']
            for k in (round(nch / 4), nch // 2, round(3 * nch / 4) - 1)]
    out = dict(I_A_per_cm=I, I_drain_terminal_electron=Id, contact_currents=cc, newton_iterations=its,
               final_update=r['final'], sections=secs,
               delta_cont=max(abs(abs(s) - abs(Id)) / abs(Id) for s in secs),
               delta_KCL=abs(sum(cc.values())) / abs(Id), delta_terminal=abs(abs(I) - abs(Id)) / abs(Id))
    ds.circuit_alter(name='V1', value=0.)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', required=True, choices=['a', 'b', 'c'])
    ap.add_argument('--L', type=float, required=True); ap.add_argument('--level', type=int, required=True)
    ap.add_argument('--lres', type=float, default=0.0); ap.add_argument('--vg', type=float, default=1.25)
    ap.add_argument('--vd', type=float, default=0.05); ap.add_argument('--out', required=True)
    a = ap.parse_args(); out = os.path.abspath(a.out); os.makedirs(out, exist_ok=True)
    geo = D4.Geo4(a.model, a.L, a.level, a.lres)
    rec = runlog.header(dict(vars(a), pm=D.pm_checklist(),
                             note='C-4 l_res diagnostic (registered); not FB-1, not Package E',
                             exploratory=(a.model == 'c' and a.lres not in D4.L_RES_FACTORS)))
    t0 = time.time()
    rec['mesh'] = D4.build(geo)
    P = replace(m0.M0Params.f0(), L=a.L * 1e-6)
    g8 = i_gca_h(geo, a.vg, a.vd, 8); g6 = i_gca_h(geo, a.vg, a.vd, 6)
    im0 = m0.reference_current(a.vg, a.vd, P)['ID'] / (P.W * 1e6) * 1e4
    s = state(geo, a.vg, a.model, a.vd)
    s.update(I_GCA_h=g8, I_GCA_h_quad_est=abs(g8 - g6), I_M0=im0,
             eps_h=g8 / s['I_A_per_cm'] - 1, eps_M0=im0 / s['I_A_per_cm'] - 1)
    rec['result'] = s
    rec['wall_s'] = time.time() - t0; rec['end_utc'] = runlog.now()
    tag = a.model + (f'_lres{a.lres:g}' if a.model == 'c' else '')
    runlog.write_json(os.path.join(out, f'fb4_{tag}_L{a.L:g}_lev{a.level}_VG{a.vg:g}.json'), rec)
    print('FB4', tag, a.L, a.level, a.vg, json.dumps(dict(eps_h=s['eps_h'], dKCL=s['delta_KCL'], dcont=s['delta_cont'])), flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
