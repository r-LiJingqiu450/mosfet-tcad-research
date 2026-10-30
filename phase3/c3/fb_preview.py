"""Minimal finite-drain-bias DIAGNOSTIC for designing the discrepancy experiment (WS6) and
for the R-I(a)/R-I(b) reservoir-sensitivity preview (R05 C-4). NOT Package E, NOT a result.

Per (L, level, V_G), V_D = 50 mV (continuation 10 mV steps), both reservoir models:
  I_2D,h   : 2-D drain current (A/cm), R-I(a) drain psi data = 1COL(phi_n = V_D), R-I(b) Neumann.
  I_GCA,h  : (q mu/L) int_0^VD N_s,h^FV(phi) dphi from 1COL slices on the identical depth lines
             (Gauss-Legendre 8 points; 6-point value kept for the quadrature estimate).
  I_M0     : continuum M0 (phase3/m0.py) with this L.
  eps_red,h^(X) := I_GCA,h / I_2D,h^(X) - 1      (same-mesh reduction discrepancy)
  eps_red^(X)   := I_M0 / I_2D,h^(X) - 1         (continuum-reference discrepancy, R04 5.1)
  S_res(a,b)    := I_2D,h^(a) / I_2D,h^(b) - 1   (sensitivity between two specified models)
  conservation  : delta_cont (3 sections), delta_KCL; Newton attained updates.
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
from phase3.c3.t1 import Geo, onecol_sol, g_gca, set_lateral_data, equilibrate, D as _D  # noqa: F401

VD = 0.05


def i_gca_h(geo, vg, vd, n):
    x, w = np.polynomial.legendre.leggauss(n)
    phis = 0.5 * vd * (x + 1)
    return D.Q * D.MU_N / geo.L_cm * 0.5 * vd * sum(
        ww * g_gca(geo, onecol_sol(geo, vg, phin=float(p)))[1] for p, ww in zip(phis, w))


def state(geo, vg, mode, vd):
    sol0, eqrec = equilibrate(geo, vg, mode)
    its = []
    for v in np.linspace(vd / 5, vd, 5):
        if mode == 'a':
            set_lateral_data(geo, vg, sol0, onecol_sol(geo, vg, phin=float(v)))
        ds.circuit_alter(name='V1', value=float(v))
        r = D.newton(1e-12, extra=1, max_it=60)
        if not r['converged']:
            raise D.NotConverged(f'fb {mode} {vg} {v}')
        its.append(r['iterations'])
    I = -ds.get_circuit_node_value(solution='dcop', node='V1.I')
    cc = M.contact_currents()
    nx = len(geo.xs); dx = geo.xs[1] - geo.xs[0]
    secs = [M.cross_section_dc(float(geo.xs[k] + dx / 2))['I_A_per_cm']
            for k in (round((nx - 1) / 4), (nx - 1) // 2, round(3 * (nx - 1) / 4) - 1)]
    Id = cc['drain:electron']
    out = dict(I_A_per_cm=I, newton_iterations=its, final_update=r['final'],
               delta_cont=max(abs(abs(s) - abs(Id)) / abs(Id) for s in secs),
               delta_KCL=abs(cc['source:electron'] + Id) / abs(Id))
    ds.circuit_alter(name='V1', value=0.)
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--L', type=float, required=True)
    ap.add_argument('--level', type=int, required=True); ap.add_argument('--out', required=True)
    a = ap.parse_args(); out = os.path.abspath(a.out); os.makedirs(out, exist_ok=True)
    geo = Geo(a.L, a.level)
    rec = runlog.header(dict(L_um=a.L, level=a.level, VD=VD, pm=D.pm_checklist(),
                             note='diagnostic preview for experiment design; not Package E'))
    t0 = time.time()
    rec['mesh'] = D.build_mesh_explicit(geo.xs, geo.y_ox, geo.y_si, corner='tie'); D.build_physics()
    P = replace(m0.M0Params.f0(), L=a.L * 1e-6)
    rec['points'] = {}
    for vg in (0.8, 1.25, 2.0):
        r = dict(vg=vg)
        g8 = i_gca_h(geo, vg, VD, 8); g6 = i_gca_h(geo, vg, VD, 6)
        im0 = m0.reference_current(vg, VD, P)['ID'] / (P.W * 1e6) * 1e4   # A/um -> A/cm
        r.update(I_GCA_h=g8, I_GCA_h_quad_est=abs(g8 - g6), I_M0=im0)
        for mode in ('a', 'b'):
            s = state(geo, vg, mode, VD)
            s['eps_red_h'] = g8 / s['I_A_per_cm'] - 1
            s['eps_red_M0'] = im0 / s['I_A_per_cm'] - 1
            r[mode] = s
        r['S_res_ab'] = r['a']['I_A_per_cm'] / r['b']['I_A_per_cm'] - 1
        rec['points'][f'VG={vg}'] = r
        print('FB', a.L, a.level, vg, json.dumps(dict(eps_a=r['a']['eps_red_h'], eps_b=r['b']['eps_red_h'],
                                                     S=r['S_res_ab'])), flush=True)
    rec['wall_s'] = time.time() - t0; rec['end_utc'] = runlog.now()
    runlog.write_json(os.path.join(out, f'fb_L{a.L:g}_lev{a.level}.json'), rec)
    return 0


if __name__ == '__main__':
    sys.exit(main())
