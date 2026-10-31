"""TC-1: pre-registered mechanism test of Proposition 2 (phase3/reports/
C4_Small_Bias_Perturbation_Analysis_2026-10-30.md section 4). NOT FB-1, NOT Package E.

One process = one (L, level); all four end treatments aa, bb, ab, ba are run on the SAME
corrected R-I(b) mesh (deck4 model 'b': dual bottom corners). Lateral psi:
  Dirichlet end : contact PotentialEquation psi = PsiLat (1COL slice at the local phi_n),
                  plus the triple-point oxide tie oxc_* with the same value (Research_08 sec. 5)
  Neumann end   : no psi equation on the lateral contact; oxc_* inert
V_G = 1.25 V; V_D in {0.5, 1, 2, 4, 8} mV; Newton tau = 1e-13 V and a repeat at 1e-15 V;
GCA side: 4-point Gauss-Legendre over 1COL slices (3-point kept as quadrature estimate).
Fit r(d) = c1 d + c2 d^2 + c3 d^3, weights 1/e_i^2,
  e_i = |r(tau) - r(tau/100)| + |I_GCA4 - I_GCA3| / I_GCA4 + 1e-14,
covariance scaled by max(1, chi2/dof).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np

import devsim as ds

from phase3 import runlog
from phase3.c1 import deck as D
from phase3.c3.t1 import equilibrate, onecol_sol, set_lateral_data, i_gca
from phase3.c4 import deck4 as D4

DELTAS = [0.5e-3, 1e-3, 2e-3, 4e-3, 8e-3]
ENDS = {'aa': (True, True), 'bb': (False, False), 'ab': (True, False), 'ba': (False, True)}


def set_ends(src_dirichlet, drn_dirichlet):
    for c, oc, b, on in (('source', 'oxc_s', 'source_bias', src_dirichlet), ('drain', 'oxc_d', 'drain_bias', drn_dirichlet)):
        for cc in (c, oc):
            try:
                ds.delete_contact_equation(device=D.DEV, contact=cc, name='PotentialEquation')
            except ds.error:
                pass
        if not on:
            continue
        for cc in (c, oc):
            ds.contact_node_model(device=D.DEV, contact=cc, name=f'{cc}_psi', equation=f'Potential - (PsiLat + s_couple*{b})')
            ds.contact_node_model(device=D.DEV, contact=cc, name=f'{cc}_psi:Potential', equation='1')
            kw = dict(device=D.DEV, contact=cc, name='PotentialEquation', node_model=f'{cc}_psi',
                      edge_charge_model='contactcharge_edge')
            if c == 'drain':
                ds.contact_node_model(device=D.DEV, contact=cc, name=f'{cc}_psi:{b}', equation='-s_couple')
                kw['circuit_node'] = b
            ds.contact_equation(**kw)


def fit(ds_, rs, es):
    X = np.vstack([ds_, np.square(ds_), np.power(ds_, 3)]).T
    w = 1 / np.square(es)
    Wm = np.diag(w)
    C = np.linalg.inv(X.T @ Wm @ X)
    c = C @ X.T @ Wm @ rs
    res = rs - X @ c
    chi2 = float(np.sum(w * res ** 2)); dof = len(rs) - 3
    s = max(1.0, chi2 / dof) if dof > 0 else 1.0
    sig = np.sqrt(np.diag(C) * s)
    return dict(c1=float(c[0]), c2=float(c[1]), c3=float(c[2]), s1=float(sig[0]), s2=float(sig[1]), s3=float(sig[2]),
                chi2=chi2, dof=dof, c1_over_s1=float(abs(c[0]) / sig[0]))


def run_ends(geo, vg, name, sol0, ig):
    src, drn = ENDS[name]
    ds.circuit_alter(name='V1', value=0.)
    set_lateral_data(geo, vg, sol0)
    set_ends(src, drn)
    D.converge()
    rows = []
    for d in DELTAS:
        set_lateral_data(geo, vg, sol0, onecol_sol(geo, vg, phin=d))
        ds.circuit_alter(name='V1', value=d)
        r1 = D.newton(1e-13, extra=0, max_it=60)
        if not r1['converged']:
            raise D.NotConverged(f'tc1 {name} {d}')
        I1 = -ds.get_circuit_node_value(solution='dcop', node='V1.I')
        r2 = D.newton(1e-15, extra=0, max_it=20)
        I2 = -ds.get_circuit_node_value(solution='dcop', node='V1.I')
        g4, g3 = ig[d]
        r = g4 / I1 - 1                      # epsilon = I_GCA / I_2D - 1 (FB-1 sign convention)
        e = abs(g4 / I2 - 1 - r) + abs(g4 - g3) / abs(g4) + 1e-14
        rows.append(dict(delta_V=d, I2D=I1, I2D_tight=I2, IGCA4=g4, IGCA3=g3, eps=r, e=e,
                         newton=[r1['iterations'], r2['converged']]))
    ds.circuit_alter(name='V1', value=0.)
    set_lateral_data(geo, vg, sol0)
    D.converge()
    f = fit(np.array([x['delta_V'] for x in rows]), np.array([x['eps'] for x in rows]), np.array([x['e'] for x in rows]))
    return dict(points=rows, fit=f)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=float, required=True); ap.add_argument('--level', type=int, required=True)
    ap.add_argument('--vg', type=float, default=1.25); ap.add_argument('--out', required=True)
    a = ap.parse_args(); out = os.path.abspath(a.out); os.makedirs(out, exist_ok=True)
    geo = D4.Geo4('b', a.L, a.level)
    rec = runlog.header(dict(vars(a), deltas=DELTAS, ends=ENDS, pm=D.pm_checklist(),
                             note='TC-1 mechanism test (registered); not FB-1'))
    t0 = time.time()
    rec['mesh'] = D4.build(geo)
    ds.set_parameter(device=D.DEV, name='s_couple', value=0.)
    sol0, rec['equilibrium'] = equilibrate(geo, a.vg, 'b')
    ig = {d: (i_gca(geo, a.vg, d, 4), i_gca(geo, a.vg, d, 3)) for d in DELTAS}
    rec['ends'] = {}
    for name in ('bb', 'aa', 'ab', 'ba'):
        rec['ends'][name] = run_ends(geo, a.vg, name, sol0, ig)
        f = rec['ends'][name]['fit']
        print('TC1', a.L, a.level, name, json.dumps({k: f[k] for k in ('c1', 's1', 'c2', 's2', 'c1_over_s1')}), flush=True)
    rec['wall_s'] = time.time() - t0; rec['end_utc'] = runlog.now()
    runlog.write_json(os.path.join(out, f'tc1_L{a.L:g}_lev{a.level}.json'), rec)
    return 0


if __name__ == '__main__':
    sys.exit(main())
