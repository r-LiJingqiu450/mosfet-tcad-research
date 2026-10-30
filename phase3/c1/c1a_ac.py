"""C-1(a): small-signal AC with a circuit voltage source on a custom (R-I(a)-type) contact.

Research_05 v2.1 C-1(a), with L1 (direct linear response from the exported matrix,
validated by C-1(b)) as the cross-check. FEASIBILITY ONLY: one mesh level, not the
D-1 mesh family, no uncertainty analysis. The zero-bias identity ratio reported
here is a PREVIEW of the T1 quantity on this deck, NOT a T1 result.

Per V_G in {0.8, 1.25, 2.0} V at V_D = 0 (circuit source V1 = 0, acreal = 1 V):
  L2: Y(f) = -I_V1(f) / 1 V from solve(type='ac'), f = 1e6 ... 1e-3 Hz;
      E_freq(f1, f1/10) = |ReY(f1) - ReY(f1/10)| / (0.99 G).
  L1: G = -(dI_V1/dV1) from A^-1 (rhs(V1=delta) - rhs(V1=0)) with the exported
      static matrix; also DEVSIM's own first Newton step / delta.
  Decoupling: repeat L1 and L2 with s_couple = 1 (lateral psi data follow V_D).
  Identity preview: r0 = G / G_GCA,h - 1 with G_GCA,h = (q mu/L) sum n0_i w_i.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
import scipy.sparse.linalg as spla

import devsim as ds

from phase3 import m0, runlog
from phase3.c1 import deck as D, measure as M
from phase3.c1.c1b_jacobian import REC, recording_callback, export, snapshot, restore

FREQS = [1e6, 1e4, 1e2, 10., 1., .1, .01, 1e-3]


def circuit_rows():
    """Equation rows of circuit unknowns (valid only during a solve -> callback)."""
    stock = ds.get_parameter(name='solver_callback')
    REC.clear(); X0 = snapshot()
    ds.set_parameter(name='direct_solver', value='custom')
    ds.set_parameter(name='solver_callback', value=recording_callback)
    try:
        D.solve(abs_err=1e300, rel_err=1e300, it=1)
    finally:
        ds.set_parameter(name='solver_callback', value=stock)
    restore(X0)
    return {k[1]: int(v[0]) for k, v in REC['eqnums'].items() if k[0] == 'circuit'}


def l1(rows, deltas=(1e-3, 1e-6)):
    out = []
    A, b0, *_ = export()
    lu = spla.splu(A.tocsc())
    for d in deltas:
        ds.circuit_alter(name='V1', value=d)
        _, bd, *_ = export()
        s = lu.solve(bd - b0)
        g_ext = -s[rows['V1.I']] / d
        # DEVSIM's own first undamped Newton step from the V1=0 solution
        X0 = snapshot()
        D.solve(abs_err=1e300, rel_err=1e300, it=1)
        X1 = snapshot()
        g_nwt = -(X1[('circuit', 'V1.I')][0] - X0[('circuit', 'V1.I')][0]) / d
        restore(X0)
        ds.circuit_alter(name='V1', value=0.)
        out.append(dict(delta_V=d, G_ext_S_per_cm=float(g_ext), G_newton_step_S_per_cm=float(g_nwt),
                        rhs_linear_in_delta=float(np.max(np.abs((bd - b0) / d)))))
    return out


def l2():
    rows = []
    for f in FREQS:
        ds.solve(type='ac', frequency=f, solver_type='direct')
        ir = ds.get_circuit_node_value(solution='ssac_real', node='V1.I')
        ii = ds.get_circuit_node_value(solution='ssac_imag', node='V1.I')
        vr = ds.get_circuit_node_value(solution='ssac_real', node='drain_bias')
        vi = ds.get_circuit_node_value(solution='ssac_imag', node='drain_bias')
        rows.append(dict(f_Hz=f, ReY_S_per_cm=-ir, ImY_S_per_cm=-ii, V_drain_ac=[vr, vi]))
    return rows


def efreq(rows, G):
    by = {r['f_Hz']: r['ReY_S_per_cm'] for r in rows}
    out = []
    for f1 in FREQS:
        f2 = f1 / 10
        if f2 in by:
            out.append(dict(f1=f1, f2=f2, E_freq=abs(by[f1] - by[f2]) / (0.99 * G)))
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out', required=True)
    ap.add_argument('--lateral', choices=['neumann2d', 'onecol'], default='onecol',
                    help='source of lateral psi data: 2-D Neumann run (attempt 1) or independent 1COL FV')
    ap.add_argument('--oxide-corner-pin', type=int, default=0,
                    help='1: pin the oxide node adjacent to each Si/SiO2/contact triple point to 1COL data')
    a = ap.parse_args(); out = os.path.abspath(a.out); os.makedirs(out, exist_ok=False)
    rec = runlog.header(dict(test='C-1(a)', freqs=FREQS, lateral=a.lateral, oxide_corner_pin=a.oxide_corner_pin,
                             pm=D.pm_checklist()))
    t0 = time.time()
    rec['mesh'] = D.build_mesh(corner_contacts=bool(a.oxide_corner_pin)); D.build_physics()
    if a.oxide_corner_pin:
        x, y = map(np.array, D.coords(D.OX))
        pinned = {c: sorted(set(n for e in ds.get_element_node_list(device=D.DEV, region=D.OX, contact=c) for n in e))
                  for c in ('oxc_s', 'oxc_d')}
        rec['oxide_pin_nodes_from_elements'] = pinned
    L = rec['mesh']['L_cm']
    P = m0.M0Params.f0()
    rec['points'] = {}
    for vg in (0.8, 1.25, 2.0):
        r = dict(vg=vg)
        ds.set_parameter(device=D.DEV, name='s_couple', value=0.)
        r['equilibrium'] = D.equilibrium(vg) if a.lateral == 'neumann2d' else D.equilibrium_1col(vg)
        r['lateral_uniformity'] = [M.lateral_uniformity('Potential'), M.lateral_uniformity('Electrons')]
        x, _ = M.si_arrays(); xs = np.unique(np.round(x, 15))
        r['G_GCA_h'] = M.gca_conductance(L)
        r['G_GCA_h_by_column'] = [M.gca_conductance(L, xc)['G_S_per_cm'] for xc in xs]
        G = r['G_GCA_h']['G_S_per_cm']
        rows = circuit_rows(); r['circuit_rows'] = rows
        r['L1'] = l1(rows)
        r['L2'] = l2()
        Gl1 = r['L1'][0]['G_ext_S_per_cm']
        r['E_freq'] = efreq(r['L2'], Gl1)
        r['L2_vs_L1_rel'] = [dict(f_Hz=q['f_Hz'], rel=q['ReY_S_per_cm'] / Gl1 - 1) for q in r['L2']]
        ReY01 = next(q['ReY_S_per_cm'] for q in r['L2'] if q['f_Hz'] == .1)
        r['identity_preview_NOT_T1'] = dict(r0_L1=Gl1 / G - 1, r0_L2_at_0p1Hz=ReY01 / G - 1,
                                            r0_newton_step=r['L1'][0]['G_newton_step_S_per_cm'] / G - 1)
        # decoupling: lateral psi data follow V_D with unit sensitivity
        ds.set_parameter(device=D.DEV, name='s_couple', value=1.)
        D.converge()
        l1c = l1(rows); l2c = l2()
        r['decoupling_s1'] = dict(
            G_L1_S_per_cm=l1c[0]['G_ext_S_per_cm'],
            rel_change_L1=l1c[0]['G_ext_S_per_cm'] / Gl1 - 1,
            ReY_0p1Hz=next(q['ReY_S_per_cm'] for q in l2c if q['f_Hz'] == .1),
            rel_change_L2_0p1Hz=next(q['ReY_S_per_cm'] for q in l2c if q['f_Hz'] == .1) / ReY01 - 1,
            ImY_1e6Hz_s0=next(q['ImY_S_per_cm'] for q in r['L2'] if q['f_Hz'] == 1e6),
            ImY_1e6Hz_s1=next(q['ImY_S_per_cm'] for q in l2c if q['f_Hz'] == 1e6))
        ds.set_parameter(device=D.DEV, name='s_couple', value=0.)
        D.converge()
        # context only (different discretisation; not T2): M0 G' in S/um
        r['context_M0_G0_S_per_um'] = m0.conductance0(vg, P)
        r['G_GCA_h_S_per_um'] = G * 1e-4
        rec['points'][f'VG={vg}'] = r
    rec['wall_s'] = time.time() - t0; rec['end_utc'] = runlog.now()
    runlog.write_json(os.path.join(out, 'c1a.json'), rec)
    summ = {k: dict(G_L1=v['L1'][0]['G_ext_S_per_cm'], preview=v['identity_preview_NOT_T1'],
                    Efreq=[(e['f1'], e['E_freq']) for e in v['E_freq']],
                    decoupling=v['decoupling_s1']['rel_change_L1']) for k, v in rec['points'].items()}
    print('C1A_DONE', json.dumps(summ, default=str))
    return 0


if __name__ == '__main__':
    sys.exit(main())
