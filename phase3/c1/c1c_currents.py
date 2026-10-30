"""C-1(c): contact currents vs internal SG edge-flux sums, DC and AC; signs and units.

Research_05 v2.1 C-1(c). FEASIBILITY ONLY (one mesh level). Deck: C-1 deck with
1COL lateral psi data and the oxide pin next to the triple points (best C-1(a)
configuration). V_G = 1.25 V.

DC: V1 = 1 mV and 50 mV (lateral psi data kept at the V_D = 0 values, s_couple = 0,
    so this is a current-extraction test, NOT an R-I(a) physics state at finite V_D).
    Compare: drain/source contact electron currents (get_contact_current), circuit
    current V1.I, cross-section SG edge-flux sums at mid-cell lines near L/4, L/2,
    3L/4 (and next to both contacts); delta_cont, delta_KCL (R04 7.6), sign map.
AC: V1 = 0 operating point, f = 0.1 Hz and 1 MHz; linearized edge-flux sums from the
    AC node solutions vs -I_V1(f).
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
from phase3.c1 import deck as D, measure as M


def sections(L, nx):
    dx = L / (nx - 1)
    ks = [0, round((nx - 1) / 4) , (nx - 1) // 2, round(3 * (nx - 1) / 4) - 1, nx - 2]
    return [float((k + .5) * dx) for k in ks]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out', required=True); ap.add_argument('--vg', type=float, default=1.25)
    a = ap.parse_args(); out = os.path.abspath(a.out); os.makedirs(out, exist_ok=False)
    rec = runlog.header(dict(test='C-1(c)', vg=a.vg, pm=D.pm_checklist()))
    t0 = time.time()
    rec['mesh'] = D.build_mesh(corner_contacts=True); D.build_physics()
    L, nx = rec['mesh']['L_cm'], rec['mesh']['nx']
    rec['equilibrium'] = D.equilibrium_1col(a.vg)
    xs = sections(L, nx); rec['section_x_cm'] = xs
    # ---- DC
    dc = {}
    for vd in (1e-3, 5e-2):
        ds.circuit_alter(name='V1', value=vd)
        r = dict(newton=D.converge())
        cc = M.contact_currents(); r['contacts_A_per_cm'] = cc
        secs = [dict(x_cm=x, **M.cross_section_dc(x)) for x in xs]
        r['sections'] = secs
        Id = cc['drain:electron']; Is = cc['source:electron']
        r['delta_cont_vs_drain'] = max(abs(abs(s['I_A_per_cm']) - abs(Id)) / abs(Id) for s in secs)
        r['section_sign_vs_drain'] = sorted({int(np.sign(s['I_A_per_cm']) * np.sign(Id)) for s in secs})
        disp = [v for k, v in cc.items() if 'displacement' in k and v is not None]
        r['delta_KCL'] = abs(Is + Id + sum(disp)) / abs(Id)
        r['V1I_plus_drain'] = abs(cc['V1.I'] + Id) / abs(Id)
        r['V1I_minus_drain'] = abs(cc['V1.I'] - Id) / abs(Id)
        r['I_prime_A_per_um_drain'] = Id * 1e-4
        r['sign_map'] = dict(V_D_positive=vd > 0, drain_electron_current_sign=int(np.sign(Id)),
                             source_electron_current_sign=int(np.sign(Is)),
                             V1I_sign=int(np.sign(cc['V1.I'])),
                             section_sum_sign_with_s_e_rule=int(np.sign(secs[2]['I_A_per_cm'])))
        dc[f'VD={vd}'] = r
    rec['dc'] = dc
    ds.circuit_alter(name='V1', value=0.)
    rec['back_to_VD0'] = D.converge()
    # ---- AC
    ac = {}
    for f in (0.1, 1e6):
        ds.solve(type='ac', frequency=f, solver_type='direct')
        Y = -(ds.get_circuit_node_value(solution='ssac_real', node='V1.I')
              + 1j * ds.get_circuit_node_value(solution='ssac_imag', node='V1.I'))
        secs = [M.cross_section_ac(x) for x in xs]
        ac[f'f={f}'] = dict(Y_S_per_cm=[Y.real, Y.imag],
                            sections=[dict(x_cm=x, re=s.real, im=s.imag) for x, s in zip(xs, secs)],
                            re_rel_vs_ReY=[abs(abs(s.real) - abs(Y.real)) / abs(Y.real) for s in secs],
                            re_sign_vs_ReY=sorted({int(np.sign(s.real) * np.sign(Y.real)) for s in secs}),
                            im_sections=[s.imag for s in secs], ImY=Y.imag)
    rec['ac'] = ac
    rec['wall_s'] = time.time() - t0; rec['end_utc'] = runlog.now()
    runlog.write_json(os.path.join(out, 'c1c.json'), rec)
    print('C1C_DONE', json.dumps({k: dict(dcont=v['delta_cont_vs_drain'], kcl=v['delta_KCL'],
                                          v1=v['V1I_plus_drain']) for k, v in dc.items()}, default=str))
    return 0


if __name__ == '__main__':
    sys.exit(main())
