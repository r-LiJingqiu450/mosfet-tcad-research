"""C-9 probe: phase3/c7/endprobe.py (unchanged physics and analysis) plus
  (i)  an optional Gaussian halo/pocket acceptor implant next to each metallurgical junction (deck R):
       N_halo(x, y) = N_h exp(-(x - x_h)^2 / 2 s_x^2) exp(-(y - y_h)^2 / 2 s_y^2), mirrored at the drain end,
       added to the channel acceptors (net doping = N_D - N_A - N_halo). Lombardi impurity term uses
       N_A + N_D - N_halo (approximation, documented: it under-counts total ionised impurities in the halo);
  (ii) terminal-charge outputs for a split-C-V audit: equilibrium contact (displacement) charges of all contacts
       (DEVSIM get_contact_charge, PotentialEquation) and the total electron and hole charge in the silicon
       (sum of node density x NodeVolume). The gate-to-channel charge seen by tied source+drain terminals is
       Q_SD = Q_contact(source) + Q_contact(drain) + Q_contact(oxc_s) + Q_contact(oxc_d) + Q_electrons
       (electrons are supplied by the S/D contacts; holes by the body; frozen holes follow the potential).
Zero bias, F0 physics otherwise. Output name prefix e9_."""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np

import devsim as ds

from phase3.c1 import deck as D
from phase3.c6 import wprobe as WP
from phase3.c7 import endprobe as E7

_DONOR0 = E7.donor
_LMAP0 = WP.lateral_map


def lateral_map9(L_um, a_um=0.35, h0_nm=0.5, hmax_nm=20.0):
    """C-9: allow L < 0.7 um by shrinking the graded end segment (a = 0.35 L); identical for L >= 0.7 um."""
    if L_um <= 2 * a_um:
        a_um = 0.35 * L_um
    return _LMAP0(L_um, a_um, h0_nm, hmax_nm)


WP.lateral_map = lateral_map9
HALO = dict(Nh=0.0, xh=0.0, sx=1.0, yh=0.0, sy=1.0)


def halo(geo, xs, ys):
    if HALO['Nh'] <= 0:
        return np.zeros_like(xs)
    L = geo.L_cm
    gx = lambda d: np.exp(-0.5 * ((d - HALO['xh']) / HALO['sx']) ** 2)
    gy = np.exp(-0.5 * ((ys - HALO['yh']) / HALO['sy']) ** 2)
    left = np.where(xs >= 0, gx(xs), 0.0); right = np.where(xs <= L, gx(L - xs), 0.0)
    return HALO['Nh'] * (left + right) * gy


def donor9(geo, xs, ys):
    return _DONOR0(geo, xs, ys) - halo(geo, xs, ys)


def terminal_charges():
    out = {}
    for c in ds.get_contact_list(device=D.DEV):
        try:
            out[c] = float(ds.get_contact_charge(device=D.DEV, contact=c, equation='PotentialEquation'))
        except ds.error:
            out[c] = None
    n = np.array(ds.get_node_model_values(device=D.DEV, region=D.SI, name='Electrons'))
    p = np.array(ds.get_node_model_values(device=D.DEV, region=D.SI, name='Holes'))
    vol = np.array(ds.get_node_model_values(device=D.DEV, region=D.SI, name='NodeVolume'))
    out['Q_electrons'] = float(-D.Q * np.sum(n * vol))
    out['Q_holes'] = float(D.Q * np.sum(p * vol))
    out['Q_SD'] = sum(out[c] for c in ('source', 'drain', 'oxc_s', 'oxc_d') if out.get(c) is not None) + out['Q_electrons']
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--deck', choices=['W', 'R'], required=True)
    ap.add_argument('--L', type=float, required=True); ap.add_argument('--level', type=int, default=0)
    ap.add_argument('--vg', type=float, required=True)
    ap.add_argument('--tox', type=float, default=2.0); ap.add_argument('--na', type=float, default=1e18)
    ap.add_argument('--hw', type=float, default=0.5); ap.add_argument('--lw', type=float, default=0.3)
    ap.add_argument('--ov', type=float, default=0.0); ap.add_argument('--sigma', type=float, default=0.0)
    ap.add_argument('--psibi-shift', type=float, default=0.0)
    ap.add_argument('--halo-n', type=float, default=0.0, help='halo peak (cm^-3), deck R')
    ap.add_argument('--halo-x', type=float, default=20.0, help='nm from junction'); ap.add_argument('--halo-sx', type=float, default=15.0)
    ap.add_argument('--halo-y', type=float, default=15.0, help='nm depth'); ap.add_argument('--halo-sy', type=float, default=15.0)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    HALO.update(Nh=a.halo_n, xh=a.halo_x * 1e-7, sx=a.halo_sx * 1e-7, yh=a.halo_y * 1e-7, sy=a.halo_sy * 1e-7)
    E7.donor = donor9
    name = (f"e9_{a.deck}_L{a.L:g}_lev{a.level}_VG{a.vg:g}_tox{a.tox:g}_NA{a.na:.0e}_ov{a.ov:g}_sg{a.sigma:g}"
            f"_halo{a.halo_n:.0e}.json")
    path = os.path.join(a.out, name)
    rec = dict(version='probe9-2026-10-31', start_utc=E7.now(), args=vars(a), halo=dict(HALO), N_D=E7.ND, lombardi=E7.LOMB)
    t0 = time.time()
    E7.set_bulk(a.na)
    geo = E7.Geo(a)
    rec['geometry'] = geo.describe()
    rec['mesh'] = E7.build_mesh(geo)
    D.build_physics()
    sol0, rec['equilibrium'] = E7.equilibrium(geo, a)
    rec['terminal_charges'] = terminal_charges()
    rec['result'] = E7.analyse(geo, sol0)
    rec['wall_s'] = time.time() - t0
    rec['end_utc'] = E7.now()
    json.dump(rec, open(path, 'w'), indent=1)
    c, m, tq = rec['result']['const_mu'], rec['result']['lombardi_mu'], rec['terminal_charges']
    print(f"dL {c['dL_net_nm']:.3f} | mu-dL {m['dL_net_nm']:.3f} | Q_SD {tq['Q_SD']:.4e} Qe {tq['Q_electrons']:.4e} | wall {rec['wall_s']:.0f}s")


if __name__ == '__main__':
    main()
