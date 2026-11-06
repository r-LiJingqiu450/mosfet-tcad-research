#!/usr/bin/env python3
"""WP2 path P1: zero-bias terminal conductances of the resolved-well deck R (pre-registration sections 3-5).

For body B1 (ohmic body) the zero-bias electron network has three terminals (source, drain, body):
  G_DD = d I_D / d V_D    reference: edge energy phi^T L phi (drain excited, others 0); second form -(I_S + I_B)
  G_SD = -d I_S / d V_D   reference: -I_S (evaluated where phi = 0); second form: reciprocity -I_D(source excited)
  G_DB = -d I_B / d V_D   reference: -I_B; second form G_DD - G_SD
For body B0 (electron-blocking body, as all existing decks) G_DD = G_SD = G0 and the same-mesh Rayleigh bracket
G_rows <= G0 <= G_slice is computed.  The direct drain-side sum I_D (cancellation-prone) is recorded for the atlas.
With --smallsignal the DEVSIM coupled-Jacobian linear response (G_E equilibrated SuperLU, G_M UMFPACK, kappa_G) is
added (path P3).
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import platform
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rdeck_common as RC  # noqa: E402
from rdeck_common import D, WP, ds, g_gca  # noqa: E402

VERSION = 'wp2-rdeck-probe-2026-11-04'


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec='milliseconds')


def bracket(c, n0i, n1i, xs, ys, src):
    """Same-mesh Rayleigh bounds of the two-terminal (B0) network: G_rows <= G <= G_slice."""
    xcol = np.array(sorted(set(np.round(xs, 16))))
    ci = np.searchsorted(xcol, np.round(xs, 16))
    a_, b_ = np.minimum(ci[n0i], ci[n1i]), np.maximum(ci[n0i], ci[n1i])
    btw = (b_ - a_) == 1
    Ccol = np.bincount(a_[btw], weights=c[btw], minlength=len(xcol) - 1)[:len(xcol) - 1]
    G_slice = 1.0 / float(np.sum(1.0 / Ccol))
    lat = np.abs(ys[n0i] - ys[n1i]) < 1e-18
    G_rows = 0.0
    for y in np.array(sorted(set(np.round(ys[src], 18)))):
        m = lat & (np.abs(ys[n0i] - y) < 1e-18)
        if m.any():
            G_rows += 1.0 / float(np.sum(1.0 / c[m]))
    return G_rows, G_slice


def analyse(geo, sol0, body):
    xs, ys = map(np.array, D.coords(D.SI))
    n0i, n1i, c, asym, psi, n = WP.network(geo)
    N = len(xs)
    M = RC.terminal_masks(geo, xs, ys)
    T = dict(source=np.where(M['source'])[0], drain=np.where(M['drain'])[0])
    if body == 'B1':
        T['body'] = np.where(M['body'] & ~M['source'] & ~M['drain'])[0]
    Lap = RC.laplacian(N, n0i, n1i, c)
    rd = RC.terminal_solve(Lap, n0i, n1i, c, T, 'drain')
    rs = RC.terminal_solve(Lap, n0i, n1i, c, T, 'source')
    I_S, I_D = rd['I']['source'], rd['I']['drain']
    I_B = rd['I'].get('body', 0.0)
    G_DD, G_DD_2 = rd['energy'], -(I_S + I_B)
    G_SD, G_SD_2 = -I_S, -rs['I']['drain']
    out = dict(G_DD=G_DD, G_DD_second=G_DD_2, G_DD_rel_forms=G_DD_2 / G_DD - 1,
               G_SD=G_SD, G_SD_second=G_SD_2, G_SD_rel_forms=G_SD_2 / G_SD - 1,
               G_DD_direct_drain_sum=I_D, G_DD_direct_rel=I_D / G_DD - 1,
               residual_drain_excited=rd['residual'], residual_source_excited=rs['residual'],
               sg_equilibrium_asym=asym, n_terminal_nodes={k: int(len(v)) for k, v in T.items()})
    if body == 'B1':
        G_DB = -I_B
        out.update(G_DB=G_DB, G_DB_second=G_DD - G_SD, G_DB_rel_forms=(G_DD - G_SD) / G_DB - 1 if G_DB else None,
                   kcl_closure_rel=abs(G_DD - G_SD - G_DB) / G_DD)
    else:
        G_rows, G_slice = bracket(c, n0i, n1i, xs, ys, T['source'])
        out.update(G0=G_SD, G_rows=G_rows, G_slice=G_slice, kcl_closure_rel=abs(G_DD - G_SD) / G_DD,
                   bounds_ok=bool(G_rows * (1 - 1e-12) <= G_SD <= G_slice * (1 + 1e-12)),
                   bounds_ok_direct_drain=bool(G_rows * (1 - 1e-12) <= I_D <= G_slice * (1 + 1e-12)))
    G_GCA, Ns = g_gca(geo, sol0)
    L_nm = geo.L_cm * 1e7
    out.update(G_GCA_1COL=G_GCA, Ns_cm2=Ns, dL_SD_nm=0.5 * L_nm * (1 - G_GCA / G_SD))
    return out, dict(phi_drain=rd['phi'])


def env():
    try:
        import importlib.metadata as im
        v = {p: im.version(p) for p in ('numpy', 'scipy', 'devsim')}
    except Exception:  # pragma: no cover
        v = {}
    return dict(python=sys.version.split()[0], platform=platform.platform(), packages=v,
                DEVSIM_MATH_LIBS=os.environ.get('DEVSIM_MATH_LIBS'))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--id', required=True, help='state id, e.g. AR-1')
    ap.add_argument('--na', type=float, required=True)
    ap.add_argument('--tox', type=float, required=True)
    ap.add_argument('--L', type=float, required=True)
    ap.add_argument('--vg', type=float, required=True)
    ap.add_argument('--level', type=int, required=True)
    ap.add_argument('--body', choices=['B0', 'B1'], required=True)
    ap.add_argument('--smallsignal', action='store_true', help='also DEVSIM coupled-Jacobian linear response (P3)')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, f'P1_{a.id}_{a.body}_lev{a.level}.json')
    rec = dict(version=VERSION, path='P1', start_utc=now(), args=vars(a), argv=sys.argv, environment=env())
    t0 = time.time()
    _, geo, sol0, b = RC.build(a.L, a.level, a.tox, a.na, a.vg, a.body)
    rec.update(b)
    rec['t_equilibrium_s'] = time.time() - t0
    rec['network'], _ = analyse(geo, sol0, a.body)
    if a.smallsignal:
        t1 = time.time()
        dv = RC.devsim_smallsignal()
        ref = rec['network']['G_DD']
        dv.update(rel_G_E_vs_G_DD=dv['G_E'] / ref - 1, rel_G_M_vs_G_DD=dv['G_M'] / ref - 1, t_s=time.time() - t1)
        rec['smallsignal'] = dv
    rec['wall_s'] = time.time() - t0
    rec['end_utc'] = now()
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(rec, f, indent=1)
    nw = rec['network']
    msg = (f"{a.id} {a.body} lev{a.level} psi_hat={rec['onecol']['psi_hat']:.3f} G_DD={nw['G_DD']:.10e} "
           f"G_SD={nw['G_SD']:.10e} forms(DD,SD)=({nw['G_DD_rel_forms']:.1e},{nw['G_SD_rel_forms']:.1e}) "
           f"direct={nw['G_DD_direct_rel']:.2e} wall={rec['wall_s']:.0f}s")
    if a.body == 'B1':
        msg += f" G_DB={nw['G_DB']:.4e} kcl={nw['kcl_closure_rel']:.1e}"
    else:
        msg += f" bracket_ok={nw['bounds_ok']}"
    if a.smallsignal:
        msg += f" G_E_rel={rec['smallsignal']['rel_G_E_vs_G_DD']:.2e} kappa={rec['smallsignal']['kappa_G']:.1e}"
    print(msg, flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
