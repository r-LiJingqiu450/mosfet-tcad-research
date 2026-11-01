"""Fixed-depth R-II probe for the direction-reset hypothesis.

The production endprobe snaps the well depth to the nearest existing depth
line.  This isolated probe inserts the requested depth as an exact mesh line,
so a depth trend is not confused with that snapping artifact.  It reuses the
verified F0/R-II equations and changes no production module.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import sys
import time
from types import SimpleNamespace

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

import devsim as ds
from phase3 import meshmap as MM
from phase3.c1 import deck as D
from phase3.c3.t1 import onecol_sol
from phase3.c6 import wprobe as WP
from phase3.c7 import endprobe as EP


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="milliseconds")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--L", type=float, default=2.0); p.add_argument("--level", type=int, default=0)
    p.add_argument("--vg", type=float, default=1.25); p.add_argument("--hw", type=float, required=True)
    p.add_argument("--na", type=float, default=1e16); p.add_argument("--tox", type=float, default=20.0)
    p.add_argument("--sigma", type=float, default=0.0); p.add_argument("--lw", type=float, default=0.30)
    p.add_argument("--out", required=True); p.add_argument("--no-devsim", action="store_true")
    a = p.parse_args(); os.makedirs(a.out, exist_ok=True); t0 = time.time()
    EP.set_bulk(a.na)
    aa = SimpleNamespace(deck="R", L=a.L, level=a.level, tox=a.tox, na=a.na, ov=0.0,
                         sigma=a.sigma, lw=a.lw, hw=a.hw, vg=a.vg, psibi_shift=0.0)
    geo = EP.Geo(aa)
    h_exact = a.hw * 1e-4
    geo.y_si = np.unique(np.concatenate([geo.y_si, [h_exact]]))
    geo.jw_si = int(np.flatnonzero(np.isclose(geo.y_si, h_exact, rtol=0.0, atol=1e-18))[0])
    geo.h_w_cm = float(geo.y_si[geo.jw_si])
    rec = dict(version="direction-reset-well-depth-2026-11-01", start_utc=now(), args=vars(a),
               hypothesis="fixed-depth R-II sensitivity; exact depth breakpoint", geometry=geo.describe(),
               W_dep_um=float(MM.wdep_max_cm(a.na) * 1e4), h_over_Wdep=float(a.hw / (MM.wdep_max_cm(a.na) * 1e4)))
    rec["mesh"] = EP.build_mesh(geo); D.build_physics()
    xs, ys = map(np.array, D.coords(D.SI)); nd = EP.donor(geo, xs, ys)
    ds.node_solution(device=D.DEV, region=D.SI, name="NetDop")
    ds.set_node_values(device=D.DEV, region=D.SI, name="NetDop", values=list(nd - a.na))
    D._model(D.SI, "PotentialNodeCharge", "-ElectronCharge*kahan3(Holes, -Electrons, NetDop)", ("Potential", "Electrons"))
    sol0, eq = EP.equilibrium(geo, aa); rec["equilibrium"] = eq
    WP.contact_mask = lambda region, g, x, y: EP.contact_mask(region, geo, x, y)
    rec["decomposition"] = WP.decomposition(geo, sol0)
    rec["devsim"] = None if a.no_devsim else WP.devsim_G()
    rec["wall_s"] = time.time() - t0; rec["end_utc"] = now()
    fn = os.path.join(a.out, f"well_depth_L{a.L:g}_lev{a.level}_VG{a.vg:g}_hw{a.hw:g}.json")
    with open(fn, "w", encoding="utf-8") as f: json.dump(rec, f, indent=2)
    print(json.dumps({"file": fn, "hw_um": a.hw, "h_over_Wdep": rec["h_over_Wdep"],
                      "dL_net_nm": rec["decomposition"]["dL_net_nm"],
                      "G_net": rec["decomposition"]["G_net"], "wall_s": rec["wall_s"]}, indent=2), flush=True)


if __name__ == "__main__":
    raise SystemExit(main())
