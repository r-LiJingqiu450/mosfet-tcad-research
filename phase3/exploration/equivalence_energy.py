#!/usr/bin/env python3
"""Run the existing W-deck probe with an independent energy identity.

The production probe reports the drain-side network sum and the source-side
sum.  This wrapper recomputes phi^T L phi from the same edge weights, which is
an independent algebraic check of the network conductance and exposes
cancellation in endpoint current extraction.  It does not alter production
modules.
"""
from __future__ import annotations
import json
import os
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np
import scipy.sparse as sp

from phase3.c6 import wprobe as WP
from phase3.c7 import endprobe as EP

ENERGY = []
_ORIG = WP.solve_net

def solve_net_energy(N, n0i, n1i, c, src, drn):
    G, Is, phi, res = _ORIG(N, n0i, n1i, c, src, drn)
    W = sp.coo_matrix(
        (np.concatenate([c, c]),
         (np.concatenate([n0i, n1i]), np.concatenate([n1i, n0i]))),
        shape=(N, N),
    ).tocsr()
    lap = (sp.diags(np.asarray(W.sum(axis=1)).ravel()) - W).tocsr()
    energy = float(phi @ (lap @ phi))
    edge_energy = float(np.sum(c * (phi[n0i] - phi[n1i]) ** 2))
    source = float(np.sum((lap @ phi)[src]))
    drain = float(np.sum((lap @ phi)[drn]))
    free = np.ones(N, dtype=bool); free[src] = False; free[drn] = False
    free_residual = float(np.max(np.abs((lap @ phi)[free])))
    ENERGY.append({"energy_G": energy, "edge_energy_G": edge_energy, "source_sum": source, "drain_sum": drain, "free_residual_max": free_residual})
    return G, Is, phi, res

WP.solve_net = solve_net_energy
if __name__ == "__main__":
    EP.main()
    out = Path(sys.argv[sys.argv.index("--out") + 1])
    try:
        vg = float(sys.argv[sys.argv.index("--vg") + 1])
    except (ValueError, IndexError):
        vg = None
    files = sorted(out.glob("e7_*.json"))
    if vg is not None:
        target = [q for q in files if f"VG{vg:g}_" in q.name]
    else:
        target = []
    if not target:
        raise RuntimeError(f"endprobe did not create the requested result JSON (vg={vg})")
    p = target[-1]
    d = json.loads(p.read_text())
    d.setdefault("energy_identity", {})["calls"] = ENERGY
    if ENERGY:
        d["energy_identity"]["constant_mobility"] = ENERGY[0]
    if len(ENERGY) > 1:
        d["energy_identity"]["lombardi_mobility"] = ENERGY[1]
    d["energy_identity"]["formula"] = "phi^T L phi, with L assembled from the same symmetric edge conductances"
    p.write_text(json.dumps(d, indent=1) + "\n")
    print(json.dumps(d["energy_identity"], indent=2))
