"""Bounded finite-bias probe for the resolved n+ (R-II) reference.

This is an exploratory falsification test for the zero-bias two-length
extraction result.  It does not modify the production decks.  The physical
R-II deck keeps the electrostatic potential of the resolved ohmic wells at
its equilibrium value (s_couple=0), while the drain electron quasi-Fermi
reservoir is continued to a small positive V_D.  We record terminal and
all-cross-section currents and report whether 1/(I/V_D) is affine in the
metallurgical length for a paired length test.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import os
import sys
import time
from types import SimpleNamespace

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import devsim as ds
from phase3.c1 import deck as D
from phase3.c1 import measure as M
from phase3.c7 import endprobe as EP


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="milliseconds")


def make_args(a):
    # EP.Geo and EP.equilibrium use this small attribute namespace.
    return SimpleNamespace(deck="R", L=a.L, level=a.level, tox=a.tox, na=a.na,
                           ov=a.ov, sigma=a.sigma, lw=a.lw, hw=a.hw,
                           vg=a.vg, psibi_shift=0.0)


def all_sections():
    """Evaluate every tensor-grid interior x section, not a center surrogate."""
    x, _ = map(np.array, D.coords(D.SI))
    ux = np.unique(np.round(x, 15))
    vals = []
    for xa, xb in zip(ux[:-1], ux[1:]):
        xc = 0.5 * (float(xa) + float(xb))
        rec = M.cross_section_dc(xc)
        rec["x_cm"] = xc
        vals.append(rec)
    return vals


def finite_state(geo, vd, step_count=5, tau=1e-10):
    """Continue V_D from zero and return current/conservation evidence."""
    iterations = []
    histories = []
    # A linear continuation is sufficient for the first bounded probe.  If a
    # step fails, halve it until it is below 1/128 of the target.
    v = 0.0
    target_steps = max(1, int(step_count))
    nominal = abs(vd) / target_steps
    while abs(vd - v) > 1e-18:
        h = min(nominal, abs(vd - v))
        accepted = False
        while not accepted:
            vn = v + math.copysign(h, vd - v)
            saved = D._state()
            ds.circuit_alter(name="V1", value=float(vn))
            r = D.newton(tau=tau, max_it=80, extra=1)
            histories.append({"target_V": float(vn), "step_V": float(h),
                              "converged": bool(r["converged"]),
                              "iterations": int(r["iterations"]),
                              "final": r.get("final"),
                              "history": r.get("history", [])})
            if r["converged"]:
                v = vn
                iterations.append(int(r["iterations"]))
                accepted = True
            else:
                D._restore(saved)
                ds.circuit_alter(name="V1", value=float(v))
                h *= 0.5
                if h < abs(vd) / 128.0:
                    raise RuntimeError(f"finite-bias continuation failed at {vn:g} V")

    cc = M.contact_currents()
    # The deck reports electron current with a contact sign.  Use the drain
    # contact as the terminal reference, and check source+drain KCL separately.
    idrain = float(cc["drain:electron"])
    isource = float(cc["source:electron"])
    i_circuit = float(-ds.get_circuit_node_value(solution="dcop", node="V1.I"))
    sec = all_sections()
    sec_abs = np.array([abs(float(s["I_A_per_cm"])) for s in sec], float)
    ref = max(abs(idrain), abs(i_circuit), 1e-300)
    return {
        "vd_V": float(vd),
        "I_drain_A_per_cm": idrain,
        "I_source_A_per_cm": isource,
        "I_circuit_A_per_cm": i_circuit,
        "G_eff_S_per_cm": float(i_circuit / vd),
        "kcl_rel": float(abs(isource + idrain) / ref),
        "terminal_vs_circuit_rel": float(abs(idrain - i_circuit) / ref),
        "cross_sections": sec,
        "cross_section_abs_spread_rel": float((sec_abs.max() - sec_abs.min()) / max(sec_abs.max(), 1e-300)),
        "cross_section_max_terminal_rel": float(np.max(np.abs(sec_abs - abs(idrain))) / ref),
        "newton_iterations": iterations,
        "continuation": histories,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--L", type=float, required=True)
    p.add_argument("--level", type=int, default=0)
    p.add_argument("--vg", type=float, required=True)
    p.add_argument("--vd", type=float, nargs="+", required=True)
    p.add_argument("--tox", type=float, default=2.0)
    p.add_argument("--na", type=float, default=1e18)
    p.add_argument("--sigma", type=float, default=15.0)
    p.add_argument("--ov", type=float, default=0.0)
    p.add_argument("--lw", type=float, default=0.30)
    p.add_argument("--hw", type=float, default=0.50)
    p.add_argument("--s-couple", type=float, default=0.0,
                   help="electrostatic drain-contact shift per V_D (0 or 1 are diagnostic specifications)")
    p.add_argument("--out", required=True)
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rec = {"version": "physical-r2-finite-bias-2026-11-01", "start_utc": now(),
           "args": vars(a), "model": f"R-II resolved n+ wells; s_couple={a.s_couple:g}",
           "units": {"length": "um/nm input, cm internal", "current": "A/cm", "V": "V"}}
    t0 = time.time()
    EP.set_bulk(a.na)
    aa = make_args(a)
    geo = EP.Geo(aa)
    rec["geometry"] = geo.describe()
    rec["mesh"] = EP.build_mesh(geo)
    D.build_physics()
    nd = EP.donor(geo, *map(np.array, D.coords(D.SI)))
    ds.node_solution(device=D.DEV, region=D.SI, name="NetDop")
    ds.set_node_values(device=D.DEV, region=D.SI, name="NetDop", values=list(nd - a.na))
    D._model(D.SI, "PotentialNodeCharge", "-ElectronCharge*kahan3(Holes, -Electrons, NetDop)", ("Potential", "Electrons"))
    rec["n_donor_nodes"] = int(np.sum(nd > a.na))
    sol0, eq = EP.equilibrium(geo, aa)
    rec["equilibrium"] = eq
    rec["equilibrium_psi_s_V"] = float(sol0["psi_s_V"] if isinstance(sol0, dict) and "psi_s_V" in sol0 else np.nan)
    ds.set_parameter(device=D.DEV, name="s_couple", value=float(a.s_couple))
    rec["finite_states"] = []
    for vd in a.vd:
        rec["finite_states"].append(finite_state(geo, float(vd)))
        # return to equilibrium before the next independent target
        ds.circuit_alter(name="V1", value=0.0)
        rr = D.newton(tau=1e-10, max_it=80, extra=1)
        rec["finite_states"][-1]["return_to_zero"] = {"converged": bool(rr["converged"]),
                                                          "iterations": int(rr["iterations"]),
                                                          "final": rr.get("final")}
    rec["wall_s"] = time.time() - t0
    rec["end_utc"] = now()
    fn = os.path.join(a.out, f"r2_finite_L{a.L:g}_lev{a.level}_VG{a.vg:g}.json")
    with open(fn, "w", encoding="utf-8") as f:
        json.dump(rec, f, indent=2)
    print(json.dumps({"file": fn, "wall_s": rec["wall_s"],
                      "states": [{"vd": x["vd_V"], "G": x["G_eff_S_per_cm"],
                                   "kcl": x["kcl_rel"], "cross": x["cross_section_abs_spread_rel"]}
                                  for x in rec["finite_states"]]}, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
