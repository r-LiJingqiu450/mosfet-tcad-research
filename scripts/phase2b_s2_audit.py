#!/usr/bin/env python3
"""Production S2 Newton verification across bias and mesh settings."""

from __future__ import annotations

import csv
import hashlib
import json
import platform
import sys
import warnings
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "semiconductor_sim"))

from common import constants as C  # noqa: E402
from common.grid import make_junction_grid  # noqa: E402
from stage2_pn import solver as S  # noqa: E402
from stage2_pn.postprocess import current_density, depletion_width  # noqa: E402


TRACE_SUMMARIES = []
TRACE_ACTIVE = {}


def _f(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _trace(frame, event, arg):
    name = frame.f_code.co_name
    key = id(frame)
    if event == "call" and name == "newton_solve":
        TRACE_ACTIVE[key] = {"history": []}
        return _trace
    if key not in TRACE_ACTIVE:
        return None
    meta = TRACE_ACTIVE[key]
    loc = frame.f_locals
    if event == "line" and frame.f_lineno == 654:
        meta["history"].append({
            "iteration": int(loc.get("niter", -1)),
            "normF": _f(loc.get("normF")),
            "normF_prev": _f(loc.get("normF_prev")),
            "dmove": _f(loc.get("dmove")),
            "alpha": _f(loc.get("alpha")),
            "du_max": _f(loc.get("du_max")),
            "scale": _f(loc.get("sc")),
        })
    if event == "return":
        summary = {
            "V_a": _f(loc.get("V_a")),
            "iters": int(loc.get("niter", -1)),
            "max_iter": int(loc.get("max_iter", -1)),
            "converged": bool(loc.get("converged", False)),
            "final_normF": _f(loc.get("normF")),
            "tol": _f(loc.get("tol")),
            "history": meta["history"],
        }
        TRACE_SUMMARIES.append(summary)
        del TRACE_ACTIVE[key]
    return _trace


def _scaled_residual(result):
    Na = C.NA_CM3 * 1e6
    Nd = C.ND_CM3 * 1e6
    phi = result["phi"]
    n = result["n"]
    p = result["p"]
    Va = 0.0
    # Infer the applied bias from the left carrier boundary state when possible.
    # The caller passes the actual V_a through result['_V_a'].
    Va = float(result.get("_V_a", 0.0))
    bc = S.contact_values(Va, Na, Nd)
    phin = phi - C.VT * np.log(np.maximum(n, 1.0) / C.NI)
    phip = phi + C.VT * np.log(np.maximum(p, 1.0) / C.NI)
    phin[0], phin[-1] = bc["phi_n_L"], bc["phi_n_R"]
    phip[0], phip[-1] = bc["phi_p_L"], bc["phi_p_R"]
    x = result["_x"]
    hw = np.diff(x)
    hi = 0.5 * (hw[:-1] + hw[1:])
    Na_arr = np.where(x >= 0.0, 0.0, Na)
    Nd_arr = np.where(x >= 0.0, Nd, 0.0)
    F, *_ = S._qfp_assemble(phi, phin, phip, x, Na_arr, Nd_arr,
                             hw, hi, *S.contact_carriers(bc, Na, Nd))
    return float(np.linalg.norm(F)), float(np.max(np.abs(F)))


def _metrics(x, Va, result, previous=None, trace=None):
    Na = C.NA_CM3 * 1e6
    Nd = C.ND_CM3 * 1e6
    hw = np.diff(x)
    hi = 0.5 * (hw[:-1] + hw[1:])
    J, Jt, Jn, Jp = current_density(result["phi"], result["n"], result["p"], x)
    U = S.srh_rate(result["n"], result["p"])[0]
    qU = C.Q_E * float(np.sum(U[1:-1] * hi))
    ebal = Jn[1:] - Jn[:-1] - C.Q_E * hi * U[1:-1]
    pbal = Jp[1:] - Jp[:-1] + C.Q_E * hi * U[1:-1]
    qfp_l2, qfp_inf = _scaled_residual({**result, "_x": x, "_V_a": Va})
    Na_arr = np.where(x >= 0.0, 0.0, Na)
    Nd_arr = np.where(x >= 0.0, Nd, 0.0)
    raw_F, *_ = S._newton_system(result["phi"], result["n"], result["p"],
                                  x, Na_arr, Nd_arr, hw, hi)
    raw_phi_inf = float(np.max(np.abs(raw_F[:, 0])))
    raw_n_inf = float(np.max(np.abs(raw_F[:, 1])))
    raw_p_inf = float(np.max(np.abs(raw_F[:, 2])))
    physical_scale = max(abs(qU), abs(J), 1.0)
    identity_pass = bool(qU and abs((Jn[-1] - Jn[0]) / qU - 1.0) < 0.02
                         and abs(-(Jp[-1] - Jp[0]) / qU - 1.0) < 0.02)
    update = None if previous is None else float(np.max(np.abs(result["phi"] - previous["phi"])))
    W, Wcrit = depletion_width(result["phi"], result["n"], result["p"], Na, Nd, x)
    return {
        "V_a": float(Va), "converged": bool(result["converged"]),
        "iters": int(result["iters"]), "qfp_residual_l2": qfp_l2,
        "qfp_residual_inf": qfp_inf, "phi_update_inf_V": update,
        "raw_poisson_residual_inf_C_m2": raw_phi_inf,
        "raw_electron_residual_inf_A_m2": raw_n_inf,
        "raw_hole_residual_inf_A_m2": raw_p_inf,
        "physical_balance_scale_A_m2": float(physical_scale),
        "physical_balance_pass_1pct": bool(max(raw_n_inf, raw_p_inf) <= 1e-2 * physical_scale),
        "generation_identity_pass_2pct": identity_pass,
        "J_A_m2": float(J), "Jt_min_A_m2": float(np.min(Jt)),
        "Jt_max_A_m2": float(np.max(Jt)), "Jt_median_A_m2": float(np.median(Jt)),
        "Jt_spread_abs_A_m2": float(np.max(Jt) - np.min(Jt)),
        "Jn_delta_A_m2": float(Jn[-1] - Jn[0]),
        "Jp_delta_A_m2": float(Jp[-1] - Jp[0]), "qU_A_m2": float(qU),
        "Jn_over_qU": float((Jn[-1] - Jn[0]) / qU) if qU else None,
        "minus_Jp_over_qU": float(-(Jp[-1] - Jp[0]) / qU) if qU else None,
        "electron_balance_inf": float(np.max(np.abs(ebal))),
        "hole_balance_inf": float(np.max(np.abs(pbal))),
        "balance_scale": float(max(abs(qU), abs(J), 1e-30)),
        "electron_balance_relative": float(np.max(np.abs(ebal)) / max(abs(qU), abs(J), 1e-30)),
        "hole_balance_relative": float(np.max(np.abs(pbal)) / max(abs(qU), abs(J), 1e-30)),
        "W_charge_m": float(W), "W_criterion_m": float(Wcrit),
        "trace": trace or {},
    }


def _call(x, Va, previous=None):
    kwargs = {"tol": 1e-9}
    if previous is not None:
        kwargs.update(phi=previous["phi"], n=previous["n"], p=previous["p"])
    result = S.newton_solve(x, C.NA_CM3 * 1e6, C.ND_CM3 * 1e6, Va, **kwargs)
    trace = TRACE_SUMMARIES[-1]
    row = _metrics(x, Va, result, previous=previous, trace=trace)
    return result, row


def _sequence(x, forward_target=0.5, reverse_target=-2.0):
    rows = []
    states = {}
    eq, row = _call(x, 0.0)
    row.update(section="production", branch="equilibrium", mesh_N=len(x), termination=("converged" if row["converged"] else "budget_exhausted"))
    rows.append(row)
    states["equilibrium"] = eq
    for branch, target in (("forward", forward_target), ("reverse", reverse_target)):
        state = eq
        for Va in np.arange(0.1, abs(target) + 0.0001, 0.1) * (1.0 if target > 0 else -1.0):
            state, row = _call(x, float(Va), previous=state)
            row.update(section="production", branch=branch, mesh_N=len(x), termination=("converged" if row["converged"] else "budget_exhausted"))
            rows.append(row)
    return rows


def main():
    global TRACE_SUMMARIES
    grids = {
        "coarse": dict(L=2.0e-6, d0=20e-9, transition=0.15e-6, h_max=40e-9),
        "medium": dict(L=2.0e-6, d0=5e-9, transition=0.30e-6, h_max=20e-9),
        "production": dict(L=2.0e-6),
    }
    all_rows = []
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        sys.settrace(_trace)
        try:
            for name, config in grids.items():
                x = make_junction_grid(**config)
                rows = _sequence(x)
                for row in rows:
                    row["grid"] = name
                all_rows.extend(rows)
        finally:
            sys.settrace(None)
    warning_counts = Counter(str(item.message) for item in caught)
    payload = {
        "metadata": {
            "command": "MPLCONFIGDIR=/tmp/mplconfig .venv/bin/python scripts/phase2b_s2_audit.py",
            "python": sys.version, "platform": platform.platform(),
            "solver_sha256": hashlib.sha256((ROOT / "semiconductor_sim/stage2_pn/solver.py").read_bytes()).hexdigest(),
            "official_golden_touched": False,
        },
        "parameters": {"Na_cm3": C.NA_CM3, "Nd_cm3": C.ND_CM3,
                       "tol_V": 1e-9, "forward_target_V": 0.5, "reverse_target_V": -2.0,
                       "step_V": 0.1},
        "warnings": {"count": len(caught), "unique": dict(warning_counts)},
        "records": [{k: v for k, v in row.items() if k != "trace"} for row in all_rows],
        "trace_summary": [
            {"grid": row["grid"], "branch": row["branch"], "V_a": row["V_a"],
             "iters": row["iters"], "converged": row["converged"],
             "final_normF": row["trace"].get("final_normF"),
             "final_dmove": (row["trace"].get("history") or [{}])[-1].get("dmove"),
             "history_length": len(row["trace"].get("history", []))}
            for row in all_rows
        ],
    }
    (ROOT / "verification/phase2b_s2_verification.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    fields = ["grid", "section", "branch", "mesh_N", "V_a", "converged", "iters",
              "termination", "qfp_residual_l2", "qfp_residual_inf", "phi_update_inf_V",
              "raw_poisson_residual_inf_C_m2", "raw_electron_residual_inf_A_m2",
              "raw_hole_residual_inf_A_m2", "physical_balance_scale_A_m2",
              "physical_balance_pass_1pct", "generation_identity_pass_2pct",
              "J_A_m2", "Jt_min_A_m2", "Jt_max_A_m2", "Jt_median_A_m2",
              "Jt_spread_abs_A_m2", "Jn_delta_A_m2", "Jp_delta_A_m2", "qU_A_m2",
              "Jn_over_qU", "minus_Jp_over_qU", "electron_balance_inf", "hole_balance_inf",
              "electron_balance_relative", "hole_balance_relative", "W_charge_m", "W_criterion_m"]
    with (ROOT / "verification/phase2b_s2_metrics.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in all_rows:
            writer.writerow({field: row.get(field, "") for field in fields})
    trace_fields = ["grid", "branch", "V_a", "iteration", "normF", "normF_prev",
                    "dmove", "alpha", "du_max", "scale"]
    with (ROOT / "verification/phase2b_s2_convergence_history.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=trace_fields)
        writer.writeheader()
        for row in all_rows:
            for item in row["trace"].get("history", []):
                writer.writerow({"grid": row["grid"], "branch": row["branch"], "V_a": row["V_a"], **item})
    print(json.dumps({"rows": len(all_rows), "history_rows": sum(len(row["trace"].get("history", [])) for row in all_rows),
                      "warnings": len(caught), "nonconverged": sum(not row["converged"] for row in all_rows)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
