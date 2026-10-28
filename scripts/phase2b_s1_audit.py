#!/usr/bin/env python3
"""S1 verification audit without changing the production solver."""

from __future__ import annotations

import csv
import hashlib
import json
import platform
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "semiconductor_sim"))

from common import constants as C  # noqa: E402
from common.grid import make_grid  # noqa: E402
from stage1_mos.postprocess import postprocess  # noqa: E402
from stage1_mos.run import run as s1_run  # noqa: E402
from stage1_mos.solver import _residual, solve_poisson  # noqa: E402


def _grid(n):
    return make_grid(1.0e-6, n, kind="uniform")


def _residual_metrics(phi, x):
    p0 = 0.5 * (C.NA + np.sqrt(C.NA**2 + 4.0 * C.NI**2))
    n0 = C.NI**2 / p0
    hs = np.diff(x)
    aa = 2.0 / (hs[:-1] * (hs[:-1] + hs[1:]))
    bb = -2.0 / (hs[:-1] * hs[1:])
    cc = 2.0 / (hs[1:] * (hs[:-1] + hs[1:]))
    residual = _residual(phi, p0, n0, C.NA, aa, bb, cc)
    rho = p0 * np.exp(-np.clip(phi[1:-1], -3.0, 3.0) / C.VT)
    rho -= n0 * np.exp(np.clip(phi[1:-1], -3.0, 3.0) / C.VT) + C.NA
    hi = 0.5 * (hs[:-1] + hs[1:])
    q_integral = C.Q_E * float(np.sum(rho * hi))
    surface_flux = C.EPS_SI * (phi[1] - phi[0]) / hs[0]
    right_flux = C.EPS_SI * (phi[-1] - phi[-2]) / hs[-1]
    qs_3point, _ = postprocess(phi, x)
    flux_balance = surface_flux - right_flux
    scale = (C.Q_E / C.EPS_SI) * C.NA * np.sqrt(max(len(x) - 2, 1))
    return {
        "residual_l2": float(np.linalg.norm(residual)),
        "residual_inf": float(np.max(np.abs(residual))),
        "residual_l2_relative": float(np.linalg.norm(residual) / scale),
        "residual_inf_relative": float(np.max(np.abs(residual)) / scale),
        "q_integral_interior_C_m2": float(q_integral),
        "surface_flux_C_m2": float(surface_flux),
        "right_flux_C_m2": float(right_flux),
        "flux_balance_C_m2": float(flux_balance),
        "flux_balance_rel_error": float(abs(flux_balance - q_integral) / max(abs(q_integral), 1e-300)),
        "postprocess_charge_rel_error": float(abs(qs_3point - q_integral) / max(abs(q_integral), 1e-300)),
        "postprocess_Qs_C_m2": float(qs_3point),
    }


def _solve_continuation(phi_s, x):
    previous = None
    history = []
    count = max(2, int(abs(phi_s) / 0.025) + 1)
    for value in np.linspace(0.0, phi_s, count):
        phi, iters, converged = solve_poisson(
            float(value), x, C.NA, phi_init=previous,
            tol=1e-11, max_iter=200)
        metrics = _residual_metrics(phi, x)
        history.append({"phi_s": float(value), "iters": int(iters),
                        "converged": bool(converged), **metrics})
        previous = phi.copy()
    return phi, history


def main():
    points = []
    x = _grid(801)
    psi_f = C.VT * np.log(C.NA / C.NI)
    for phi_s in (-0.5, 0.0, 0.3, 0.6, 2.0 * psi_f, 0.85):
        phi, history = _solve_continuation(phi_s, x)
        qs, es = postprocess(phi, x)
        metrics = history[-1]
        points.append({"phi_s": float(phi_s), "Qs_C_m2": float(qs),
                       "Es_V_m": float(es), **metrics,
                       "history": history})

    grid = []
    for n in (201, 401, 801, 1601, 3201):
        xx = _grid(n)
        phi, _ = _solve_continuation(2.0 * psi_f, xx)
        row = {"N": n, "h_m": float(np.max(np.diff(xx))),
               "h_min_m": float(np.min(np.diff(xx))),
               "Vg_at_2psiF": float(2.0 * psi_f - postprocess(phi, xx)[0] / C.COX),
               **_residual_metrics(phi, xx)}
        grid.append(row)

    cv = []
    for n in (201, 401, 801):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = s1_run(NGRID=n, n_samples=101, verbose=False)
        cv.append({"N": n, "n_samples": 101,
                   "Vt_numeric": float(result["Vt_numeric"]),
                   "W_num_m": float(result["W_num"]),
                   "Cfb_norm": float(result["Cfb_numeric"]),
                   "iter_max": int(result["iter_max"]),
                   "iter_mean": float(result["iter_mean"]),
                   "validation_all_passed": bool(result["validation"].all_passed()),
                   "warnings": [str(item.message) for item in caught]})

    reference = {
        "psi_F_V": float(psi_f),
        "Vt_depletion_V": float(C.VFB + 2.0 * psi_f +
                                 np.sqrt(4.0 * C.Q_E * C.EPS_SI * C.NA * psi_f) / C.COX),
        "Wmax_m": float(np.sqrt(4.0 * C.EPS_SI * psi_f / (C.Q_E * C.NA))),
    }
    result = {
        "metadata": {
            "command": "MPLCONFIGDIR=/tmp/mplconfig .venv/bin/python scripts/phase2b_s1_audit.py",
            "python": sys.version, "platform": platform.platform(),
            "solver_sha256": hashlib.sha256((ROOT / "semiconductor_sim/stage1_mos/solver.py").read_bytes()).hexdigest(),
            "postprocess_sha256": hashlib.sha256((ROOT / "semiconductor_sim/stage1_mos/postprocess.py").read_bytes()).hexdigest(),
            "official_golden_touched": False,
        },
        "parameters": {"NA_cm3": C.NA_CM3, "tox_nm": 20.0, "L_um": 1.0,
                       "tol_V": 1e-11, "points": [-0.5, 0.0, 0.3, 0.6, 2.0 * float(psi_f), 0.85]},
        "analytical_reference": reference,
        "point_checks": points,
        "grid_refinement": grid,
        "cv_grid_sensitivity": cv,
        "classification": {
            "independent_references": ["depletion approximation Vt/Wmax"],
            "discrete_identity_checks": ["finite-volume flux versus interior charge integral"],
            "postprocess_sensitivity": ["three-point surface derivative versus finite-volume flux"],
            "production_convergence_criterion": "all continuation solves converged and residual relative norm reported",
        },
    }
    (ROOT / "verification/phase2b_s1_verification.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    fields = ["section", "N", "phi_s", "Vg_at_2psiF", "Qs_C_m2", "q_integral_interior_C_m2",
              "surface_flux_C_m2", "flux_balance_rel_error", "postprocess_charge_rel_error",
              "residual_l2", "residual_inf", "residual_l2_relative", "residual_inf_relative",
              "iters", "converged", "Vt_numeric", "W_num_m", "Cfb_norm", "iter_max", "iter_mean",
              "validation_all_passed"]
    with (ROOT / "verification/phase2b_s1_metrics.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in points:
            writer.writerow({"section": "point", **{key: row.get(key, "") for key in fields}})
        for row in grid:
            writer.writerow({"section": "grid", **{key: row.get(key, "") for key in fields}})
        for row in cv:
            writer.writerow({"section": "cv", **{key: row.get(key, "") for key in fields}})
    print(json.dumps({"reference": reference,
                      "point_summary": [(r["phi_s"], r["converged"], r["residual_l2_relative"],
                                          r["flux_balance_rel_error"], r["postprocess_charge_rel_error"]) for r in points],
                      "grid_summary": [(r["N"], r["residual_l2_relative"], r["flux_balance_rel_error"]) for r in grid],
                      "cv": cv}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
