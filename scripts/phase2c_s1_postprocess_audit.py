#!/usr/bin/env python3
"""S1 surface charge, C-V, threshold, and depletion post-processing audit."""

from __future__ import annotations

import csv
import hashlib
import json
import platform
import sys
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "semiconductor_sim"))

from common import constants as C  # noqa: E402
from common.grid import make_grid  # noqa: E402
from stage1_mos.postprocess import postprocess  # noqa: E402
from stage1_mos.run import run as s1_run  # noqa: E402
from stage1_mos.solver import _residual, solve_poisson  # noqa: E402


def _residual_arrays(phi, x):
    p0 = 0.5 * (C.NA + np.sqrt(C.NA**2 + 4.0 * C.NI**2))
    n0 = C.NI**2 / p0
    hs = np.diff(x)
    aa = 2.0 / (hs[:-1] * (hs[:-1] + hs[1:]))
    bb = -2.0 / (hs[:-1] * hs[1:])
    cc = 2.0 / (hs[1:] * (hs[:-1] + hs[1:]))
    residual = _residual(phi, p0, n0, C.NA, aa, bb, cc)
    p = p0 * np.exp(-np.clip(phi[1:-1], -3.0, 3.0) / C.VT)
    n = n0 * np.exp(np.clip(phi[1:-1], -3.0, 3.0) / C.VT)
    hi = 0.5 * (hs[:-1] + hs[1:])
    q_integral = C.Q_E * float(np.sum((p - n - C.NA) * hi))
    q_flux = C.EPS_SI * (phi[1] - phi[0]) / hs[0]
    q_right = C.EPS_SI * (phi[-1] - phi[-2]) / hs[-1]
    q_three, _ = postprocess(phi, x)
    scale = (C.Q_E / C.EPS_SI) * C.NA * np.sqrt(max(len(x) - 2, 1))
    return {
        "Q_three_C_m2": float(q_three), "Q_flux_C_m2": float(q_flux),
        "Q_integral_C_m2": float(q_integral), "Q_right_flux_C_m2": float(q_right),
        "flux_integral_rel_error": float(abs((q_flux - q_right) - q_integral) / max(abs(q_integral), 1e-300)),
        "three_vs_integral_rel_error": float(abs(q_three - q_integral) / max(abs(q_integral), 1e-300)),
        "residual_l2": float(np.linalg.norm(residual)),
        "residual_inf": float(np.max(np.abs(residual))),
        "residual_l2_relative": float(np.linalg.norm(residual) / scale),
        "residual_inf_relative": float(np.max(np.abs(residual)) / scale),
    }


def _solve_continuation(phi_s, x):
    previous = None
    history = []
    count = max(2, int(abs(phi_s) / 0.025) + 1)
    for value in np.linspace(0.0, phi_s, count):
        phi, iters, converged = solve_poisson(float(value), x, C.NA,
                                               phi_init=previous, tol=1e-11,
                                               max_iter=200)
        history.append((float(value), phi.copy(), int(iters), bool(converged)))
        previous = phi.copy()
    return history[-1][1], history


def _scan(n, n_samples=251):
    x = make_grid(1e-6, n, kind="uniform")
    phi_s = np.linspace(-0.35, 0.90, n_samples)
    rows = []
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        previous = None
        for value in phi_s:
            phi, iters, converged = solve_poisson(float(value), x, C.NA,
                                                   phi_init=previous, tol=1e-11,
                                                   max_iter=200)
            previous = phi.copy()
            m = _residual_arrays(phi, x)
            h = np.diff(x)
            q_three = m["Q_three_C_m2"]
            q_flux = m["Q_flux_C_m2"]
            q_int = m["Q_integral_C_m2"]
            rows.append({"phi_s": float(value), "converged": bool(converged),
                         "iters": int(iters), "Vg_three": float(C.VFB + value - q_three / C.COX),
                         "Vg_flux": float(C.VFB + value - q_flux / C.COX),
                         "Vg_integral": float(C.VFB + value - q_int / C.COX), **m})
    return x, rows, [str(item.message) for item in caught]


def main():
    points = {}
    grids = {}
    for n in (201, 401, 801, 1601, 3201):
        x, scan, warns = _scan(n)
        grids[str(n)] = {"N": n, "h_m": float(np.max(np.diff(x))),
                         "warnings": warns, "all_converged": all(r["converged"] for r in scan)}
        for target in (-0.5, 0.3, 0.85):
            # The C-V scan starts at -0.35 V; solve -0.5 V separately so the
            # requested accumulation point is not silently replaced by -0.35 V.
            phi_target, history_target = _solve_continuation(target, x)
            target_metrics = _residual_arrays(phi_target, x)
            q_three = target_metrics["Q_three_C_m2"]
            q_flux = target_metrics["Q_flux_C_m2"]
            q_int = target_metrics["Q_integral_C_m2"]
            grids[str(n)].setdefault("points", {})[str(target)] = {
                "phi_s": float(target), "converged": bool(history_target[-1][3]),
                "iters": int(history_target[-1][2]),
                "Vg_three": float(C.VFB + target - q_three / C.COX),
                "Vg_flux": float(C.VFB + target - q_flux / C.COX),
                "Vg_integral": float(C.VFB + target - q_int / C.COX),
                **target_metrics,
            }
        # Existing production threshold convention: nearest sampled 2 psi_F.
        psi_f = C.VT * np.log(C.NA / C.NI)
        idx = int(np.argmin(np.abs(np.array([r["phi_s"] for r in scan]) - 2.0 * psi_f)))
        grids[str(n)]["threshold"] = {"phi_s_sampled": scan[idx]["phi_s"],
                                       "Vt_three": scan[idx]["Vg_three"],
                                       "Vt_flux": scan[idx]["Vg_flux"],
                                       "Vt_integral": scan[idx]["Vg_integral"],
                                       "W_three_m": abs(scan[idx]["Q_three_C_m2"]) / (C.Q_E * C.NA),
                                       "W_flux_m": abs(scan[idx]["Q_flux_C_m2"]) / (C.Q_E * C.NA),
                                       "W_integral_m": abs(scan[idx]["Q_integral_C_m2"]) / (C.Q_E * C.NA)}
        grids[str(n)]["scan"] = scan

    # C-V sensitivity from each charge definition on the N=201 and N=801 scans.
    for n in (201, 801):
        scan = grids[str(n)]["scan"]
        ps = np.array([r["phi_s"] for r in scan])
        for method in ("three", "flux", "integral"):
            q = np.array([r[f"Q_{method if method != 'three' else 'three'}_C_m2"] for r in scan])
            cs = -np.gradient(q, ps)
            vg = np.array([r[f"Vg_{method}"] for r in scan])
            i = int(np.argmin(np.abs(vg - C.VFB)))
            grids[str(n)].setdefault("cv", {})[method] = {"Cfb_norm": float((C.COX * cs[i] / (C.COX + cs[i])) / C.COX),
                                                            "index": i}

    result = {
        "metadata": {
            "command": "MPLCONFIGDIR=/tmp/mplconfig .venv/bin/python scripts/phase2c_s1_postprocess_audit.py",
            "python": sys.version, "platform": platform.platform(),
            "solver_sha256": hashlib.sha256((ROOT / "semiconductor_sim/stage1_mos/solver.py").read_bytes()).hexdigest(),
            "postprocess_sha256": hashlib.sha256((ROOT / "semiconductor_sim/stage1_mos/postprocess.py").read_bytes()).hexdigest(),
            "official_golden_touched": False,
        },
        "parameters": {"NA_cm3": C.NA_CM3, "L_um": 1.0, "tox_nm": 20.0,
                       "tol_V": 1e-11, "n_samples": 251},
        "reference": {"psi_F_V": float(C.VT * np.log(C.NA / C.NI)),
                      "Vt_depletion_V": float(C.VFB + 2.0 * C.VT * np.log(C.NA / C.NI) +
                                              np.sqrt(4.0 * C.Q_E * C.EPS_SI * C.NA * C.VT * np.log(C.NA / C.NI)) / C.COX),
                      "Wmax_m": float(np.sqrt(4.0 * C.EPS_SI * C.VT * np.log(C.NA / C.NI) / (C.Q_E * C.NA)))},
        "grids": grids,
        "classification": {
            "Q_three": "postprocess.py three-point one-sided derivative",
            "Q_flux": "first interface finite-volume flux, with right boundary flux retained",
            "Q_integral": "interior control-volume space-charge integral",
            "threshold": "nearest sampled phi_s=2 psi_F; method-specific V_G mapping",
            "independent_check": "finite-volume flux balance, not a second physical model",
        },
    }
    (ROOT / "verification/phase2c_s1_postprocess.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    fields = ["N", "phi_s", "converged", "iters", "Vg_three", "Vg_flux", "Vg_integral",
              "Q_three_C_m2", "Q_flux_C_m2", "Q_integral_C_m2", "Q_right_flux_C_m2",
              "flux_integral_rel_error", "three_vs_integral_rel_error", "residual_l2_relative", "residual_inf_relative"]
    with (ROOT / "verification/phase2c_s1_postprocess.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for n, block in grids.items():
            for row in block["scan"]:
                writer.writerow({"N": int(n), **{field: row.get(field, "") for field in fields if field != "N"}})

    # Two plots: charge methods at selected surface potentials and threshold/C-V sensitivity.
    plot_dir = ROOT / "verification/plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    targets = (-0.5, 0.3, 0.85)
    for method, label in (("Q_three_C_m2", "three-point"), ("Q_flux_C_m2", "boundary flux"), ("Q_integral_C_m2", "volume integral")):
        ax.plot([int(n) for n in grids], [grids[n]["points"][str(targets[1])][method] for n in grids], "o-", label=label)
    ax.set_xlabel("N"); ax.set_ylabel("Q_s (C/m²) at φ_s=0.3 V"); ax.grid(True, alpha=.3); ax.legend()
    fig.tight_layout(); fig.savefig(plot_dir / "phase2c_s1_charge_methods.png", dpi=180); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ns = [int(n) for n in grids]
    for key, label in (("Vt_three", "three-point"), ("Vt_flux", "boundary flux"), ("Vt_integral", "volume integral")):
        ax.plot(ns, [grids[n]["threshold"][key] for n in grids], "o-", label=label)
    ax.axhline(result["reference"]["Vt_depletion_V"], color="k", ls="--", label="depletion reference")
    ax.set_xlabel("N"); ax.set_ylabel("V_T (V)"); ax.grid(True, alpha=.3); ax.legend()
    fig.tight_layout(); fig.savefig(plot_dir / "phase2c_s1_threshold_methods.png", dpi=180); plt.close(fig)
    print(json.dumps({"grids": list(grids), "plots": [str(plot_dir / "phase2c_s1_charge_methods.png"), str(plot_dir / "phase2c_s1_threshold_methods.png")]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
