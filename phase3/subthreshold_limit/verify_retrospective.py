#!/usr/bin/env python3
"""Verify the C10 retrospective fixed-coefficient chi fit.

This reads the recorded C10 ``all_states`` table and the coefficients produced
by the exploratory Phase 3 fit.  It does not run a device simulation, refit
coefficients, or modify any historical result.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


EPS_SI_OVER_EPS_OX = 11.7 / 3.9


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--c10",
        type=Path,
        default=Path("results/phase3/C10/analysis/c10_result.json"),
    )
    ap.add_argument(
        "--coefficients",
        type=Path,
        default=Path("results/phase3/subthreshold_limit/subthreshold_limit_summary.json"),
    )
    ap.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/phase3/subthreshold_limit/retrospective"),
    )
    ap.add_argument(
        "--attachment",
        type=Path,
        help="Optional downloaded retrospective JSON; recorded for provenance only.",
    )
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    c10 = json.loads(args.c10.read_text())
    coeff_doc = json.loads(args.coefficients.read_text())
    coeffs = {
        float(row["NA_cm^-3"]): row
        for row in coeff_doc["predictions"]
    }
    selected = []
    for state in c10["all_states"]:
        if (
            not state["file"].startswith("results/phase3/C10/")
            or state["level"] != 0
            or abs(state["hw"] - 0.5) > 1e-12
            or abs(state["shift"]) > 1e-12
        ):
            continue
        if state["na"] not in coeffs or state["tox"] not in (2.0, 20.0):
            continue
        # Keep the corrected source-side value.  The drain-side quantity is
        # explicitly invalid in the deep-subthreshold cancellation cases.
        wd_nm = float(state["Wd_nm"])
        chi = float(state["psibi_hat_eff"] - state["psi_hat"])
        eta = EPS_SI_OVER_EPS_OX * float(state["tox"]) / wd_nm
        f_obs = float(state["dL"]) / wd_nm
        coeff = coeffs[state["na"]]
        f_pred = float(coeff["F_inf_fit"]) + float(coeff["A_fit"]) / chi
        rel_pct = 100.0 * (f_pred / f_obs - 1.0)
        selected.append(
            {
                "source_file": state["file"],
                "NA_cm3": state["na"],
                "tox_nm": state["tox"],
                "L_um": state["L"],
                "VG_V": state["vg"],
                "level": state["level"],
                "hw_requested_um": state["hw"],
                "hw_actual_um": state["hw_real_um"],
                "alpha_inf": state["alpha"],
                "chi": chi,
                "eta": eta,
                "Wd_nm": wd_nm,
                "dL_source_corrected_nm": state["dL"],
                "F_observed": f_obs,
                "F_predicted_fixed_coefficients": f_pred,
                "absolute_error_F": f_pred - f_obs,
                "relative_error_percent": rel_pct,
                "source_current_within_bounds": state["bounds_ok_src"],
                "drain_current_within_bounds": state["bounds_ok_drain"],
            }
        )

    selected.sort(key=lambda r: (r["NA_cm3"], r["VG_V"]))
    if len(selected) != 10:
        raise RuntimeError(f"expected 10 selected C10 states, found {len(selected)}")

    fields = list(selected[0])
    with (args.output_dir / "retrospective_check.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(selected)

    result = {
        "status": "RETROSPECTIVE_REPRODUCTION",
        "new_device_simulations": 0,
        "formula": "F_prediction=F_inf+A/chi; chi=psibi_hat_eff-psi_hat; F_observed=dL_source/Wd",
        "constants": {"eps_si_over_eps_ox": EPS_SI_OVER_EPS_OX},
        "inputs": {
            "c10_path": str(args.c10),
            "c10_sha256": sha256(args.c10),
            "coefficients_path": str(args.coefficients),
            "coefficients_sha256": sha256(args.coefficients),
            "attachment_name": "chi_fit_C10_retrospective_check.json",
            "attachment_sha256": sha256(args.attachment) if args.attachment else None,
        },
        "selection": "C10 all_states, level 0, requested hw=0.5 um, NA/tox matching the two fitted devices; no refitting",
        "rows": selected,
        "summary": {
            "n": len(selected),
            "max_abs_relative_error_percent": max(abs(r["relative_error_percent"]) for r in selected),
            "min_relative_error_percent": min(r["relative_error_percent"] for r in selected),
            "max_relative_error_percent": max(r["relative_error_percent"] for r in selected),
            "source_bound_failures": sum(not r["source_current_within_bounds"] for r in selected),
            "drain_bound_failures": sum(not r["drain_current_within_bounds"] for r in selected),
        },
        "caveats": [
            "This is a retrospective check on historical C10 data, not a new or blind experiment.",
            "Eta changes with gate bias; failure of fixed coefficients does not disprove a more general F(eta, chi).",
            "The calculation uses the cancellation-free source-side end correction; invalid drain-side values are not substituted.",
        ],
    }
    (args.output_dir / "retrospective_check.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result["summary"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
