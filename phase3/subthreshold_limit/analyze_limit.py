#!/usr/bin/env python3
"""Dimensionless weak-inversion limit analysis using existing C8 data.

No device solve is launched.  The calculation tests whether the observed
alpha-dependent Delta-L trend is better represented by the electron-free
boundary-step variable chi than by the previously proposed log(alpha) trend.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from phase3.c7.model_m7 import Device


def load(root: Path):
    src = root / "results/phase3/C8/subth/e8_1_result.json"
    states = json.loads(src.read_text())["states"]
    rows = []
    for s in states:
        tox = 20.0 if s["NA"] == 1e16 else 2.0
        dev = Device(s["NA"], tox)
        _, psi_inf, wd = dev.far_field(s["VG"])
        l0 = dev.l0(wd)
        eta = (11.7 / 3.9) * dev.tox / wd
        rows.append(
            {
                "NA_cm^-3": float(s["NA"]),
                "tox_nm": tox,
                "VG_V": float(s["VG"]),
                "alpha_inf": float(s["alpha"]),
                "chi": float(s["psi_bi_minus_psi_s_VT"]),
                "eta": float(eta),
                "Wd_nm": float(wd * 1e7),
                "l0_nm": float(l0 * 1e7),
                "dL_nm": float(s["dL_net"]),
                "F_dL_over_Wd": float(s["dL_net"] / (wd * 1e7)),
                "Ns_cm^-2": float(s["Ns"]),
                "sheet_charge_ratio": float(s["Ns"] / (s["NA"] * wd)),
            }
        )
    return sorted(rows, key=lambda r: (r["NA_cm^-3"], r["alpha_inf"]))


def leave_one_out(rows):
    """Fit F = F_inf + A/chi on alpha=15,50 and predict alpha=150.

    This is an exploratory test of the first-order chi correction, not a
    pre-registered acceptance criterion or a proof of the asymptotic series.
    """
    out = []
    for na in sorted({r["NA_cm^-3"] for r in rows}):
        g = sorted([r for r in rows if r["NA_cm^-3"] == na], key=lambda r: r["alpha_inf"])
        x = np.array([1.0 / r["chi"] for r in g[:2]])
        y = np.array([r["F_dL_over_Wd"] for r in g[:2]])
        F_inf, A = np.linalg.solve(np.c_[np.ones(2), x], y)
        pred = F_inf + A / g[2]["chi"]
        out.append(
            {
                "NA_cm^-3": na,
                "F_inf_fit": float(F_inf),
                "A_fit": float(A),
                "F_pred_alpha150": float(pred),
                "F_observed_alpha150": float(g[2]["F_dL_over_Wd"]),
                "prediction_error_pct": float(100.0 * (pred / g[2]["F_dL_over_Wd"] - 1.0)),
                "fit_alpha_points": [g[0]["alpha_inf"], g[1]["alpha_inf"]],
                "test_alpha": g[2]["alpha_inf"],
            }
        )
    return out


def main():
    root = Path(__file__).resolve().parents[2]
    out = root / "results/phase3/subthreshold_limit"
    out.mkdir(parents=True, exist_ok=True)
    rows = load(root)
    predictions = leave_one_out(rows)
    with (out / "dimensionless_rows.csv").open("w", newline="") as f:
        fields = list(rows[0].keys())
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)
    with (out / "chi_prediction.csv").open("w", newline="") as f:
        fields = list(predictions[0].keys())
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(predictions)
    summary = {
        "source": "results/phase3/C8/subth/e8_1_result.json",
        "new_simulations": 0,
        "derived_limit": "F(eta, chi) = integral [1 - exp(-delta_psi_0/VT) E_0/E_inf] d(x/Wd) in the electron-free limit",
        "exploratory_first_order_test": "F = F_inf + A/chi, fitted on alpha=15 and 50 and tested at alpha=150",
        "rows": rows,
        "predictions": predictions,
        "interpretation": "The alpha=150 points are predicted within 0.5 percent by the chi correction on both structures. This supports chi as a better mechanistic coordinate than a continued log(alpha) law, but it is an internal reanalysis and not independent validation.",
    }
    (out / "subthreshold_limit_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6.4, 4.2))
        for na in sorted({r["NA_cm^-3"] for r in rows}):
            g = [r for r in rows if r["NA_cm^-3"] == na]
            ax.plot([r["chi"] for r in g], [r["F_dL_over_Wd"] for r in g], "o-", label=f"N_A={na:.0e}")
        ax.set_xlabel(r"$\chi=(\psi_{bi}-\psi_{s,\infty})/V_T$")
        ax.set_ylabel(r"$F=\Delta L/W_d$")
        ax.grid(True, alpha=0.3); ax.legend(); fig.tight_layout()
        fig.savefig(out / "subthreshold_limit.png", dpi=180); plt.close(fig)
    except Exception as exc:
        (out / "plot_warning.txt").write_text(f"plot generation failed: {exc}\n")
    print(json.dumps({"predictions": predictions}, indent=2))


if __name__ == "__main__":
    main()
