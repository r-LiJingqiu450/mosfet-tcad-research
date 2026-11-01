#!/usr/bin/env python3
"""Use existing C-6 zero-bias states to quantify extraction impact.

The script performs no new device solve.  It reads the archived multi-length
states and applies the exact affine relation

    1/G(L, Vg) = (L - 2 DeltaL(Vg)) / sigma_inf(Vg).

This is an observable-impact calculation, not an independent physical model.
"""
from __future__ import annotations

import csv
import glob
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


def read_c6(root: Path):
    states = []
    for path in glob.glob(str(root / "results/phase3/C6/w1/w1_L*_lev0_VG*_hw0.5.json")):
        d = json.loads(Path(path).read_text())
        a, z = d["args"], d["decomposition"]
        if a["L"] not in (2.0, 4.0, 8.0) or a["vg"] not in (1.0, 1.25, 2.0):
            continue
        states.append(
            dict(
                source_json=str(Path(path)),
                VG_V=float(a["vg"]),
                L_um=float(a["L"]),
                G_GCA_S_per_cm=float(z["G_GCA"]),
                G_net_S_per_cm=float(z["G_net"]),
                dL_net_nm=float(z["dL_net_nm"]),
                conductance_ratio=float(z["G_net"] / z["G_GCA"]),
                conductance_correction_pct=float(100.0 * (z["G_net"] / z["G_GCA"] - 1.0)),
            )
        )
    if len(states) != 9:
        raise RuntimeError(f"expected 9 C6 states (three lengths at three Vg), found {len(states)}")
    return sorted(states, key=lambda r: (r["VG_V"], r["L_um"]))


def fit_affine(states):
    fits = []
    for vg in sorted({r["VG_V"] for r in states}):
        group = [r for r in states if r["VG_V"] == vg]
        L = np.array([r["L_um"] for r in group]) * 1e-4
        invG = np.array([1.0 / r["G_net_S_per_cm"] for r in group])
        design = np.c_[np.ones(len(L)), L]
        intercept, slope = np.linalg.lstsq(design, invG, rcond=None)[0]
        fit = design @ np.array([intercept, slope])
        dL_um = -intercept / (2.0 * slope)
        fits.append(
            dict(
                VG_V=vg,
                n_lengths=len(group),
                lengths_um=(L * 1e4).tolist(),
                DeltaL_fit_nm=float(dL_um * 1e7),
                sigma_inf_S=float(1.0 / slope),
                max_relative_fit_residual=float(np.max(np.abs((fit - invG) / invG)),),
                rms_relative_fit_residual=float(np.sqrt(np.mean(((fit - invG) / invG) ** 2))),
            )
        )
    return fits


def read_well_depth(root: Path):
    out = []
    for hw in (0.5, 1.0):
        p = root / f"results/phase3/direction_reset/well_depth/well_depth_L2_lev0_VG1.25_hw{hw:g}.json"
        d = json.loads(p.read_text())
        z = d["decomposition"]
        out.append(dict(hw_um=hw, dL_net_nm=float(z["dL_net_nm"]), G_net_S_per_cm=float(z["G_net"])))
    d0, d1 = out
    return {
        "states": out,
        "delta_dL_nm": d1["dL_net_nm"] - d0["dL_net_nm"],
        "relative_dL_change_pct_of_deeper": 100.0 * (d1["dL_net_nm"] - d0["dL_net_nm"]) / d1["dL_net_nm"],
        "relative_conductance_change_pct_of_shallow": 100.0 * (d1["G_net_S_per_cm"] - d0["G_net_S_per_cm"]) / d0["G_net_S_per_cm"],
    }


def main():
    root = Path(__file__).resolve().parents[2]
    out = root / "results/phase3/next_investment/extraction_impact"
    out.mkdir(parents=True, exist_ok=True)
    states = read_c6(root)
    fits = fit_affine(states)
    wells = read_well_depth(root)

    with (out / "state_rows.csv").open("w", newline="") as f:
        fields = list(states[0].keys())
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(states)
    with (out / "affine_fit_rows.csv").open("w", newline="") as f:
        fields = list(fits[0].keys())
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(fits)

    summary = {
        "question": "Does a bias-dependent intrinsic DeltaL create a measurable channel-length/extraction bias?",
        "equation": "1/G(L,VG) = (L - 2*DeltaL(VG))/sigma_inf(VG)",
        "source": "existing C6 W1 states; no new device solve",
        "fits": fits,
        "L2_observable_rows": [r for r in states if r["L_um"] == 2.0],
        "well_depth_comparison": wells,
        "decision": "The multi-length intercept is gate-bias dependent and the conductance impact is large enough to be a practical extraction question; this does not establish a new physical model.",
    }
    (out / "extraction_impact_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        vg = [r["VG_V"] for r in fits]
        dl = [r["DeltaL_fit_nm"] for r in fits]
        ratio = [next(r["conductance_ratio"] for r in states if r["VG_V"] == v and r["L_um"] == 2.0) for v in vg]
        fig, ax1 = plt.subplots(figsize=(6.4, 4.2))
        ax1.plot(vg, dl, "o-", label=r"$\Delta L$ fit")
        ax1.set_xlabel("V_G (V)"); ax1.set_ylabel(r"$\Delta L$ (nm)"); ax1.grid(True, alpha=0.3)
        ax2 = ax1.twinx(); ax2.plot(vg, ratio, "s--", color="tab:red", label=r"$G/G_{GCA}$ at L=2 µm")
        ax2.set_ylabel(r"$G/G_{GCA}$")
        fig.tight_layout(); fig.savefig(out / "extraction_impact.png", dpi=180); plt.close(fig)
    except Exception as exc:
        (out / "plot_warning.txt").write_text(f"plot generation failed: {exc}\n")
    print(json.dumps({"fits": fits, "well_depth_comparison": wells}, indent=2))


if __name__ == "__main__":
    main()
