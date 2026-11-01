#!/usr/bin/env python3
"""Summarise the preregistered finite n+ well-depth probe.

This script only reads the raw JSON files emitted by well_depth_probe.py.  It
does not regenerate or alter any production result or golden file.
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import math
from pathlib import Path


def load_rows(raw_dir: Path):
    rows = []
    for path in glob.glob(str(raw_dir / "well_depth_L2_lev0_VG1.25_hw*.json")):
        data = json.loads(Path(path).read_text())
        dec = data["decomposition"]
        eq = data["equilibrium"]["final"]
        rows.append(
            {
                "source_json": str(Path(path)),
                "hw_um": float(data["args"]["hw"]),
                "W_dep_um": float(data["W_dep_um"]),
                "h_over_Wdep": float(data["h_over_Wdep"]),
                "dL_net_nm": float(dec["dL_net_nm"]),
                "G_net_S_per_cm": float(dec["G_net"]),
                "net_rel_residual": float(dec["net_rel_residual"]),
                "kcl_rel": float(dec["kcl_rel"]),
                "G_GCA_S_per_cm": float(dec["G_GCA"]),
                "iterations": int(eq["iterations"]),
                "converged": bool(eq["converged"]),
                "wall_s": float(data["wall_s"]),
            }
        )
    rows.sort(key=lambda r: r["hw_um"])
    if len(rows) != 4:
        raise RuntimeError(f"expected four preregistered cases, found {len(rows)}")
    for i, row in enumerate(rows):
        if i:
            prev = rows[i - 1]
            row["successive_dL_rel_change"] = abs(row["dL_net_nm"] - prev["dL_net_nm"]) / abs(row["dL_net_nm"])
            row["successive_G_rel_change"] = abs(row["G_net_S_per_cm"] - prev["G_net_S_per_cm"]) / abs(row["G_net_S_per_cm"])
        else:
            row["successive_dL_rel_change"] = math.nan
            row["successive_G_rel_change"] = math.nan
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", type=Path, default=Path("results/phase3/direction_reset/well_depth"))
    ap.add_argument("--devsim-json", type=Path, default=Path("results/phase3/direction_reset/well_depth_devsim/well_depth_L2_lev0_VG1.25_hw1.json"))
    ap.add_argument("--out-dir", type=Path, default=Path("results/phase3/direction_reset/well_depth"))
    args = ap.parse_args()
    rows = load_rows(args.raw_dir)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    csv_path = args.out_dir / "well_depth_rows.csv"
    fields = list(rows[0].keys())
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    devsim = json.loads(args.devsim_json.read_text())
    d = devsim["decomposition"]
    ds = devsim["devsim"]
    cross = {
        "source_json": str(args.devsim_json),
        "hw_um": float(devsim["args"]["hw"]),
        "G_network_S_per_cm": float(d["G_net"]),
        "G_devsim_E_S_per_cm": float(ds["G_E"]),
        "G_devsim_M_S_per_cm": float(ds["G_M"]),
        "network_vs_devsim_E_rel": abs(d["G_net"] - ds["G_E"]) / abs(d["G_net"]),
        "network_vs_devsim_M_rel": abs(d["G_net"] - ds["G_M"]) / abs(d["G_net"]),
        "devsim_condition_number_proxy": float(ds["kappa_G"]),
        "devsim_E_dwr": float(ds["E_dwr"]),
        "devsim_E_round": float(ds["E_round"]),
    }
    changes = [r["successive_dL_rel_change"] for r in rows[1:]]
    summary = {
        "specification": {
            "NA_cm^-3": 1.0e16,
            "tox_nm": 20.0,
            "L_um": 2.0,
            "VG_V": 1.25,
            "level": 0,
            "mobility_cm2_per_Vs": 1350.0,
            "depths_um": [r["hw_um"] for r in rows],
            "W_dep_um": rows[0]["W_dep_um"],
        },
        "rows": rows,
        "cross_check": cross,
        "acceptance": {
            "criterion": "successive |dL(h)-dL(previous)|/|dL(h)| <= 0.01 once h/W_dep >= 3, and network residual < 1e-10",
            "tested_pair": "0.50 -> 1.00 um",
            "dL_relative_change_at_h_over_Wdep_3": changes[-1],
            "network_residual_at_h_over_Wdep_3": rows[-1]["net_rel_residual"],
            "saturation_pass": bool(changes[-1] <= 0.01 and rows[-1]["net_rel_residual"] < 1e-10),
        },
        "interpretation": "The fixed specification shows monotone approach to saturation but fails the preregistered 1% criterion at 3.29 W_dep; a semi-infinite replacement is therefore not established by this test.",
    }
    (args.out_dir / "well_depth_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        xs = [r["h_over_Wdep"] for r in rows]
        fig, ax1 = plt.subplots(figsize=(6.4, 4.2))
        ax1.plot(xs, [r["dL_net_nm"] for r in rows], "o-", label=r"$\Delta L_{net}$")
        ax1.set_xscale("log")
        ax1.set_xlabel(r"$h_w/W_{dep}$")
        ax1.set_ylabel(r"$\Delta L_{net}$ (nm)")
        ax1.grid(True, which="both", alpha=0.3)
        ax2 = ax1.twinx()
        ax2.plot(xs, [100 * r["successive_dL_rel_change"] if math.isfinite(r["successive_dL_rel_change"]) else math.nan for r in rows], "s--", color="tab:red", label="successive change")
        ax2.set_ylabel("successive change (%)")
        ax2.axhline(1.0, color="tab:red", alpha=0.5, linestyle=":")
        fig.tight_layout()
        fig.savefig(args.out_dir / "well_depth_saturation.png", dpi=180)
        plt.close(fig)
    except Exception as exc:  # plotting is supplementary to the machine-readable result
        (args.out_dir / "plot_warning.txt").write_text(f"plot generation failed: {exc}\n")
    print(json.dumps(summary["acceptance"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
