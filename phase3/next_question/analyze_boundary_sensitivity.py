"""Compare the two explicit R-II electrostatic contact specifications.

The finite-bias probe uses ``Potential = PsiLat + s_couple*V_D`` at the
resolved well contacts.  This script compares s_couple=0 and 1 at the same
states; it does not decide which is physically correct.  A small difference
is evidence of local gauge/boundary insensitivity, not a proof of equivalence.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--s0-dir", type=Path, required=True)
    ap.add_argument("--s1-dir", type=Path, required=True); ap.add_argument("--out-dir", type=Path, required=True)
    a = ap.parse_args(); out = a.out_dir.resolve(); out.mkdir(parents=True, exist_ok=True)
    rows = []
    for vd_i, vd in enumerate((0.005, 0.020)):
        # s0 files carry 1, 5, 20 mV; s1 files carry 5, 20 mV.
        j0 = vd_i + 1; j1 = vd_i
        g0 = []; g1 = []; k0 = []; k1 = []; c0 = []; c1 = []
        for L in (1, 2, 4):
            d0 = json.loads((a.s0_dir / f"r2_finite_L{L}_lev0_VG1.6_VD1_5_20mV.json").read_text())
            d1 = json.loads((a.s1_dir / f"r2_finite_L{L}_lev0_VG1.6.json").read_text())
            x0 = d0["finite_states"][j0]; x1 = d1["finite_states"][j1]
            g0.append(float(x0["G_eff_S_per_cm"])); g1.append(float(x1["G_eff_S_per_cm"]))
            k0.append(float(x0["kcl_rel"])); k1.append(float(x1["kcl_rel"]))
            c0.append(float(x0["cross_section_abs_spread_rel"])); c1.append(float(x1["cross_section_abs_spread_rel"]))
        g0 = np.array(g0); g1 = np.array(g1); L = np.array((1, 2, 4), float) * 1e-4
        rows.append({"VG_V": 1.60, "VD_V": vd,
                     "max_abs_G_relative_difference": float(np.max(np.abs(g1 / g0 - 1.0))),
                     "max_KCL_s0": max(k0), "max_KCL_s1": max(k1),
                     "max_cross_spread_s0": max(c0), "max_cross_spread_s1": max(c1),
                     "slope_s0": float(np.linalg.lstsq(np.c_[np.ones(3), L], 1/g0, rcond=None)[0][1]),
                     "slope_s1": float(np.linalg.lstsq(np.c_[np.ones(3), L], 1/g1, rcond=None)[0][1])})
    with (out / "boundary_sensitivity_rows.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    summary = {"status": "SPECIFICATION_SENSITIVITY_SMALL_IN_TESTED_WINDOW",
               "specifications": {"s0": "Potential=PsiLat", "s1": "Potential=PsiLat+V_D"},
               "rows": rows,
               "interpretation": "The observed differences are below the tested current/conservation uncertainty scale; this does not select the physical specification or prove gauge equivalence."}
    (out / "boundary_sensitivity_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
