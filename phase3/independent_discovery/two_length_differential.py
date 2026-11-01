#!/usr/bin/env python3
"""Evaluate the two-length zero-bias differential extraction identity.

For a fixed gate bias in the zero-bias equilibrium model,

    G(L) = sigma_inf / (L - 2 DeltaL),
    Q_e(L) = q * (N_inf * L + 2 E),

so finite differences in 1/G and Q_e remove both end terms.  This script
only post-processes independently generated endprobe JSON files; it does not
change the production solver or claim the identity outside its assumptions.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

Q_E = 1.602176634e-19  # C
MU_REFERENCE = 1350.0  # cm^2/(V s), constant-mobility deck parameter


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def channel_charge(record: dict) -> float:
    """Return q integral(N_s(x) dx) over the metallurgical channel [0, L]."""
    c = record["result"]["column_sheet_density"]
    x = np.asarray(c["x_cm"], dtype=float)
    n = np.asarray(c["N_cm2"], dtype=float)
    L = float(record["args"]["L"]) * 1e-4
    mask = (x >= -1e-14) & (x <= L + 1e-14)
    if mask.sum() < 3:
        raise ValueError(f"insufficient channel points for L={L:g} cm")
    return float(Q_E * np.trapezoid(n[mask], x[mask]))


def load(path: Path) -> dict:
    d = json.loads(path.read_text())
    if d["args"]["deck"] != "R" or d["args"]["sigma"] != 15.0:
        raise ValueError(f"unexpected configuration in {path}")
    return d


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()
    inp = args.input_dir.resolve()
    out = args.out_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    files = sorted(inp.glob("e7_R_*.json"))
    if len(files) != 6:
        raise RuntimeError(f"expected six graded probe JSON files, found {len(files)}")

    records = [load(p) for p in files]
    by_vg = {}
    for d, p in zip(records, files):
        vg = float(d["args"]["vg"])
        L = float(d["args"]["L"])
        by_vg.setdefault(vg, {})[L] = (d, p)

    rows = []
    for vg in sorted(by_vg):
        pair = by_vg[vg]
        if set(pair) != {1.0, 2.0}:
            raise RuntimeError(f"V_G={vg:g} does not have exactly L=1 and 2 um")
        d1, p1 = pair[1.0]
        d2, p2 = pair[2.0]
        L1, L2 = 1e-4, 2e-4
        G1 = float(d1["result"]["const_mu"]["G_net"])
        G2 = float(d2["result"]["const_mu"]["G_net"])
        Q1 = channel_charge(d1)
        Q2 = channel_charge(d2)
        dR_dL = (1.0 / G2 - 1.0 / G1) / (L2 - L1)
        dQ_dL = (Q2 - Q1) / (L2 - L1)
        mu = 1.0 / (dR_dL * dQ_dL)
        dl1 = float(d1["result"]["const_mu"]["dL_net_nm"])
        dl2 = float(d2["result"]["const_mu"]["dL_net_nm"])
        rows.append(
            {
                "VG_V": vg,
                "L1_um": 1.0,
                "L2_um": 2.0,
                "G_L1_S_per_cm": G1,
                "G_L2_S_per_cm": G2,
                "Qe_L1_C_per_cm": Q1,
                "Qe_L2_C_per_cm": Q2,
                "d_invG_dL_cm_per_S_per_cm": dR_dL,
                "dQe_dL_C_per_cm2": dQ_dL,
                "mu_two_length_cm2_Vs": mu,
                "mu_reference_cm2_Vs": MU_REFERENCE,
                "relative_error_vs_reference": mu / MU_REFERENCE - 1.0,
                "dL_L1_nm": dl1,
                "dL_L2_nm": dl2,
                "abs_dL_pair_difference_nm": abs(dl2 - dl1),
                "source_L1": str(p1.relative_to(Path.cwd().resolve()))
                if p1.is_relative_to(Path.cwd().resolve())
                else str(p1),
                "source_L2": str(p2.relative_to(Path.cwd().resolve()))
                if p2.is_relative_to(Path.cwd().resolve())
                else str(p2),
            }
        )

    fields = list(rows[0])
    csv_path = out / "two_length_differential_rows.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "status": "PASS_WITHIN_ZERO_BIAS_MODEL_ASSUMPTIONS",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "configuration": {
            "deck": "R",
            "N_A_cm-3": 1e18,
            "tox_nm": 2.0,
            "lateral_grading_sigma_nm": 15.0,
            "level": 0,
            "mobility": "constant mu=1350 cm^2/(V s)",
            "bias": "V_D -> 0; V_G in {1.32, 1.60, 2.00} V",
        },
        "identity": {
            "conductance": "G(L)=sigma_inf/(L-2 DeltaL)",
            "charge": "Q_e(L)=q*(N_inf*L+2E)",
            "estimator": "mu=[d(1/G)/dL]^-1/[dQ_e/dL]",
            "scope": "conditional on affine zero-bias response and the stated charge definition; not an independent physical law",
        },
        "input_sha256": {str(p.relative_to(Path.cwd())): sha256(p) for p in files},
        "rows": rows,
        "max_abs_relative_error": max(abs(r["relative_error_vs_reference"]) for r in rows),
        "max_abs_dL_pair_difference_nm": max(r["abs_dL_pair_difference_nm"] for r in rows),
        "failed_or_unseen_limits": [
            "L=0.5 um is not evaluated here because the inherited endprobe lateral mesh asserts L>0.7 um; the failed invocation is preserved in the log.",
            "No finite-drain-bias, transient, quantum, or Fermi-Dirac claim is made.",
        ],
    }
    (out / "two_length_differential_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: summary[k] for k in ("status", "git_head", "max_abs_relative_error", "max_abs_dL_pair_difference_nm")}, indent=2))


if __name__ == "__main__":
    main()
