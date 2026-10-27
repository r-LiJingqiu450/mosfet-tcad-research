"""Run the reproducible S1 surface-potential and grid-convergence experiment."""
from __future__ import annotations

import csv
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common.poisson_boltzmann import PoissonBoltzmannConfig, solve_poisson_boltzmann
OUT = ROOT / "results" / "s1_baseline"

def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    potentials = np.linspace(-0.2, 0.2, 9)
    rows = []
    for points in (101, 201, 401):
        config = PoissonBoltzmannConfig(points=points)
        for surface_potential in potentials:
            solution = solve_poisson_boltzmann(float(surface_potential), config)
            rows.append({
                "points": points,
                "surface_potential_V": float(surface_potential),
                "surface_charge_C_per_m2": solution.surface_charge,
                "residual_inf": solution.residual_inf,
                "iterations": solution.iterations,
                "converged": solution.converged,
            })
    with (OUT / "results.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    np.savez_compressed(OUT / "results.npz", rows=np.array([
        [r["points"], r["surface_potential_V"], r["surface_charge_C_per_m2"], r["residual_inf"]]
        for r in rows
    ]))
    summary = {
        "experiment": "S1 intrinsic Poisson-Boltzmann MOS-capacitor baseline",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "configuration": {"temperature_K": 300.0, "intrinsic_density_m-3": 1.0e16,
                          "relative_permittivity": 11.7, "domain_debye_lengths": 20.0,
                          "surface_potentials_V": potentials.tolist(), "grid_points": [101, 201, 401]},
        "rows": len(rows),
        "all_converged": all(r["converged"] for r in rows),
        "max_residual_inf": max(r["residual_inf"] for r in rows),
        "surface_charge_odd_symmetry_check": max(
            abs(group[i]["surface_charge_C_per_m2"] + group[-1 - i]["surface_charge_C_per_m2"])
            for points in (101, 201, 401)
            for group in [[r for r in rows if r["points"] == points]]
            for i in range(len(group) // 2)
        ),
        "numpy": np.__version__,
        "matplotlib": matplotlib.__version__,
    }
    (OUT / "results.json").write_text(json.dumps(summary, indent=2) + "\n")

    fig, ax = plt.subplots(figsize=(6, 4))
    for points in (101, 201, 401):
        subset = [r for r in rows if r["points"] == points]
        ax.plot([r["surface_potential_V"] for r in subset],
                [r["surface_charge_C_per_m2"] for r in subset], marker="o", label=f"N={points}")
    ax.set(xlabel="Surface potential (V)", ylabel="Surface charge (C m$^{-2}$)",
           title="S1 Poisson–Boltzmann surface-charge relation")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "surface_charge.png", dpi=160)
    plt.close(fig)

    report = f"""# S1 baseline experiment

This run solves the intrinsic one-dimensional Poisson–Boltzmann boundary-value problem with a damped Newton method. The surface potential is prescribed at the semiconductor interface and the far boundary is neutral. The reported surface charge uses a second-order one-sided derivative.

- Commit: `{summary['git_commit']}`
- Grids: 101, 201, and 401 points over 20 Debye lengths
- Surface-potential sweep: -0.20 V to +0.20 V in 0.05 V increments
- Converged rows: {sum(r['converged'] for r in rows)}/{len(rows)}
- Maximum dimensionless residual: `{summary['max_residual_inf']:.3e}`
- Odd-symmetry absolute charge mismatch: `{summary['surface_charge_odd_symmetry_check']:.3e} C m^-2`

The grid sweep is a numerical convergence check, not an independent physical validation. The model is intrinsic and omits oxide capacitance, fixed charge, doping, mobility, channel current, CLM, and MOSFET/GCA behavior. S3 cannot be reproduced from the current repository because no S3 implementation or reference data is present.
"""
    (OUT / "report.md").write_text(report)
    manifest = {"files": sorted([p.name for p in OUT.iterdir()] + ["logs/run.log"]),
                "command": "MPLCONFIGDIR=/tmp/mplconfig .venv/bin/python experiments/run_s1_baseline.py"}
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
