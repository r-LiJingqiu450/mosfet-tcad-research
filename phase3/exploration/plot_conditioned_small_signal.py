#!/usr/bin/env python3
"""Plot the already-generated conditioning audit CSVs; no solver is run."""
from __future__ import annotations
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/phase3/exploration"
with (OUT / "conditioned_small_signal_rows.csv").open(newline="") as f:
    rows = list(csv.DictReader(f))
with (OUT / "coupled_small_signal_rows.csv").open(newline="") as f:
    coupled = list(csv.DictReader(f))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.1), constrained_layout=True)
for r in rows:
    ax1.scatter(float(r["alpha"]), float(r["abs_err_lin"])*100, c="#2166ac", s=20)
ax1.set_xscale("log")
ax1.set_xlabel(r"$\alpha$ (C10 dimensionless state parameter)")
ax1.set_ylabel(r"$|\Delta L_{lin}/\Delta L_{2D}-1|$ (%)")
ax1.set_title("H-lin error across 24 C11 states")
ax1.grid(True, which="both", alpha=0.25)

vgs = [float(r["vg"]) for r in coupled]
net = [float(r["network_vs_independent_rel_dL"])*100 for r in coupled]
ss = [float(r["abs_G_M_minus_G_E_over_max"])*100 for r in coupled]
width = 0.18
x = list(range(len(vgs)))
ax2.bar([i-width/2 for i in x], net, width=width, label="network vs FV dL", color="#1b9e77")
ax2.bar([i+width/2 for i in x], ss, width=width, label="coupled G_M vs G_E", color="#d95f02")
ax2.set_xticks(x, [f"V_G={v:.3f} V" for v in vgs])
ax2.set_yscale("log")
ax2.set_ylabel("relative discrepancy (%) [log scale]")
ax2.set_title("Independent network vs coupled small-signal")
ax2.legend(fontsize=8)
ax2.grid(True, axis="y", alpha=0.25)
fig.suptitle("Conditioned small-signal exploratory audit (archived C10/C11 states)")
fig.savefig(OUT / "conditioned_small_signal_audit.png", dpi=180)
fig.savefig(OUT / "conditioned_small_signal_audit.pdf")
