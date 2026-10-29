#!/usr/bin/env python3
"""Make Phase 2C audit plots from frozen JSON evidence."""

from pathlib import Path
import json
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "verification/plots"
OUT.mkdir(parents=True, exist_ok=True)


def main():
    s2 = json.loads((ROOT / "verification/phase2c_s2_conservation.json").read_text())
    records = s2["records"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for grid in ("coarse", "medium", "production"):
        rows = sorted((r for r in records if r["grid"] == grid), key=lambda r: r["V_a"])
        axes[0].plot([r["V_a"] for r in rows], [r["direct"]["Jn_over_qU"] for r in rows], "o-", label=f"{grid} direct")
        axes[0].plot([r["V_a"] for r in rows], [r["high_precision_qfp"]["Jn_over_qU"] for r in rows], "--", label=f"{grid} high precision")
        axes[1].plot([r["V_a"] for r in rows], [r["direct"]["total_current_spread_relative"] for r in rows], "o-", label=grid)
    axes[0].axhline(1.0, color="k", ls="--")
    axes[0].set_xlabel("V_a (V)"); axes[0].set_ylabel("ΔJ_n / (q∫U dx)"); axes[0].grid(True, alpha=.3)
    axes[0].legend(fontsize=7, ncol=2)
    axes[1].set_yscale("log"); axes[1].set_xlabel("V_a (V)"); axes[1].set_ylabel("total-current spread / scale"); axes[1].grid(True, alpha=.3)
    axes[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "phase2c_s2_conservation.png", dpi=180); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for grid in ("coarse", "medium", "production"):
        rows = sorted((r for r in records if r["grid"] == grid), key=lambda r: r["V_a"])
        ax.plot([r["V_a"] for r in rows], [r["max_direct_vs_high_returned_Jn"] for r in rows], "o-", label=f"{grid} Jn")
        ax.plot([r["V_a"] for r in rows], [r["max_direct_vs_high_returned_Jp"] for r in rows], "--", label=f"{grid} Jp")
    ax.set_yscale("log"); ax.set_xlabel("V_a (V)"); ax.set_ylabel("direct SG vs high precision returned-state flux (A/m²)"); ax.grid(True, alpha=.3); ax.legend(fontsize=8, ncol=2)
    fig.tight_layout(); fig.savefig(OUT / "phase2c_s2_formula_error.png", dpi=180); plt.close(fig)
    print(str(OUT))


if __name__ == "__main__":
    main()
