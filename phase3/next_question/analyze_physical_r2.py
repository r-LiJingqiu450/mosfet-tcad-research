"""Post-process the bounded physical R-II finite-bias probe.

The script has no solver dependency.  It tests the conditional question
whether R(L)=1/(I/V_D) is affine in the metallurgical channel length and
whether its V_D -> 0 slope agrees with the independently produced zero-bias
network conductance.  The quadratic extrapolation is a diagnostic, not a
claimed physical law.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

Q_E = 1.602176634e-19
L_UM = np.array([1.0, 2.0, 4.0])


def load_json(path: Path):
    d = json.loads(path.read_text())
    if not d["finite_states"]:
        raise ValueError(f"no finite states in {path}")
    return d


def zero_bias_reference(root: Path, vg: float):
    vals = []
    for L in (1, 2):
        p = root / "results/phase3/independent_discovery/graded_probe" / (
            f"e7_R_L{L}_lev0_VG{vg:g}_tox2_NA1e+18_hw0.5_ov0_sg15_dpsi0.json")
        d = json.loads(p.read_text())
        c = d["result"]["column_sheet_density"]
        x = np.asarray(c["x_cm"], float)
        n = np.asarray(c["N_cm2"], float)
        lm = L * 1e-4
        m = (x >= -1e-14) & (x <= lm + 1e-14)
        q = Q_E * np.trapezoid(n[m], x[m])
        vals.append((float(d["result"]["const_mu"]["G_net"]), float(q)))
    dR = (1.0 / vals[1][0] - 1.0 / vals[0][0]) / 1e-4
    dQ = (vals[1][1] - vals[0][1]) / 1e-4
    return dict(G_L1=vals[0][0], G_L2=vals[1][0], dR_dL=dR, dQ_dL=dQ,
                mu=1.0 / (dR * dQ), Q_L1=vals[0][1], Q_L2=vals[1][1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    ap.add_argument("--input-dir", type=Path, default=None)
    ap.add_argument("--out-dir", type=Path, required=True)
    a = ap.parse_args()
    root = a.root.resolve()
    inp = (a.input_dir or root / "results/phase3/next_question/physical_r2_finite_bias").resolve()
    out = a.out_dir.resolve(); out.mkdir(parents=True, exist_ok=True)
    cases = {1.32: "VD0p1_5_20mV", 1.60: "VD1_5_20mV"}
    # Add the 0.1-mV point for VG=1.60 from its separate run.
    rows = []
    summaries = {}
    for vg, suffix in cases.items():
        ds = [load_json(inp / f"r2_finite_L{int(L)}_lev0_VG{vg:g}_{suffix}.json") for L in L_UM]
        by_vd = {}
        for d in ds:
            for s in d["finite_states"]:
                by_vd.setdefault(float(s["vd_V"]), []).append(s)
        if vg == 1.60:
            low = [load_json(inp / f"r2_finite_L{int(L)}_lev0_VG1.6_VD0p1mV.json") for L in L_UM]
            for d in low:
                s = d["finite_states"][0]
                by_vd.setdefault(float(s["vd_V"]), []).append(s)
        ref = zero_bias_reference(root, vg)
        srows = []
        for vd in sorted(by_vd):
            states = by_vd[vd]
            if len(states) != 3:
                raise ValueError(f"VG={vg}, VD={vd}: expected 3 lengths, got {len(states)}")
            G = np.array([float(s["G_eff_S_per_cm"]) for s in states])
            R = 1.0 / G
            A = np.c_[np.ones(3), L_UM * 1e-4]
            coef, *_ = np.linalg.lstsq(A, R, rcond=None)
            resid = R - A @ coef
            pair12 = (R[1] - R[0]) / 1e-4
            pair24 = (R[2] - R[1]) / 2e-4
            max_kcl = max(float(s["kcl_rel"]) for s in states)
            max_xsec = max(float(s["cross_section_abs_spread_rel"]) for s in states)
            row = dict(VG_V=vg, VD_V=vd, R_L1=R[0], R_L2=R[1], R_L4=R[2],
                       slope_affine_Sinv_per_cm=coef[1], intercept_Sinv=coef[0],
                       pair12_slope=pair12, pair24_slope=pair24,
                       slope_pair_rel_diff=abs(pair24 - pair12) / abs(coef[1]),
                       affine_max_rel_residual=max(abs(resid) / abs(R)),
                       max_KCL_rel=max_kcl, max_cross_section_spread_rel=max_xsec,
                       zero_bias_slope=ref["dR_dL"],
                       slope_rel_vs_zero_bias=coef[1] / ref["dR_dL"] - 1.0)
            rows.append(row); srows.append(row)
        # A low-order extrapolation is only used to check the limit; report its
        # order and data count so it cannot be mistaken for a theorem.
        x = np.array([r["VD_V"] for r in srows]); y = np.array([r["slope_affine_Sinv_per_cm"] for r in srows])
        deg = min(2, len(x) - 1)
        coeff = np.polyfit(x, y, deg)
        summaries[str(vg)] = dict(zero_bias_reference=ref, fit_degree=deg,
                                  fit_coefficients=coeff.tolist(),
                                  extrapolated_slope=float(np.polyval(coeff, 0.0)),
                                  extrapolated_rel_error=float(np.polyval(coeff, 0.0) / ref["dR_dL"] - 1.0),
                                  rows=srows)
    fields = list(rows[0])
    with (out / "physical_r2_finite_bias_rows.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)
    summary = {"status": "EXPLORATORY_CONDITIONAL_EVIDENCE", "configuration": {
        "deck": "R-II resolved n+ wells", "N_A_cm-3": 1e18, "sigma_nm": 15.0,
        "tox_nm": 2.0, "level": 0, "VG_V": [1.32, 1.60], "L_um": [1, 2, 4],
        "VD_V": [1e-4, 1e-3, 5e-3, 2e-2], "mobility": "constant 1350 cm2/(V s)",
        "potential_reservoir": "equilibrium PsiLat, s_couple=0; drain quasi-Fermi bias"},
              "equation": "R(L,VD)=1/(I(L,VD)/VD)=a(VD)+b(VD)L; b(0) compared with zero-bias network",
              "cases": summaries}
    (out / "physical_r2_finite_bias_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    # Compact figure, if matplotlib is available.
    try:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(9, 3.5))
        for vg, ss in summaries.items():
            rr = ss["rows"]; x = [r["VD_V"] * 1e3 for r in rr]; y = [r["slope_affine_Sinv_per_cm"] for r in rr]
            ax[0].plot(x, y, "o-", label=f"V_G={vg} V")
            ax[1].semilogy(x, [r["affine_max_rel_residual"] for r in rr], "o-", label=f"V_G={vg} V")
        ax[0].set(xlabel="V_D (mV)", ylabel="affine slope b (cm/S per cm)")
        ax[1].set(xlabel="V_D (mV)", ylabel="max relative residual")
        for aa in ax: aa.grid(True, alpha=.3); aa.legend()
        fig.tight_layout(); fig.savefig(out / "physical_r2_finite_bias.png", dpi=180); plt.close(fig)
    except Exception as e:
        (out / "plot_warning.txt").write_text(f"plot unavailable: {e}\n")
    print(json.dumps({"rows": len(rows), "out": str(out),
                      "extrapolated_rel_errors": {k: v["extrapolated_rel_error"] for k, v in summaries.items()}}, indent=2))


if __name__ == "__main__":
    main()
