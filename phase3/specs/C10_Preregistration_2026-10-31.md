# C-10 Pre-registration — subthreshold end correction: similarity, bias range, oxide scaling, Fjeldly–Shur comparison

**Date:** 2026-10-31. **Branch:** `phase3/C10-subthreshold-theory-2026-10-31`, from `c41b5bf`.

**Status:** written and committed **before** any C-10 2-D run.

**Design numbers:** `phase3/c10/c10_design.py` writes `results/phase3/C10/design.json` and the job lists `jobs_{sim,sweep,tox,mesh}.json`, 22 states in all.

**Physics:** F0, zero bias, deck W (Ward truncation, abrupt junction), constant μ. The quantity is the conductance end correction ΔL per end (`dL_net_nm`).

**Definitions:**
- ψ̂ = (ψ_s∞ − ψ_b0)/V_T;
- L_D = (ε_Si V_T/(qN_A))^{1/2};
- t̂ = (ε_Si/ε_ox) t_ox/L_D;
- ĥ = h_w/L_D;
- ψ̂_bi = (ψ_bi − ψ_b0)/V_T.

## Prior-art input read before registration

Fjeldly & Shur, IEEE T-ED 40(1), 137–145 (1993). Full text supplied by the user; §§I–IV read.
- **Their effective length.** The subthreshold current uses the exact integral (2), ∫ exp(−ψ/V_th) dx. The integral is approximated by the flat-region length, L → L − x_s − x_d (eq. 8).
- **Their end lengths.** x_s and x_d are conventional 1-D abrupt-junction depletion widths: x_s ≈ [2ε_s(−ψ_s⁰)/(qN_a)]^{1/2} (eq. 9), where −ψ_s⁰ = ψ_bi − ψ_s∞ at V_DS = 0.
- **What is not validated.** The paper does not compare x_s with 2-D simulation.

Existing C-8 states (deck W, α∞ = 15–150) give ΔL_2D/x_s = 1.81–2.12 (device 1) and 2.61–3.36 (device 2). This was computed before registration and is disclosed here.

## Hypotheses and decision rules

**P10-1 — Similarity theorem [P] (C-10 note §2).**
- Device 2' is (1e18, t_ox = 2 nm, h_w = 47.27 nm, ψ_bi shifted by V_T ln(1e16/1e18) = −0.11905 V). It is the electron-free similarity image of device 1 = (1e16, 20 nm, h_w = 469.05 nm). The design matches t̂ = 1.4675, ψ̂_bi = 36.84 and ĥ = 11.56 vs 11.47.
- Both are run at V_G = 0.7835 V and 0.6214 V, which give the same ψ̂ = 21.02 and 16.0.
- **SUPPORTED** if |(ΔL_2'/L_D2)/(ΔL_1/L_D1) − 1| ≤ 2 % at both biases. Otherwise **REJECTED**, and the residual is reported.

**P10-2 — Factor decomposition (descriptive; no decision rule).** At V_G = 0.7835 V, device 2 is run as is (A), with ψ̂_bi matched (B), with ĥ matched (C), and with both matched (D = device 2'). ΔL/L_D is reported against device 1 for each.

**P10-3 — Oxide scaling at fixed W_d.** Device 2 with t_ox = 6 nm and 1 nm is run at the ψ̂ of the existing α = 150 state (30.50) and at 21.02. The ratio ΔL(t_ox)/ΔL(2 nm) is formed at the same ψ̂. The baselines are the existing state (1.08 V) and run A.

| Hypothesis | Predicted ratio |
|---|---|
| **H-Wd** (C-8/C-9: ΔL ≈ c·W_d, independent of t_ox) | [0.95, 1.05] |
| **H-l0** (C-10: ΔL ∝ l₀ at fixed ψ̂) | l₀-ratio × (1 ± 0.10) |

The l₀-ratios are 1.2069 and 1.2201 for 6 nm, and 0.9281 and 0.9186 for 1 nm.

**Verdict rule:** decided on the two 6-nm states, where the intervals are disjoint. The 1-nm states are reported, because their intervals overlap.

**P10-4 — Is the subthreshold "saturation" a true limit?** Bias sweeps:
- device 2 (L = 1 µm): ψ̂ = 12, 15, 18, 21.02, 24, 27, plus the existing 30.50, 32.70, 35.11;
- device 1 (L = 2 µm): ψ̂ = 12, 15, 18, plus the existing 21.02, 23.22, 25.62.

| Outcome | Condition, for each device |
|---|---|
| **H-sat** (true electron-free saturation) | every ΔL within ±5 % of ΔL(α = 150) over ψ̂ ∈ [12, ψ̂(α = 150)] |
| **H-max** (broad maximum; ΔL falls in deep subthreshold) | ΔL(ψ̂ = 12) ≤ 0.92 · max ΔL over the sweep |

Anything else is **INCONCLUSIVE**.

**P10-5 — Asymptotic identity [C].** On every new state with α ≥ 100:

|λ_fit·[Ein(B_fit/V_T) + r] / ΔL_net − 1| ≤ 5 %

with λ_fit and B_fit from the surface-potential tail fit (`phase3/c10/tail_analysis.py`), and r = C_ox V_T/Q_dep,∞.

**P10-6 — Fjeldly–Shur eq. (9) as the zero-bias end correction.**
- **H-FS** (|ΔL_2D/x_s − 1| ≤ 20 %) is tested on every new state.
- **Stated expectation:** it fails, with ratio > 1.5 everywhere.

**Mesh check (no decision rule):** level 1 vs level 0 for device 1 at α = 150 and for device 2, run A.

**Failures:** any non-converged or failed run is kept and reported. No state is dropped or re-defined after its result is seen.
