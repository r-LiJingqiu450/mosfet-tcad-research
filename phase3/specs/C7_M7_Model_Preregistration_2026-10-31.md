# C-7 — Model M7: Pre-registration of Predictions for Unseen States

**Date:** 2026-10-31. Committed **before any 2-D run of the test states.**

| Item | Location |
|---|---|
| Model | `phase3/c7/model_m7.py` (frozen at this commit) |
| Predictions | `results/phase3/C7/m7/m7_predictions.json` |
| Test jobs | `results/phase3/C7/jobs_m7_test.json` (deck W, abrupt, L = 2 µm, level 0) |

## Model M7 (summary; full statement in the module docstring)

M7 has no fitted parameters. It has four ingredients:
1. exact 1-D far field;
2. the exact depletion-approximation natural length l₀ of the two-layer strip;
3. a generalised quasi-2-D surface-potential equation, ε_Si W_eff ψ'' = Q_s(ψ) − C_ox(ψ_g − ψ), with W_eff fixed by l₀;
4. Ward's thin-inversion-layer corner-charge formula with his 90° depletion-field factor.

The result is in Ward's slice form. Current spreading is not modelled.

## Training (in-sample) comparison

This was seen before freezing. **It is not validation.**

| V_G (V) | 2-D ΔL_net (nm) | M7 (nm) | r = 2-D/M7 | Ward §6.4 (nm) |
|---|---|---|---|---|
| 1.0 | 369.9 | 275.0 | 1.345 | 142.7 |
| 1.25 | 153.7 | 136.5 | 1.126 | 98.1 |
| 1.5 | 92.9 | 94.4 | 0.984 | 70.9 |
| 2.0 | 53.4 | 59.4 | 0.899 | 46.6 |
| 3.0 | 29.9 | 32.9 | 0.910 | 27.9 |

**Known deficiency** (from the profile comparison, training data). M7's surface-potential excess is:
- too high within about 60 nm of the junction (2 nm: 224 vs 169 mV; 10 nm: 190 vs 108 mV, at 1.0 V);
- too low in the far tail (200 nm: 22 vs 40 mV; 400 nm: 3 vs 12 mV).

The 2-D profile is two-scale: an oxide-fringe-like layer of a few nm, plus a long depletion tail whose local decay length (124–170 nm) is longer than l₀ (113 nm). M7 has one intermediate scale. Near threshold the tail dominates, so M7 under-predicts.

M7 was **not** revised after this comparison.

## Test states (never simulated before)

| State | t_ox (nm) | N_A (cm⁻³) | V_G (V) | Inversion level (target N_s) | M7 ΔL (nm) | Ward ΔL (nm) |
|---|---|---|---|---|---|---|
| S0 | 20 | 1e16 | 1.1 | Unseen bias (N_s = 4.5e10) | 195.8 | 126.1 |
| S1 | 10 | 1e16 | 0.8585 | Threshold (6.2e9) | 276.0 | 138.3 |
| S1 | 10 | 1e16 | 1.0291 | Moderate (1.6e11) | 148.6 | 92.1 |
| S1 | 10 | 1e16 | 1.435 | Strong (9.0e11) | 65.8 | 40.6 |
| S2 | 40 | 1e16 | 1.2831 | Threshold | 278.0 | 156.3 |
| S2 | 40 | 1e16 | 1.692 | Moderate | 129.0 | 112.2 |
| S2 | 40 | 1e16 | 3.1303 | Strong | 55.6 | 59.2 |
| S3 | 20 | 4e16 | 1.3536 | Threshold | 155.2 | 81.4 |
| S3 | 20 | 4e16 | 1.6162 | Moderate | 81.2 | 65.7 |
| S3 | 20 | 4e16 | 2.3671 | Strong | 40.9 | 39.1 |

## Decision rules (fixed now)

r = ΔL_2D,net / ΔL_M7. The 2-D values are at mesh level 0. The C-6 mesh uncertainty is about 0.5 nm per end, which is negligible against these bands.

| Rule | Scope | Outcome |
|---|---|---|
| **R1** | The 6 moderate/strong states of S1–S3 | **SUPPORTED** if every \|r − 1\| ≤ 0.15. **REFUTED** if any \|r − 1\| > 0.25. Otherwise PARTIAL |
| **R2** | The 3 threshold states | M7's known under-prediction is expected to persist, with r ∈ [1.1, 1.6] |
| **R2** | S0 at 1.1 V | r ∈ [1.0, 1.4] |
| **R2** | All four | **REFUTED** if any r lies outside [0.9, 1.8] |
| **R3** | All 10 states, M7 vs Ward | M7 is **better than Ward** if Σ\|ln r_M7\| < Σ\|ln r_Ward\| |
| **R4** | Physics trends M7 predicts | See below |

**R4 trends.**
- At threshold, ΔL is nearly independent of t_ox: M7 gives ΔL(t_ox = 40)/ΔL(t_ox = 10) = 1.007. The trend is **confirmed** if the 2-D ratio lies in [0.8, 1.25].
- At threshold, raising N_A from 1e16 to 4e16 roughly halves ΔL: M7 gives 155.2/275.0 = 0.56. The 2-D ratio uses the C-6 reference of 369.9 nm for N_A = 1e16 and is **confirmed** if it lies in [0.4, 0.7].

Whatever the outcome, no parameter of M7 will be adjusted to these test data in this phase. A revised model would be a new model with new predictions.
