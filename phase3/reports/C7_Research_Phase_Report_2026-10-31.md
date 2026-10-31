# C-7 Research Phase Report — Realism of the End Correction, a Failed Model, and a Scaling Law

**Date:** 2026-10-31. **Branch:** `phase3/C7-realism-model-2026-10-31`, from C-6 `08d61cb`.

**Evidence labels:**
- **[P]** proved;
- **[V]** verified numerical result;
- **[PR]** pre-registered test, with the commit made before the runs;
- **[O]** post-hoc observation;
- **[H]** hypothesis;
- **[I]** interpretation.

All 2-D results are zero-bias, F0 physics (Boltzmann, frozen holes, no generation/recombination), with constant μ unless marked "Lombardi". **About 55 CPU-minutes were used, within the approved ~2 CPU-hours.**

## 1. Strongest new finding — an end-correction scaling law [PR][O]

The intrinsic source/drain end correction ΔL (per end, L-independent; C-6 P4–P6) collapses onto a single curve:

  **ΔL = l₀ · κ(α∞)**

where:
- **l₀** is the exact decay length of the lowest end-layer mode in the depletion approximation. It is the lowest root of tan(kW_d) = −(ε_Si/ε_ox) tan(k t_ox), with l₀ = 1/k and W_d the 1-D depletion width. It has no fitted parameter.
- **α∞** is Ward's inversion parameter of the far-field channel ([W88] p. 139): the ratio of depletion field to inversion field, from the exact 1-D solution. α∞ ≈ 5 is threshold; α∞ < 0.3 is strong inversion.
- **κ** is a dimensionless function. Its values (from `kappa_table_with_B.csv`) are:

  | α∞ | κ |
  |---|---|
  | 0.2 | ≈ 0.3 |
  | 0.6 | ≈ 0.87 |
  | 0.9 | ≈ 1.3 |
  | 3 | ≈ 2.8 |
  | 5.2 | ≈ 3.3 |
  | 7.6 | ≈ 3.7 |

**Evidence (all deck W: abrupt junction, Ward truncation, which matches resolved wells to ≤ 2 %):**

1. **Development set.** 17 states, N_A = 1e16–4e16, t_ox = 10–40 nm. At fixed α∞ the spread is:
   - α∞ = 5.23: κ = 3.03–3.30 (three oxides);
   - α∞ = 0.90: κ = 1.18–1.34.
2. **Addendum A [PR], SUPPORTED.** Scaled device N_A = 1e18, t_ox = 2 nm (100× doping, 10× thinner oxide):

   | State | Pre-registered band (nm) | 2-D result (nm) |
   |---|---|---|
   | Threshold | 38.6–47.6 | **43.8** |
   | Moderate | 14.7–19.4 | **17.2** |

3. **Addendum B [PR], PARTIAL under the registered rule.** Two never-simulated devices (N_A = 3e17 with t_ox = 5 nm; N_A = 1e17 with t_ox = 10 nm), four α∞ values each, predicted from the interpolated development curve:
   - for α∞ ∈ {3.0, 1.2, 0.6}, all six errors lie in **+6.4 % to −7.0 %**;
   - at α∞ = 0.2 the errors are **−15 % and −13 %**: above the 12 % "supported" limit, below the 20 % "refuted" limit.

   [I] At those two points the Boltzmann surface density is 2.8e20 and 8.6e19 cm⁻³, and ψ_s∞ is at or above ψ_bi. That is degenerate-inversion territory, where F0 physics is no longer valid.

   The two new devices give nearly identical κ at every α∞: 2.82/2.85, 1.61/1.59, 0.87/0.86, 0.29/0.30.

**Meaning.** In weak-to-moderate inversion the end correction is a fixed multiple of an electrostatic eigen-length.
- At threshold ΔL ≈ 3.3 l₀ ≈ 1.1–1.5 W_d per end, about 44 nm even for N_A = 1e18 and t_ox = 2 nm.
- κ depends on bias only through α∞.
- **[H]** In weak inversion κ grows about logarithmically in α∞ (slope dκ/d ln α∞ ≈ 1.0 between α∞ = 5.2 and 7.6), as expected from integrating an exponential surface-potential tail against Boltzmann charge. There is no derivation yet. This is the open theory problem (§9).

## 2. Hypotheses: supported, modified, rejected

| Hypothesis | Test | Outcome |
|---|---|---|
| Near-threshold ΔL survives realistic junctions (registered S1/S0) | 24 states, deck R | **SUPPORTED (S1).** §3 |
| M7 (parameter-free: quasi-2-D surface potential + Ward corner charge) predicts ΔL within ±15 % in moderate/strong inversion on unseen structures (R1) | 10 unseen states [PR] | **REFUTED.** Max deviation 44 % (N_A = 4e16, moderate) |
| M7 under-predicts at threshold by r ∈ [1.1, 1.6] (R2) | Same | **As predicted:** r = 1.15–1.55; 1.30 at the unseen bias |
| M7 is better than Ward §6.4 (R3) | Same | **SUPPORTED.** Σ\|ln r\| = 2.19 vs 5.52 |
| M7 trends (R4): ΔL_th independent of t_ox; halved by 4× N_A | Same | **t_ox trend REFUTED** (2-D ratio t_ox 40/10 = 1.32). **N_A trend confirmed** (0.65) |
| Ward's corner-charge relation is accurate given the true surface potential (C-6 finding 6, formed on C-6 data) | Applied unchanged to 9 unseen-structure states [O] | **SUPPORTED near threshold:** −0 % to −3 % of the 2-D slice value for α∞ ≥ 1.38; −6 to −8 % at α∞ = 0.9. **Fails in strong inversion:** −17 to −28 % |
| Scaling law ΔL = l₀ κ(α∞) | Addendum A [PR]; Addendum B [PR] | **SUPPORTED** (A); **PARTIAL** (B: ±7 % for α∞ ≥ 0.6; −13/−15 % at α∞ = 0.2) |
| C-6 statement: Ward's strong-inversion agreement is "partly a cancellation" | Surface-layer slice diagnostic | **Modified.** At 2.0 V the surface-layer (y ≤ 20 nm) slice gives 40.8 nm, the full slice 62.2 nm and the network 52.4 nm. Charge below 20 nm is *partly* effective. The slice/spreading split is an exact identity, not a unique causal decomposition. I withdraw "fortuitous" and keep only "the two terms have opposite signs" |
| P7: the end term is equivalent to a constant threshold shift | C-6 | Remains **rejected** |

**Failure diagnosis of M7 [V].** M7's surface-potential profile is single-scale:
- too high within about 60 nm of the junction;
- too short in the depletion tail (e.g., 1.0 V, 400 nm: 3 mV vs 12 mV in 2-D).

Replacing only the profile by the true 2-D ψ_s(x) brings ΔL to within 3 % near threshold for every structure (`m7_failure_diag.json`). The model's weak point is therefore **the prediction of ψ_s(x)**, specifically the nonlinear depletion tail. It is not Ward's charge relation. This is the regime Ward himself declared not constructible by his method ([W88] p. 138). M7 is not revised in this phase. A successor needs new predictions on new states.

## 3. Does the near-threshold end correction survive more realistic assumptions? Yes

All values are ΔL per end (nm); 1e16/20 nm device, deck R (resolved n⁺ wells), L = 2 µm.

| Case | V_G = 1.0 V (threshold) | 1.25 V | 2.0 V |
|---|---|---|---|
| Abrupt, no overlap | 370.2 | 153.4 | 52.4 |
| Overlap 10 nm / 30 nm | 370.2 / 370.2 | 153.5 / 153.5 | 52.6 / 52.7 |
| Lateral erfc grading σ = 10 / 30 nm | 353.6 / 337.6 | 136.1 / 116.8 | 31.9 / **−2.2** |
| Overlap 30 + grading 30 | **337.8** | 120.0 | 13.2 |
| Same, with Lombardi μ | **341.4** | 130.0 | 28.7 |

| Check, combined case at V_G = 1.0 V | ΔL (nm) |
|---|---|
| Level 1 | 337.73 |
| L = 4 µm | 337.76 |

**ψ_bi ± 50 mV** (a stand-in for degeneracy or band-gap narrowing in the well): ≤ 1 % change.

**Scaled device** (N_A = 1e18, t_ox = 2 nm):

| Case | ΔL at threshold (nm) |
|---|---|
| Ward truncation (deck W) | 43.8 |
| Resolved wells (deck R) | 43.0 |
| Resolved wells, σ = 5 nm | 38.8 |
| Level 1 | 43.77 |

**Registered verdict: S1 (survives).** For the combined case ΔL = 337.8 nm (constant μ) and 341.4 nm (Lombardi), both ≥ 76 nm. ΔL(1 V)/ΔL(2 V) = 25.6 (constant μ) and 11.9 (Lombardi).

**[I]**
- Gate overlap does not matter: the n⁺ surface is pinned whether or not the gate covers it.
- Grading and mobility change the **strong-inversion** value strongly. Grading even makes it negative, the ordinary "lateral-diffusion" series-resistance sign. They change the **near-threshold** value by only −9 % to +1 %.
- So the near-threshold end correction is an intrinsic electrostatic property. The strong-inversion value is junction-geometry-dominated, which is what conventional ΔL extraction measures.

## 4. Predictive model — status

| Model | Status |
|---|---|
| **M7** (physics-derived, parameter-free) | Failed its main pre-registered criterion. Diagnosis in §2 |
| **The κ-law** (ΔL = l₀ κ(α∞)) | Predictive within ±7 % on two unseen devices for α∞ ≥ 0.6, and inside the pre-registered bands on a 100×-scaled device |

The κ-law has one empirical universal curve, and its inputs are exact and closed-form. **It is semi-empirical:** κ is tabulated, not derived.

Its domain:
- abrupt junctions, Ward truncation;
- F0 physics;
- α∞ ≳ 0.5;
- grading of σ ≤ 30 nm lowers ΔL by ≤ 11 % at threshold.

## 5. What the literature already establishes

**Full text (E1):**
- **Ward 1988; Ward–Odeh–Cohen 1990** (C-5): asymptotic end theory, unevaluated; the inner regions are left open.
- **Huang, Liu, Jeng, Ko, Hu 1993, UCB/ERL M93/57** (BSIM3), read through a tool summary of an OCR text:
  - threshold roll-off θ(L) ∝ exp(−L/2l_t) + 2exp(−L/l_t), with a characteristic length l_t built from t_ox and X_dep;
  - η and D_VT0 fitted to data;
  - no 1/L term;
  - no discussion of effective channel length near threshold.
- **Cunha, Schneider, Galup-Montoro et al., NSTI Nanotech 2005**, read through a tool summary:
  - L and W are extracted from g_m/I_D in weak-to-moderate inversion at small V_DS;
  - ΔL ≈ 32 nm (NMOS, 0.35 µm), treated as bias-independent.

**Abstract only (E2):**
- **Taur, IEEE T-ED 2000:** 2-D-model interpretation of extracted L as current "injection points".
- **Narayanan et al., SSE 1995:** a 2-D numerical study; the junction-based definition of L is "questionable".
- **Ward, SIAM J. Appl. Math. 1992:** small-drain-bias I–V of buried-channel FETs with 2-D source/drain regions.

**Metadata only (E3):** Reyna–Ward 1991; Ward–Reyna–Odeh 1994; Liu et al., T-ED 1993; Young 1989; Yau 1974; Terada–Muta 1979; Chern et al. 1980; Ng–Brews 1990.

**What is therefore established:**
- the quasi-2-D natural length (Young; Liu/BSIM) **for the exponential threshold roll-off**;
- the idea that extracted L is an electrical, model-dependent quantity (Taur; Narayanan);
- Ward's asymptotic end framework.

**Not found in any source read** (absence in E2/E3 sources is not evidence):
- an L-independent linear-regime end correction scaling with the depletion eigen-length;
- its dependence on the inversion level α∞;
- any test of Ward's §6.4.

## 6. What remains potentially original

| Contribution | Status |
|---|---|
| The κ(α∞) scaling law for the intrinsic near-threshold end correction, with the exact two-layer eigen-length as its scale | `NOVELTY_POTENTIAL`, pending full-text reading of the E2/E3 items. The highest risk is in the quasi-2-D and channel-length-extraction literature |
| First quantitative test of Ward §6.4. Its corner-charge relation validated out of sample near threshold; its phenomenological profile identified as the failure | `NOVELTY_POTENTIAL` (E1 basis) |
| Exact zero-bias framework (C-6 P4–P6) as a verification method for reduced models | Methodological; not novel as mathematics |
| **[H]** The intrinsic electrostatic ΔL in weak/moderate inversion (≈ 44 nm per end at N_A = 1e18, t_ox = 2 nm) is the same order as the ΔL values extracted in weak/moderate inversion in practice (e.g., 32 nm in Cunha et al.) | Untested hypothesis with direct significance for extraction practice |

## 7. Most important uncertainties

1. **Physics scope.**
   - F0 assumptions: Boltzmann statistics, invalid for scaled devices in strong inversion and at the α∞ = 0.2 test points; no quantum confinement, which matters at t_ox = 2 nm and N_A = 1e18; frozen holes.
   - Lombardi mobility only in post-processing, without distance damping.
   - Uniform channel doping: no halo, retrograde or vertically graded junctions. The vertical junction is a box.
2. **κ has no derivation.** The law is semi-empirical; its logarithmic large-α∞ behaviour is only conjectured.
3. **Literature.** Paywalled full texts (Taur 2000, Narayanan 1995, Liu 1993, Ward 1992, Yau 1974, extraction-method papers) were not obtainable here.
4. **Numerics.**
   - Mesh changes are ≤ 0.2 nm (level 0 → 1 → 2, C-6; level 1 checks here). The C-6 level-2 sequence was not asymptotic, so the working uncertainty is about 0.5 nm per end.
   - Network G was cross-checked against DEVSIM's full Jacobian at selected states, always within its adjoint estimate.
   - A **defect found in C-6 output**: the `column_sheet_density` field written by `phase3/c6/wprobe.py` used misaligned depth weights. No C-6 conclusion used it. The C-7 code (`endprobe.py`) maps the weights correctly, with an assertion.
5. **Definitions.** ΔL is measured from the metallurgical junction. With graded junctions part of the change is definitional.
6. **2-D and zero bias only.**

## 8. Is there a credible path to an SCI publication?

**Yes, for a focused paper.** Working title: "Intrinsic source/drain end correction of the long-channel MOSFET conductance near threshold: exact zero-bias formulation, test of asymptotic theory, and a scaling law."

The evidence base is consistent:
- verified numerics;
- pre-registered tests, including failures;
- robustness to overlap, grading and mobility;
- out-of-sample validation across 100× in doping and 20× in oxide thickness.

**Venue:**
- *Solid-State Electronics* or *IEEE J-EDS* would be credible now, after the literature gate.
- **IEEE T-ED** becomes credible if two further pieces succeed:
  - (a) the consequences for standard channel-length/threshold extraction are demonstrated quantitatively;
  - (b) the law holds with realistic doping (halo/retrograde) and quantum/degeneracy corrections, or is bounded by them.

## 9. Next investment (recommended order)

1. **Theory of κ(α∞), analysis only.** Formulate the weak-inversion end problem with the nonlinear depletion boundary. Derive κ asymptotically for α∞ → ∞ (the logarithmic law) and for α∞ → 0. A successful derivation would convert the semi-empirical law into a predictive theory, and is the core of the paper.
2. **Extraction-impact analysis, with no new simulation needed.** By the exact affine law (C-6 P6), 1/G(L, V_G) = (L − 2ΔL(V_G))/σ∞(V_G). Feeding this into Terada–Muta, shift-and-ratio and g_m/I_D L-extraction analytically shows how much the intrinsic ΔL(V_G) biases each extracted L_eff, R_SD and V_T(L).
3. **Literature gate (LG-1).** I need the paywalled PDFs listed in §5.
4. **One small realism campaign**, needing approval if it exceeds about 1 CPU-hour:
   - halo/pocket and retrograde profiles;
   - degeneracy (Fermi–Dirac) and a quantum-confinement correction (density gradient) for the scaled device.

   **Abandon the T-ED framing if** these move the near-threshold κ by more than 30 %.

**Not recommended:** large sweeps, finite-bias campaigns, short-channel devices.

## 10. Artifacts

| Kind | Location |
|---|---|
| Code | `phase3/c7/endprobe.py`, `model_m7.py`, `m7_predictions.py`, `m7_failure_diag.py`, `c7_analysis.py`, `c7_scaling.py`, `kappa_predict.py` |
| Registrations | `phase3/specs/C7_Realism_Test_Registration_2026-10-31.md` (commit 29cca27, before any realism run); `phase3/specs/C7_M7_Model_Preregistration_2026-10-31.md` (M7: 264e968; Addendum A: d6d5e39; Addendum B: 6c0cef2, each before its runs) |
| Data | `results/phase3/C7/{dev, realism, m7, scaled, kappa_test}/`: every state JSON, logs, RUNS/RUNNER_END files and runner outputs |
| Analysis | `results/phase3/C7/analysis/`: c7_summary.json, realism_table.csv, m7_test_table.csv, m7_failure_diag.json, kappa_table.csv (development table as pre-registered; extra columns added later, κ and α∞ identical), kappa_table_with_B.csv, C7_figure.png, C7_kappa.png |
| Archives | See the delivery note (`C7_PACKAGE_VERIFICATION.md`) |

**Failures:** none in C-7. All 49 C-7 states completed. One launch error happened before any simulation: the first R smoke run's shell redirect pointed at a not-yet-created directory, and it was re-launched. The C-6 failures and cancellations are preserved in `results/phase3/C6/w1/CANCELLED.txt`.
