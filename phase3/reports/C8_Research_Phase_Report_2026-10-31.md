# C-8 Research Phase Report — Theory of the End Correction and Its Consequences for Extraction

**Date:** 2026-10-31. **Branch:** `phase3/C8-theory-extraction-2026-10-31`, from C-7 `2cfa8e0`.

**Evidence labels:**
- **[P]** proved or exact;
- **[D]** derived approximation;
- **[V]** verified numerics;
- **[PR]** pre-registered test, committed in `8650c76` before any C-8 run;
- **[O]** post-hoc observation on development data;
- **[H]** hypothesis;
- **[E1/E2/E3]** literature level.

**Physics:** F0 (Boltzmann, frozen holes, constant μ), zero bias, abrupt junctions, Ward truncation (deck W), unless stated.

**Compute:** 27 new states, about 10 CPU-minutes. All completed.

## 1. Answers to the phase questions

### 1.1 Strongest contribution now supported by evidence

The end correction is an exact, L-independent length ΔL(V_G) (C-6 P4–P6). It now has a **two-limit physical structure** with one empirical crossover. It also has **closed-form, verified consequences** for three standard characterization methods:
- Terada–Muta length/resistance extraction;
- linear-extrapolation threshold;
- split-C–V mobility.

Purely electrostatic, intrinsic artefacts of measurable size appear in all three for a device with N_A = 1e18 and t_ox = 2 nm. The consequence results are the most practically significant output of the project so far.

### 1.2 Does the scaling law have a credible theoretical foundation?

**Partly.**

**Weak-inversion limit (α∞ → ∞) [PR + D].**
- The argument (pre-registered): electrons are electrostatically passive except in a screened region near the junction, and that region vanishes as N∞ → 0. The end correction must therefore saturate at an **electron-free electrostatic limit**.
- Pre-registered test E8-1: **H-sat SUPPORTED, H-log REFUTED.** From α∞ = 15 to 150:

  | Structure | Δκ (H-log predicted +2.30) | ΔL (nm), α∞ = 15 / 50 / 150 |
  |---|---|---|
  | 1e16 / 20 nm | +0.40 | 410.7 / 418.3 / 416.3 |
  | 1e18 / 2 nm | +0.33 | 49.0 / 50.0 / 49.9 |

- ΔL in nanometres is constant to about 2 %. The small residual rise of κ comes from l₀ shrinking with W_d.
- In this limit ΔL/W_d = 1.43 → 1.61 (1e16) and 1.45 → 1.59 (1e18), even though the two devices differ 100× in doping and 10× in oxide.
- The C-7 development trend κ ≈ ln α∞ + 1.6 was a **crossover effect, not an asymptote**.

**Strong-inversion limit (α∞ → 0) [D + O].**
- The ratio of ΔL_2D to Yau's deep-junction charge-sharing estimate, ΔL_Yau = (W_d/2)·Q_B/(qN∞), depends only on α∞ across all structures:

  | α∞ | ΔL_2D/ΔL_Yau |
  |---|---|
  | 0.9 | 0.47–0.50 |
  | 0.6 | 0.63–0.64 |
  | 0.5 | 0.70–0.72 |
  | 0.27 | 0.87–0.88 |
  | 0.14 | 1.02 |

- So in strong inversion the conductance end correction tends to the classical **charge-sharing** result. The excess inversion charge at each end equals the depletion charge in half a depletion width.
- Equivalently, ΔL·Q_inv → Q_B·W_d/2.
- **Status:** post-hoc on development data. The approach is still visibly incomplete at α∞ = 0.14. Smaller α∞ could not be tested within F0 physics, because the channel surface potential reaches ψ_bi (the junction step closes) and Boltzmann statistics fail.

**The crossover** between the two limits (α∞ ≈ 0.5–5) is empirical. κ(α∞) is still a tabulated curve between two physically identified limits. It is **not a derived law**.

**Why κ collapses across structures [D].** Both limits scale with W_d:
- the subthreshold limit is ≈ 1.5–1.6 W_d;
- the charge-sharing limit is (W_d/2)·Q_B/Q_inv;
- l₀ ∝ W_d up to a weak t_ox correction (W_d/l₀ ≈ 2.7 for t_ox ≪ W_d).

The ratio of depletion to inversion field, α∞, is the natural interpolation variable. Normalising by l₀ rather than W_d is therefore not fundamental.

### 1.3 Has a genuinely predictive model emerged?

**Partially.**

**Predictive, verified:**
- the exact affine law;
- the two limits:
  - weak inversion: ΔL_sub ≈ 1.5–1.6 W_d, tested on two structures 100× apart in doping;
  - strong inversion: the Yau deep-junction asymptote, approached but not reached in the testable range;
- the closed-form extraction consequences of §2.

**Not predictive without data:** the crossover. The C-7 κ-interpolation passed 6/8 unseen states within ±7 %. It failed its registered criterion at α∞ = 0.2 (−13 % and −15 %). **Those two failures stand as recorded.**
- **[I]** Post hoc, those states have ψ_bi − ψ_s∞ = −1.0 V_T and +0.13 V_T: the junction potential step has closed. They are also where Boltzmann statistics are invalid (surface density ≥ 1e20 cm⁻³).
- In the Yau-normalised variable they lie close to the other structures at similar α∞ (0.88 and 0.92, against 0.87–0.88 at α∞ = 0.27).
- This explains them but does not rescue the registered test. A model of the strong-inversion branch needs (ψ_bi − ψ_s∞)/V_T as a second parameter.

### 1.4 Physical assumptions that limit applicability

- Boltzmann statistics with no quantum confinement. Invalid for scaled devices in strong inversion; relevant to the α∞ = 0.2 failures.
- Uniform channel doping (no halo or retrograde profile); abrupt vertical well boundary.
- Grading: σ ≤ 30 nm lowers near-threshold ΔL by ≤ 11 % (C-7).
- Constant μ. Lombardi post-processing changes near-threshold ΔL by +1 % (C-7). The extraction artefacts below are computed with constant μ precisely so that they isolate electrostatics.
- Zero drain bias; 2-D.
- The junction step ψ_bi − ψ_s∞ must stay positive, by several V_T.

### 1.5 Is the contribution distinguishable from existing literature?

| Item | Status |
|---|---|
| Quasi-2-D natural length for exponential threshold roll-off (Young 1989; Liu et al. 1993; BSIM3, E1 via summary) | Established. Different object |
| Yau charge sharing (1974, E3) | Established as a threshold model. **Our strong-inversion limit recovers it as a conductance end correction.** That is a recognition, not a novelty |
| Effective channel length as an electrical, model-dependent quantity (Taur 2000, E2; Narayanan et al. 1995, E2) | Established |
| Multi-length correction of split-C–V for ΔL and parasitics (Romanjek et al. 2004, E2 title only) | Established practice |
| Low-inversion split-C–V anomalies (Thomas et al., SSE 2009, E1 via tool summary) | Reported as **unexplained** ("further work is needed"); attributed partly to trapping. The intrinsic end-charge mechanism is not mentioned |
| ΔL as an exact L-independent length with a saturating subthreshold limit ≈ 1.5 W_d, crossing over to charge sharing | Not found in any source read (E1/E2). **Potentially original** |
| The three closed-form extraction artefacts (§2), which link spurious R_sd, algebraic V_T(L) roll-off and near-threshold split-C–V mobility error to the gate-dependent intrinsic ΔL | Not found in any source read (E1/E2). **Potentially original** |

All novelty statements are bounded by the unread paywalled texts: Taur 2000, Romanjek 2004, Narayanan 1995, Liu 1993, Yau 1974, Terada–Muta 1979.

### 1.6 Is further investment justified?

**Yes, for a narrower and more practical paper than C-7 envisaged.**

The scaling-law framing is weakened:
- its log branch was refuted;
- its crossover is empirical;
- its strong-inversion limit is classical.

The extraction-consequence framing is strengthened:
- exact formulas;
- quantitative, significant magnitudes for a modern-like device;
- a constructive remedy (§2.4);
- a candidate explanation for a reported unexplained anomaly (§2.3).

### 1.7 Recommended direction

Pursue **"Intrinsic source/drain electrostatics as a hidden bias in MOSFET characterization: exact zero-bias theory and its consequences for Terada–Muta, linear-threshold and split-C–V extraction"**, in three stages:
1. Literature gate on the six paywalled texts.
2. One realism stage (needs approval if it exceeds about 1 CPU-hour):
   - halo and retrograde profiles;
   - Fermi–Dirac statistics with a quantum correction;
   - mobility models for the split-C–V case.
   This decides whether the artefact magnitudes survive.
3. Only then a manuscript. Targets: *IEEE T-ED* if stages 1 and 2 hold, otherwise *Solid-State Electronics*.

## 2. Consequences for extraction methods (E8-2, scaled device N_A = 1e18, t_ox = 2 nm)

### 2.1 Terada–Muta [P][V]

The exact intersection of the 1/G–L lines at V₁ < V₂ is

R_sd,TM = 2[ΔL(V₁) − ΔL(V₂)]/[σ∞(V₂) − σ∞(V₁)] > 0,  ΔL_TM > ΔL(V₁).

The pre-registered signs are confirmed. The true series resistance is zero; the method reports:

| V₁ / V₂ (V) | True ΔL (nm) | ΔL_TM (nm) | R_sd,TM (Ω·µm) |
|---|---|---|---|
| 1.32 / 1.52 | 39.8 / 19.2 | 42.0 | 142 |
| 1.40 / 1.60 | 28.8 / 15.6 | 33.2 | 74 |
| 1.40 / 2.00 | 28.8 / 8.0 | 30.9 | 34 |
| 1.60 / 2.00 | 15.6 / 8.0 | 19.8 | 18 |

These are the same order as real S/D resistances in such technologies, about 10²  Ω·µm.

### 2.2 Linear-extrapolation threshold [D][V]

To first order in ΔL/L:

ΔV_T,lin(L) = (2σ∞²/(s²L))·dΔL/dV_G, with s = dσ∞/dV_G.

| Check | Value |
|---|---|
| Formula at L = 1 µm | −7.23 mV |
| Direct extraction at L = 1 µm | −7.43 mV |
| ΔV_T·L, constant from 0.5 to 10 µm | −7.2 to −7.6 mV·µm |
| L = 0.25 µm | −32.5 mV |

This is an **algebraic 1/L roll-off at lengths where the exponential quasi-2-D roll-off is negligible** (L/l₀ ≥ 19, e^{−L/2l₀} < 1e-4). The pre-registered prediction is confirmed. With constant μ, g_m increases monotonically, so the "maximum-g_m" point is the upper edge of the sweep (2.0 V). The formula covers any extrapolation point.

### 2.3 Split-C–V mobility [D][V]

The gate-modulated end charge length ΔL_C, from ∫C_gc dV_G with V_lo = 1.08 V (N∞ = 9.7e7 cm⁻²), greatly exceeds ΔL near threshold:

| V_G (V) | N∞ (cm⁻²) | ΔL_C (nm) | ΔL (nm) |
|---|---|---|---|
| 1.32 | 1.4e11 | 268 | 40 |
| 1.40 | 5.5e11 | 89 | 29 |
| 1.60 | 2.2e12 | 25 | 16 |

The locally inverted end regions dominate C_gc where the far channel is barely inverted. Mobility error relative to an infinitely long device measured the same way:

| N∞ (cm⁻²) | L = 1 µm | L = 2 µm | L = 10 µm |
|---|---|---|---|
| 1.4e11 | −29 % | −18 % | −4.3 % |
| 5.5e11 | −10 % | −5.5 % | −1.2 % |
| 2.2e12 | −1.8 % | −0.9 % | −0.2 % |

The pre-registered sign is confirmed (μ underestimated).

**[H]** This is a candidate intrinsic contributor to the low-inversion mobility drop that Thomas et al. (2009, L = 10 µm) report as partly unexplained. It is untested against their data.

### 2.4 Constructive consequence [P]

By the affine law, d(1/G)/dL = 1/σ∞ and d(Q_gc)/dL = qN∞ exactly; the end terms cancel. A **two-length differential** analysis therefore removes all three artefacts exactly at zero bias, in the spirit of existing multi-length split-C–V practice (Romanjek 2004, E2). This gives that practice an exact justification and predicts when single-length methods fail.

## 3. Experiment record (C-8)

| Experiment | States | Status | Result |
|---|---|---|---|
| E8-1 deep-subthreshold κ [PR] | 6 | All OK | **H-sat SUPPORTED, H-log REFUTED** |
| E8-2 scaled-device V_G sweep | 21 | All OK | Analysis complete; pre-registered signs (a)–(c) **confirmed** |
| C-7 κ-law registered test (carried over) | — | — | **PARTIAL; failures at α∞ = 0.2 stand** (explained post hoc, not rescued) |
| Strong-inversion Yau asymptote | 0 new | — | Post-hoc on development data; incomplete approach; **not independently tested** |

No failures, timeouts or cancellations in C-8.

## 4. Artifacts

| Kind | Location |
|---|---|
| Code | `phase3/c8/e8_subth_eval.py`, `phase3/c8/e8_extraction.py`. Simulations use `phase3/c7/endprobe.py`, unchanged |
| Registration | `phase3/specs/C8_Preregistration_2026-10-31.md` |
| Data | `results/phase3/C8/subth/` (6 state JSONs, `preregistration_data.json`, `e8_1_result.json`, logs, RUNS files); `results/phase3/C8/sweep/` (21 states, logs); `results/phase3/C8/extraction/` (`e8_2_result.json`, `e8_2_figure.png`, stdout) |
| Archives and verification | See the delivery note |
