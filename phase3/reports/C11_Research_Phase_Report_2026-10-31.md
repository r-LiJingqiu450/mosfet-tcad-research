# C-11 Research Phase Report — Is the Subthreshold End Correction Predictable? An Independent Check, and a Decision on the Research Line

**Date:** 2026-10-31. **Branch:** `phase3/C11-connection-coefficient-2026-10-31`, from C-10 `da913b5`.

**Pre-registration:** `phase3/specs/C11_Preregistration_2026-10-31.md`, commit `d34f58a`, before evaluation.

**Labels:**
- [P] proved; [C] conditional identity; [V] numerical verification within the F0 model;
- [PR] pre-registered; [PH] post hoc; [O] observation; [I] interpretation;
- [E1] full text read; [E2] abstract only; [E3] title or metadata only; [T] textbook.

**Compute:** about 10 CPU-minutes.

| Item | Count / cost |
|---|---|
| New DEVSIM states | 2 (coupled small-signal check) |
| Independent finite-volume solves | 2 |
| Predictor evaluations | 24 × 2 (1-D only) |

---

## 1. Directions considered and the one pursued

| Direction | Decision |
|---|---|
| (a) Predict the connection coefficient B, and hence the subthreshold end correction, from 1-D quantities | **Pursued.** It is the cheapest decisive test of whether the C-10 structure has predictive value |
| (b) Numerical-reliability methodology (network decoupling, variational bounds, cancellation-free currents) | **Partly pursued**, only where it bears on the credibility of results: an independent solver and a coupled small-signal check. Not pursued as a paper topic, because of the substantial prior art the user flagged |
| (c) Return to the S2/S3 numerical problems | Considered. Not pursued: these are implementation-level issues with established remedies |
| (d) Other ideas: a rigorous bound statement for FS93-type integrals (done, §5); a free-boundary analysis of the junction–gate corner; undoped multi-gate bodies, where the electrostatics is linear | (d1) done; (d2) judged high-risk; (d3) **recommended next** (§8) |

## 2. Pre-registered test: two parameter-free predictors

**Test set:** the 24 deep-subthreshold states of C-10 (α∞ ≥ 100):
- 3 oxide thicknesses, 2 dopings, 2 junction depths, ψ_bi shifts, ψ̂ = 12–33;
- ΔL from the source-side current, independently validated in §4.

| Predictor | Errors (all 24 states) | Verdict |
|---|---|---|
| **H-lin**: λ and B from the 1-D linearised transverse eigenproblem (Boltzmann holes, two-layer) and the ε-weighted projection of the junction-face data; ΔL = λ·Ein(B/V_T) | −13.3 % to −24.2 %, median |error| 19.7 %. 0/24 within ±5 %, 24/24 outside ±10 % | **REJECTED** |
| **H-q2D**: classical nonlinear quasi-2-D surface-potential equation (model M7 step 3), with the surface current integral and the field factor | −12.0 % to −35.1 %, median |error| 21.9 %. 0/24 within ±5 % | **REJECTED** |

**Both predictors underestimate the end correction in every state.** No predictor was tested on new devices, because the addendum was conditional on a pass.

## 3. Why the reduced models fail (post hoc; [V] within F0)

**1. The decay length λ itself is predicted correctly.**
- The 1-D eigenproblem gives λ_lin within ≈ 1 % of the depletion-approximation l₀.
- The asymptotic local decay length of the 2-D surface tail (u < 1e-4 V) converges to λ_lin within ≈ 2 %. That residual is consistent with the lateral mesh spacing, h/λ ≈ 0.5.

**2. The amplitude is not.** On the device-2 family, the true far-field amplitude B_∞ is **2.5–5.1 times** the linear projection B_lin. The nonlinear near-junction region amplifies the evanescent mode.
- ψ̂ = 12: 151 V_T vs 29.5 V_T.
- ψ̂ = 30.5: 47 V_T vs 17.7 V_T.
- Device-1 far-tail fits are invalid and not used: at L/2 ≈ 10λ, the opposite tail contaminates the 1e-4 V window.

**3. The zone that sets ΔL is pre-asymptotic.**
- Where u ≈ 0.3–3 V_T, the local decay length of the 2-D tail is 10–35 % above λ_lin. It peaks at about 2λ_lin near u ≈ 7 V_T, then falls towards λ_lin only for u ≲ 0.03 V_T (figure panel a).
- [I] The cause is the nonlinear response of the depletion edge: near the junction the depletion region is deeper, merging with the junction's lateral depletion, so the local scale length is longer.
- A linear eigen-projection cannot represent this. A quasi-2-D equation with fixed W_eff represents it only through Q_s(ψ_s) at the surface, which is too weak.

**4. Correction to C-10.**
- C-10 attributed λ_fit > l₀ (by 2–10 %) to the "soft depletion edge". That attribution was wrong. The soft edge accounts for only ≈ 1 %. The excess comes from the pre-asymptotic fit window (u < 0.5 V_T).
- The C-10 result ΔL = λ_fit·Ein(B_fit/V_T) to ±3 % therefore remains a **descriptive identity with window-fitted parameters**. It is not the leading-order asymptotic law.
- A note to this effect has been added to the C-10 report.

**Consequence.** **No parameter-free reduced model tested here predicts the subthreshold end correction.** A predictive theory would need the nonlinear inner (junction–gate corner) problem: a free-boundary problem (in the depletion approximation, an obstacle-type problem) that has no known closed form.
- The end correction is a function of five dimensionless parameters (C-10).
- At present it can only be tabulated from 2-D solutions.

## 4. Independent check of the numerical formulation

The user cautioned that source-side and swapped-boundary checks are not independent. Two genuinely independent checks were made.

**(i) A separately written 2-D finite-volume solver.** `phase3/c11/fd2d_indep.py`: scipy, its own mesh, Newton solver, homotopy and network, sharing no code with the DEVSIM decks.

| ψ̂ | ΔL network | ΔL slice | ΔL rows |
|---|---|---|---|
| 21 | +0.15 % | −0.11 % | −0.11 % |
| 12 | −0.04 % | −0.13 % | −0.25 % |

- The far-field surface potential agrees to 3.5e-5 V.
- **The C-10 deep-subthreshold values, including the E10-1-corrected ones, are independently confirmed to ≤ 0.25 %.**
- During this check a bug was found and fixed in the new code: its rows bound included rows below h_w that do not touch the source. The table uses the corrected code.

**(ii) DEVSIM's coupled small-signal response** (`devsim_G`: full drift–diffusion Jacobian plus circuit, with its DWR and rounding estimates).

| ψ̂ | G_DEVSIM vs network | DEVSIM's own error estimate (E_dwr) | Condition number κ_G |
|---|---|---|---|
| 21 | −4.6 % | 10 % | 2.4e15 |
| 12 | **wrong sign** | 130 % | 3.8e16 |

- **The coupled small-signal solve is unusable in deep subthreshold**, and its own error estimator says so.
- The equilibrium-decoupled network stays accurate, now confirmed by (i).
- [O] This is a practical TCAD caveat for small-signal analysis at very low currents. Extended-precision options in commercial TCAD suggest it is known in practice. It is not claimed as new.

## 5. The Fjeldly–Shur current integral: assumptions and error bounds

**Statement [P].** The assumptions are:
- zero drain bias (linear response);
- F0 physics (Boltzmann statistics, Einstein relation, holes in equilibrium);
- constant μ;
- the exact 2-D equilibrium potential;
- abrupt Dirichlet junction faces (deck W).

Under these assumptions:
1. **FS93 eq. (2) is the rows bound.** Their eq. (2) applies the one-dimensional current integral along each horizontal line touching the source and sums over those lines. It equals the conductance of the network with all vertical links removed, G_rows.
2. **The Pao–Sah/GCA-type integral is the slice bound.** The sheet-density integral 1/G = ∫dx/σ(x) is the conductance with all vertical lines shorted, G_slice.
3. **Bracketing.** By Rayleigh monotonicity (Thomson/Dirichlet principles; classical), G_rows ≤ G ≤ G_slice, and hence ΔL_rows ≤ ΔL ≤ ΔL_slice.

**Measured bracket width [V], subthreshold (C-10 states, confirmed by §4(i)):**
- ≤ 1.2 % of ΔL;
- 5.4 % at ψ̂ = 12 for device 2.

**FS93's further step** — replacing the integral by L − x_s − x_d with 1-D depletion widths (eqs. 8–9) — **is not a bound**. It underestimates ΔL by a factor of 1.33–3.41 (C-10).

**Novelty.** None in the principle (variational bounds are classical). The value is the explicit error bound for a widely used approximation.

## 6. Literature (C-11 additions)

**Search results [E3].** Targeted searches returned many analytic subthreshold models (double-gate, junctionless, bulk) that use 2-D potential solutions (eigenfunction series, parabolic or quasi-2-D) together with the subthreshold current integral, typically aimed at short-channel barrier lowering. This is consistent with the ingredients being classical.

**What was not found.** No source read gives the long-channel zero-bias end correction in the form λ·Ein(B/V_T). The search was not exhaustive.

**Recommended first read for the proposed next direction:** Liang & Taur, IEEE T-ED 2004, a 2-D series solution for double-gate MOSFETs. Not read in this phase.

## 7. Final scientific decision

1. **Directions considered and pursued:** see §1. Pursued: predictive B (a), and the parts of (b) that affect credibility.

2. **Strongest new finding:**
   - **Model failure.** In bulk MOSFETs, the zone that sets the zero-bias subthreshold end correction lies in a nonlinear, pre-asymptotic transition region, not in the linear evanescent tail. Two parameter-free reduced models therefore fail by 12–35 %: a linear eigen-projection and the classical quasi-2-D equation. The true far-field amplitude is 2.5–5× the linear projection.
   - **Validation.** The C-10 numbers are independently validated (≤ 0.25 %).
   - **Small-signal failure.** The coupled small-signal formulation fails in deep subthreshold, while the decoupled network does not.

3. **Rejected or unresolved hypotheses:**
   - Rejected: H-lin and H-q2D (pre-registered).
   - Unresolved: a predictive theory of B; whether the λ·Ein form exists in the literature.
   - Still standing from C-10: the failed P10-1 and P10-5 tests. The C-10 "soft-edge" attribution is withdrawn (§3.4).

4. **Distinguishable from the literature?**
   - Concept: FS93.
   - Scale length: Frank–Taur–Wong.
   - Quasi-2-D models: Young/Liu-type.
   - Strong-inversion limit: Yau charge sharing.
   - The quantitative failure magnitudes and the pre-asymptotic mechanism appear to be unreported refinements, but they are narrow.

5. **Independent predictive or mathematical value?**
   - **Predictive: no.** No parameter-free model succeeds.
   - **Mathematical: modest.** The bound statement of §5 is a classical principle applied to a known approximation.

6. **What would invalidate the conclusions:**
   - An existing nonlinear reduced model (for example, one with an explicit junction-depletion bulge) that predicts ΔL to ±5 %. That would make the failure finding uninteresting.
   - A published λ·Ein end-correction formula, which would remove the remaining novelty.
   - Errors in the 2-D reference. Now unlikely: independent solver to ≤ 0.25 %, mesh check to ≤ 0.07 %.

7. **Credible publication direction?**
   - **Not for the bulk end-correction line.** I recommend closing it as a publication direction.
   - A numerical-methods paper is conceivable but faces substantial prior art, and I do not recommend it without a demonstrated gap.

8. **What I recommend next:** test whether the end-correction structure becomes predictive where the electrostatics is linear: **undoped double-gate or gate-all-around bodies**.
   - **Hypothesis:** in subthreshold the body obeys Laplace's equation, so the depletion-edge nonlinearity that defeats the bulk models (§3) is absent. The amplitude B should then follow from an eigen-projection of the source/drain boundary data.
   - **Pre-registered prediction to test:** ΔL = λ_DG·Ein(B_proj/V_T), within ±5 %, on a small set of 2-D double-gate solves. This needs a new DEVSIM double-gate deck (≈ 1–2 CPU-hours including verification).
   - **Main risk:** electron screening at the n⁺ face (C-10 measured a ≈ 7 % effect in bulk).
   - **Literature gate first:** Liang & Taur 2004 and its successors. If they already give ΔL in this form, stop.
   - **If it succeeds,** it would be a genuinely predictive, closed-form subthreshold electrical-length result for modern multi-gate devices, suitable for IEEE T-ED or J-EDS after the literature check.
   - **If it fails or is prior art,** end the end-correction programme and choose a new problem.

## 8. Record

**Code** (`phase3/c11/`):

| File | Purpose |
|---|---|
| `bpredict.py` | H-lin |
| `c11_eval.py` | H-lin and H-q2D on the test set; far-tail B_∞ |
| `fd2d_indep.py` | Independent 2-D solver |
| `c11_figures.py` | C-11 figure |

**Data** (`results/phase3/C11/`):

| File or directory | Contents |
|---|---|
| `c11_eval.json` | Predictor errors, B_lin, B_∞, local λ data |
| `indep/` | Independent-solver results |
| `devsim_check/` | Coupled small-signal check: 2 states with logs |
| `c11_figure.png/.pdf` | Figure |

**Failed and corrected work:**
- H-lin and H-q2D rejected.
- The C-10 soft-edge attribution corrected.
- A bug in the new independent code's rows bound, found and fixed before recording.
- The device-1 far-tail fits are invalid (contaminated by the opposite tail) and were not used.
