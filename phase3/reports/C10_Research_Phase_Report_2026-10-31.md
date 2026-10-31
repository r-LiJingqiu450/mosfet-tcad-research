# C-10 Research Phase Report — The Zero-Bias Subthreshold End Correction against Fjeldly–Shur (1993): Similarity, Bias Range, Oxide Scaling, and a Numerical Issue

**Date:** 2026-10-31. **Branch:** `phase3/C10-subthreshold-theory-2026-10-31`, from C-9 `c41b5bf`.

**Companion theory note:** `phase3/reports/C10_Theory_Subthreshold_EndCorrection_2026-10-31.md`.

**Pre-registration:** `phase3/specs/C10_Preregistration_2026-10-31.md`, commit `76b708b`, before any C-10 run.

**Labels:**
- [P] proved; [C] conditional identity; [V] numerical verification within the F0 model;
- [PR] pre-registered; [PH] post hoc; [O] observation; [I] interpretation;
- [E1] full text read; [E2] abstract only; [E3] title or metadata only; [T] textbook, not re-read.

**Model:** F0, zero bias, deck W unless stated, constant μ, ΔL per end.

**Runs** (≈ 15 CPU-minutes in total):

| Set | States | Outcome |
|---|---|---|
| Pre-registered | 22 | All OK. Job lists `results/phase3/C10/jobs_{sim,sweep,tox,mesh}.json`; runner records in `results/phase3/C10/runs/` |
| Post-hoc network check (`phase3/c10/probe10.py`) | 2 | OK |
| Post-hoc exact-similarity image (`phase3/c10/probe10ni.py`) | 2 | OK |

---

## 1. What this phase found

**1. Fjeldly & Shur (1993) [E1] already contain the concept.**
- Their subthreshold current is the exact current integral (their eq. 2). It is evaluated with L replaced by L − x_s − x_d (their eq. 8).
- x_s and x_d are conventional 1-D abrupt-junction depletion widths (their eqs. 9–10).
- At V_DS → 0 this is the same object as our conductance end correction. Our ΔL lies inside the Rayleigh bracket [ΔL_rows, ΔL_slice], and in subthreshold that bracket is ≤ 1.2 % wide, except 5.4 % at ψ̂ = 12 for device 2 (theory note §1).
- **The existence of a depletion-scale subthreshold end correction is therefore prior art.**

**2. Their quantitative estimate is wrong for the zero-bias long-channel end correction.**
- On all 28 subthreshold states examined (22 new, 6 from C-8), the exact 2-D value is **1.33–3.41 times** their x_s.
- The scaling is also wrong:
  - x_s grows monotonically as V_G falls, whereas ΔL has a broad maximum;
  - x_s does not depend on t_ox, whereas ΔL rises by 30 % from 2 to 6 nm at fixed W_d.
- Fjeldly & Shur call their estimates "simplified and partly empirical". They do not compare x_s with a 2-D solution.

**3. The exact end correction is a scale-length logarithm, not a depletion width [C, V, PH].**
- On every state with α∞ ≥ 100 (24 states: 3 oxide thicknesses, 2 dopings, 2 junction depths, ψ̂ = 12–33),

  ΔL = λ·Ein(B/V_T) ≈ λ[ln(B/V_T) + γ]

  holds to −1.7 to +2.8 % (median +0.1 %).
- λ is the far-field decay length: the 2-D scale length (Frank–Taur–Wong [T]), with λ_fit 2–10 % above the depletion-approximation l₀ because of the soft depletion edge. B is the far-field amplitude of the S/D surface-potential tail.
- In words: the electrical channel starts where the extrapolated tail falls to V_T·e^{−γ} ≈ 0.56 V_T.
- The pre-registered form, which includes a surface-field term r, overestimates ΔL by 1.2–5.7 %. It **failed** its 5 % tolerance on one state.

**4. The C-8/C-9 subthreshold claims do not survive.**
- "ΔL_sub ≈ 1.5–1.6 W_d, the same for two devices 100× apart in doping" is **not universal**.
  - At fixed W_d, ΔL/W_d ranges from 1.31 (t_ox = 1 nm) to 2.38 (6 nm).
  - At equal normalised band bending, the two C-8 devices differ by 14 %. Their equality at equal α∞ was a coincidental compensation between different ψ̂ (21.0 vs 30.5), junction step, junction depth and face carrier density.
- "Saturation" is **not a true limit**.
  - ΔL has a broad maximum near α∞ ≈ 50.
  - It falls 17 % (device 2) and 8 % (device 1) by ψ̂ = 12.
  - The C-8 flatness over α∞ = 15–150 (only 60 mV of surface potential) is the top of that maximum. The scale-length-logarithm structure explains it quantitatively.

**5. The similarity structure is now fully identified [P, V].**
- Five dimensionless parameters govern ΔL/L_D: ψ̂, ψ̂_bi, t̂, ĥ and the junction-face carrier ratio n_face/N_A.
- The pre-registered image omitted the last one, missed by 6–7.5 %, and was **rejected**.
- The exact image (post hoc) agrees to 0.6–0.8 %.
- The two C-8 devices share t̂ exactly. They differ in ψ̂_bi, ĥ and n_face/N_A, and C-8 compared them at different ψ̂.

**6. Numerical issue E10-1, caught by the Rayleigh bounds.**
- The network conductance taken as the drain-side current loses all accuracy in deep subthreshold, by up to 120×. The cause is cancellation in 1 − φ next to the drain.
- The source-side current is cancellation-free:
  - it agrees with an independent swapped-boundary solve to ~1e-14;
  - it satisfies the bounds in every state.
- Historical C-7 to C-9 results are affected by ≤ 3e-4 relative (≤ 0.11 nm without halo). No earlier conclusion changes.

---

## 2. Direct comparison with Fjeldly & Shur (1993)

**Source:** IEEE T-ED 40(1), 137–145, full text supplied by the user [E1]. §§I–IV were read in full (pp. 137–143).

| Aspect | Fjeldly & Shur 1993 | This work (zero bias, long channel, 2-D) |
|---|---|---|
| Definition | Subthreshold current integral, eq. (2). L → L − x_s − x_d, eq. (8) | Conductance end correction; equal to the end correction of the same integral within the Rayleigh bracket (≤ 1.2 % wide; 5.4 % at ψ̂ = 12, device 2) [P/V] |
| End length | x_s = [2ε_s(ψ_bi − ψ_s∞)/(qN_A)]^{1/2}: a 1-D abrupt-junction depletion width, eq. (9) | ΔL = λ·Ein(B/V_T): a 2-D scale length times a log of the tail amplitude (±3 % on 24 states) [PH] |
| Magnitude (28 states) | x_s | ΔL/x_s = 1.33–3.41 |
| Bias trend | x_s grows monotonically as V_G falls | Broad maximum near α∞ ≈ 50; −17 % (device 2) by ψ̂ = 12 |
| Oxide | x_s does not depend on t_ox | ΔL ∝ λ^{≈1.4} at fixed W_d: +30 % for t_ox 2 → 6 nm |
| 2-D length scale | λ (quasi-2-D, their eq. 22), used only for DIBL (eqs. 21–31) | λ governs ΔL itself |
| Validation | DIBL ΔV_T against measured data with fitted χ, σ. No 2-D check of x_s | Direct 2-D solves, mesh-verified (level 0 → 1: ≤ 0.07 %) |

**Verdict.** FS93 establishes the concept and the depletion scale. The defensible addition here is quantitative: the correct structure and magnitude of the zero-bias long-channel end correction, and the size of the error of the 1-D depletion-width estimate. That is a refinement of a 1993 compact-model approximation, not a new phenomenon.

---

## 3. Pre-registered tests

| Test | Result | Verdict |
|---|---|---|
| **P10-1** similarity image (t̂, ĥ, ψ̂_bi matched) | ΔL/L_D ratio vs device 1 = 1.075 (ψ̂ = 21.02) and 1.063 (16.0) | **REJECTED** (limit ±2 %). Cause identified post hoc: n_face/N_A (§4) |
| **P10-2** factor decomposition at ψ̂ = 21.02 (descriptive) | ΔL/L_D vs device 1 (10.18): A, device 2 as is, +13.7 %; B, ψ̂_bi matched, +8.8 %; C, ĥ matched, +11.8 %; D, both, +7.1 %; [PH] n_i image (all matched) +0.4 % | Most of the gap comes from the junction step and face screening. Junction depth ĥ 123 → 12 contributes only −1.6 to −1.9 points |
| **P10-3** oxide scaling at fixed W_d | 6 nm / 2 nm = 1.309 (ψ̂ 30.50), 1.302 (21.02); l₀-ratio 1.207/1.220; λ_fit ratio 1.175/1.181 | **H-l0 SUPPORTED, H-Wd REJECTED** |
| | 1 nm / 2 nm = 0.825 and 0.840; l₀-ratio 0.928/0.919 | Reported only: below both registered intervals at ψ̂ = 30.5. The oxide dependence is stronger than λ alone, because B also grows with t̂ |
| **P10-4** saturation, device 2 | ΔL(ψ̂ = 12)/max = 0.834; ΔL(12) is −16 % from ΔL(α = 150) | **H-max SUPPORTED, H-sat REJECTED** |
| **P10-4** saturation, device 1 | 0.923 (criterion ≤ 0.92); −7.3 % from α = 150 | **INCONCLUSIVE** by the registered rule. H-sat fails, and H-max misses its threshold by 0.3 % |
| **P10-5** asymptotic identity with the r term | Errors +1.2 to +5.7 % (median +2.7 %) | **REJECTED** at the 5 % tolerance (one state at 5.72 %). Post hoc without r: −1.7 to +2.8 % (§4) |
| **P10-6** FS93 eq. (9) as ΔL | ΔL/x_s = 1.33–3.41 | **H-FS REJECTED** (as stated in advance) |
| Mesh, level 0 → 1 | −0.02 % (device 1), −0.07 % (device 2) | Discretisation error negligible |

**Disclosure on the measured quantity (E10-1, §5).** The registration named `dL_net_nm`, the drain-side network current. In 16 of the 22 new states that value violates the Rayleigh bounds and is numerically invalid; the worst is ×12 at ψ̂ = 12.

All verdicts above use the cancellation-free source-side value of the same physical quantity. The registered drain-side values are kept in `results/phase3/C10/analysis/c10_result.json` as `dL_drain_prereg`, with their bound flags. With the drain-side values, P10-4 for device 2 would also read H-max, through spurious ×12 values. No verdict is rescued by the correction.

---

## 4. Post-hoc findings (not pre-registered)

**Exact similarity image** (`probe10ni.py`).
- The device-1 image uses N_A = 1e18 cm⁻³, n_i = 1e12 cm⁻³ and N_D = 1e22 cm⁻³, so that every dimensionless parameter of device 1 is matched, including n_face/N_A = 1e4. This is a mathematical test, not a device.
- It reproduces device 1's ΔL/L_D to +0.8 % (ψ̂ 21.02) and +0.6 % (16.0).
- The P10-1 failure is therefore fully explained by the omitted face-screening parameter.
- **Physical meaning:** in deck W (Ward truncation), electrons next to the junction face screen part of the junction step, and how much depends on N_D/N_A. In deck R the n⁺ region is resolved instead.

**Identity without the r term.**
- λ_fit·Ein(B/V_T) matches ΔL to −1.7 to +2.8 %, median +0.1 %, over all 24 states with α∞ ≥ 100.
- κ_fit − ln(B/V_T) = 0.44–0.65, median 0.574 (γ = 0.577).
- The surface-field factor r and the near-junction deviation of the tail roughly cancel. This is an empirical observation, not a derivation.

**E8-1 explained.**
- Over α∞ = 15–150, ψ_s changes by only 60 mV.
- l₀ shrinks by 3–4 % (Δκ ≈ +0.14), and ln(B/V_T) grows by ≈ 0.1–0.2.
- The structure therefore predicts Δκ ≈ +0.3; observed +0.33/+0.40.
- The C-7 "H-log" law (κ ≈ ln α∞ + 1.6) failed because ln α∞ is linear in ψ_s, whereas the true dependence is logarithmic in ψ_bi − ψ_s.

**Connection coefficient is not universal.**
- β = B/(ψ_bi − ψ_s∞) ranges from 1.3 to 4.7.
- It rises with t̂ (1.3–1.6 at 1 nm, 2.3–3.0 at 2 nm, 3.8–4.7 at 6 nm) and is lower for device 1. It has a shallow minimum in ψ̂.
- An a-priori theory of B is the missing step for prediction. Because ΔL depends on B only through a logarithm, a 20 % error in B costs ≈ 0.2λ (≈ 5 %).

---

## 5. Numerical issue E10-1: cancellation in the drain-side network current

**Error.**
- `phase3/c6/wprobe.solve_net` fixes φ = 0 on the source and φ = 1 on the drain.
- It returns G as the drain current Σ c_ij(1 − φ_j). Next to the drain, 1 − φ_j is far below the double-precision spacing near 1, while c_ij spans more than 30 orders of magnitude in deep subthreshold.
- The source current Σ c_ij φ_j (φ_j ≈ 0, all terms positive) has full relative accuracy.

**Detection.** The pre-registered Rayleigh bracket G_rows ≤ G ≤ G_slice was violated in 16 of the 22 new states, by up to 120× (ΔL 496 nm against a bracket of 40.7–43.0 nm at ψ̂ = 12).

**Fix and verification** (`phase3/c10/probe10.py`).
- G is taken as the source-side current.
- An independent swapped-boundary solve agrees to ≤ 3e-14 (ψ̂ = 12 and t_ox = 6 nm states, both networks).
- Every corrected value lies inside its bracket.

**Historical audit** (all 284 state records of C-7 to C-9): |G_drain/G_source − 1| ≤ 3e-6 (C-7), ≤ 4e-5 (C-8), ≤ 3e-4 (C-9). The largest ΔL change is 0.11 nm without halo and 0.75 nm on the −2021 nm Lombardi halo value. No earlier conclusion changes.

**Lesson.** Together with E9-1, this is the second implementation error in two phases that the exact-identity checks (affine law, Rayleigh bounds) exposed.

---

## 6. What remains, and what does not

**Rejected or withdrawn from earlier phases:**
- the "1.5–1.6 W_d" universality;
- "saturation" as a limit;
- the claim that the subthreshold end correction is new in concept (prior art: FS93).

**Remaining candidate contributions, with novelty status:**

| Candidate | Evidence | Novelty risk |
|---|---|---|
| (a) The scale-length-logarithm law ΔL = λ·Ein(B/V_T) for the zero-bias subthreshold end correction, with the quantitative failure of the 1-D depletion-width estimate (×1.3–3.4, wrong t_ox and bias scaling) | [PH/O] on 24 states; [C] derivation | **Moderate–high.** The ingredients (scale length, exponential tails, the subthreshold integral) are classical. Whether post-1993 subthreshold models already evaluate the integral this way (exponential-integral forms) is **unchecked**. The connection coefficient B is not yet predicted |
| (b) Certified evaluation of zero-bias conductance and end corrections in drift–diffusion TCAD | Exact network form (P4); Rayleigh brackets as a-posteriori certificates (P5); cancellation-free near-zero-side currents with a reciprocity check (E10-1); the affine law as a code-verification target; two documented real failures (E9-1, E10-1) caught | **Moderate.** Rayleigh bounds are classical in network and homogenisation theory. Their use as certificates for TCAD small-signal quantities, and the contact-current cancellation, need a literature gate (weighted-residual terminal-current methods, small-signal analysis literature) |
| (c) The five-parameter similarity structure, including junction-face screening | [P] + [V] to 0.6–0.8 % | **Low as a standalone contribution** (dimensional analysis). Useful inside (a) |

---

## 7. Literature status

| Work | Level | Relevance |
|---|---|---|
| Fjeldly & Shur, IEEE T-ED 40, 137 (1993) | **E1** (user-supplied PDF, §§I–IV) | Same concept; 1-D x_s estimate, off by ×1.3–3.4 here; quasi-2-D λ used only for DIBL; DIBL validated against experiment with fitted parameters |
| Guo, Chung & Hsu, IEEE T-ED (1994), DOI 10.1109/16.324592 | **E3.** The abstract was not retrievable by my tools; the publisher page and an NCTU repository record rendered without content | User characterisation only |
| "Extraction of MOSFET threshold voltage, series resistance, effective channel length, and inversion layer mobility from small-signal channel conductance measurement", IEEE T-ED (2001), DOI 10.1109/16.974720 | **E3** (title from IEEE Xplore and search) | Establishes zero-dc-drain-bias small-signal conductance as an extraction route: the measured quantity is our G(L) at zero bias |
| Chen et al., IEICE Trans. Electron. E87-C, 726 (2004) | E2 (repository abstract) | S-parameter, zero-drain-bias extraction of V_T, L_eff and μ |
| Frank, Taur & Wong, IEEE EDL (1998), generalised scale length | [T] | Our l₀ = λ_FTW/π |
| Cunha et al., Nanotech (2005); US 6,275,972 | E1 via summary (C-9) | Weak-inversion L extraction: a single ΔL, with no interpretation of its gate dependence |

---

## 8. Limitations

- F0 physics: Boltzmann statistics, no quantum confinement, constant μ in the subthreshold analysis.
- Zero drain bias, 2-D, deck W (Ward truncation).
  - Deck R differs by 1–2 % at threshold (C-9).
  - The face-screening parameter shows that the junction model matters at the 7 % level.
- Abrupt box junctions. Halos and LDD dominate in real devices (C-9: a strong halo reverses ΔL).
- B is measured, not predicted.
- The post-hoc identity (±3 %) is empirical in its constant. The derived form with r overestimates.
- The tail fits use 7–16 points per state, with maximum log residuals 0.007–0.063.
- The literature gate for candidates (a) and (b) beyond FS93 has not been done.

---

## 9. Is a publication campaign justified?

- **Not yet, and not for IEEE T-ED** in the present form.
- The subthreshold result is a quantitative correction of a 1993 approximation within an idealised model. Its practical significance is limited by real junction profiles.
- **Conditionally justified** for a journal such as *J. Computational Electronics* or IEEE J-EDS, if a literature gate clears it, as **one** paper that combines:
  - (b), certified zero-bias TCAD quantities, as the methodological core, with E9-1 and E10-1 as documented failure cases the certificates caught;
  - (a), the subthreshold end-correction law with the FS93 comparison, as the physical application.

## 10. Recommendation for the next phase

1. **Literature gate first** (no simulation):
   - Do post-1993 subthreshold or short-channel models evaluate the current integral with exponential tails (exponential integral, "Ein", logarithmic end lengths)? This covers Liu et al. 1993, conformal-mapping subthreshold models (Klös and co-workers) and later compact models.
   - Is certified, bounded, or cancellation-free terminal-current evaluation known in TCAD (weighted-residual currents, small-signal analysis literature)?
2. If (a) survives, attempt an a-priori model for B, e.g. an inner nonlinear corner problem solved once per (t̂, ψ̂_bi, n_face/N_A). Test it, pre-registered, on unseen structures including deck R and a graded junction.
3. If only (b) survives, write the methodology paper with the end-correction physics as its validated example.
4. If neither survives, stop the end-correction line and return to the project's open numerical issues (S2/S3 low-current cancellation, conservation, stopping criteria). E10-1 shows that the certification tools developed here apply directly to those issues.

---

## 11. Record

**Failed, corrected, rejected and inconclusive work:**
- P10-1 rejected;
- P10-5 rejected at its 5 % tolerance;
- P10-4 inconclusive for device 1;
- H-sat rejected for both devices;
- H-Wd and H-FS rejected;
- E10-1: an implementation error, corrected and audited;
- three tooling-only failures (C-9, earlier);
- no failed 2-D run in C-10.

**Code:** `phase3/c10/` — `c10_design.py`, `tail_analysis.py`, `c10_analysis.py`, `c10_figures.py`, `probe10.py` (cancellation-free conductance), `probe10ni.py` (exact similarity image).

**Data:**
- `results/phase3/C10/`: `runs/` (22 states with runner records), `netcheck/` (2 states and diagnostics), `ni_image/` (2 states);
- `tail/` (tail fits of C-7, C-8 and C-9 states);
- `analysis/`: `c10_result.json` and `c10_figure.png/.pdf`.
