# C-9 Research Phase Report — Independent Audit of C-8, Halo Realism, and Originality Reassessment

**Date:** 2026-10-31. **Branch:** `phase3/C9-audit-realism-2026-10-31`, from C-8 `536e1ac`.

**Labels:**
- **[P]** exact or proved;
- **[C]** conditional analytical identity;
- **[V]** numerical verification within the stated model;
- **[PR]** pre-registered (commit `b83892f`, before the runs);
- **[O]** observation;
- **[I]** interpretation;
- **[E1–E3]** literature level.

**Device:** N_A = 1e18 cm⁻³, t_ox = 2 nm, N_D = 1e20 cm⁻³, F0 physics, zero bias. **Every length was solved directly in 2-D.** No result in this phase uses affine extrapolation.

**Runs:** 121 states requested.

| Attempt | Outcome |
|---|---|
| 1 | 77 OK; **44 FAILED** (exit 1, AssertionError): every L = 0.5 µm job. The mesh generator required L > 0.7 µm. Logs preserved in `results/phase3/C9/runs/logs_failed_attempt1/` and recorded in `FAILURES.txt` |
| Fix | Shrink the graded end segment for L < 0.7 µm only; meshes for L ≥ 1 µm are unchanged |
| 2 | 44/44 OK; the 77 completed jobs were skipped by checkpoint |

About 25 CPU-minutes in total.

## 1. Strongest finding of this phase

Inside the model:
- the C-8 mechanisms are numerically real;
- the affine law holds on direct 2-D solves down to L = 0.5 µm;
- the equilibrium charge proxy matches the terminal-charge split-C–V within the registered tolerance.

But:
- **Their substance is prior art.** Gate-voltage-dependent L_eff and R_sd from 2-D charge sharing (Guo–Chung 1994, as characterised by the user); multi-length split-C–V correction (Romanjek et al.).
- **Their practical sign depends on halo design.** A strong halo (2e18 on a 1e18 channel) **reverses** the end correction to −550 to −740 nm near threshold. The extraction artefacts then become reverse-short-channel-type, which is well known.

The C-8 extraction-consequences paper **does not have a defensible novelty claim**. §6 gives what remains.

## 2. Audit of the C-8 claims

| C-8 claim | C-9 test | Outcome |
|---|---|---|
| **Affine length law** (ΔL independent of L) | H9-1 [PR]: direct 2-D at L = 0.5, 1, 2 µm, decks W and R | **SUPPORTED.** \|ΔL(0.5) − ΔL(2)\| ≤ 0.10 nm at every V_G ≥ 1.24 V (L/l₀ ≥ 38). The Terada–Muta lines through the direct solves deviate from straight lines by ≤ 1e-4 (no halo) and ≤ 2e-2 (2e18 halo, Lombardi μ). **Scope:** it is a property of the 2-D solutions *in this range*, with end-layer overlap decaying like e^{−L/ℓ}. It is not a theorem for arbitrary L. The C-8 formulas built on it are **conditional identities [C]**: exact given the affine form, and the affine form is verified here |
| **Split-C–V proxy** (C-8 used the equilibrium electron excess) | H9-2 [PR]: terminal charge Q_SD = S/D contact displacement charge + electrons, integrated from V_lo = 1.08 V; L = 1 µm, V_G = 1.32 V | **EQUIVALENT** within the registered 20 %. Deck W: −33.0 % (terminal) vs −29.3 % (electrons); the S/D-boundary displacement charge is 5 % of the measured charge at 1.32 V. Deck R: identical, −33.7 %. **Qualifications:** quasi-static equilibrium charges only (low-frequency limit). No terminal partitioning beyond "electrons belong to S/D". No overlap, interface traps or frequency dispersion. It is a **model-terminal equivalence, not experimental equivalence** |
| **Split-C–V mobility error** (−29 % at L = 1 µm, 1.4e11 cm⁻²) | Direct 2-D at L = 1 and 2 µm, terminal charge | **Confirmed in the model:** −33.7 % (L = 1), −21.0 % (L = 2) at 1.32 V |
| **Two-length cancellation** | Direct deck R, L = 1 and 2 µm, terminal charge, constant μ | **Exact within numerics:** μ error 0.61 % at 1.24 V, ≤ 0.07 % for V_G ≥ 1.32 V. **Not novel:** multi-length differential techniques are established prior art |
| **Linear-extrapolation V_T roll-off** (C-8: −7.43 mV at L = 1 µm, constant μ, extrapolation at the sweep edge 2.0 V) | H9-4 [PR]: maximum-g_m with Lombardi field-dependent mobility, giving an **interior** g_m peak (≈ 1.48 V); direct 2-D | **SUPPORTED but CORRECTED.** ΔV_T = −3.21 mV (L = 1 µm), −6.89 mV (0.5 µm), −1.56 mV (2 µm), scaling as 1/L (ratio 2.14). It does not change between upper sweep limits 1.8 and 2.0 V. **The C-8 magnitude (−7.4 mV) was an artefact of the extraction definition** (edge extrapolation with constant μ). The 1/L form survives; the magnitude is about 2.3× smaller under the standard procedure |
| **Terada–Muta spurious R_sd** | Regression on direct L = 0.5/1/2 µm | **Confirmed** in sign and order of magnitude. No halo: 23–126 Ω·µm (constant μ) and 87–364 Ω·µm (Lombardi). **2e18 halo: negative R_sd**, −66 to −3084 Ω·µm (constant μ). **Not novel in substance** (Guo–Chung 1994 per user; Liou et al. 1998 review, E2) |

## 3. Physical realism: halo (pocket) implants

Gaussian pockets, x_h = 20 nm, σ = 15 nm, at 15 nm depth. ΔL per end in nm; L-independent within 1.2 nm.

| V_G (V) | No halo | Halo 5e17 | Halo 2e18 |
|---|---|---|---|
| 1.08 | 49.2 | 26.6 | −552.6 |
| 1.24 | 47.5 | 23.5 | −740.9 |
| 1.32 | 39.1 | 18.6 | −609.3 |
| 1.40 | 27.9 | 13.8 | −243.6 |
| 1.48 | 20.5 | 10.4 | −84.5 |

Effect on the maximum-g_m threshold shift (Lombardi μ), from the H9-4 analysis:

| Device | ΔV_T(0.5 µm) | ΔV_T(1 µm) | ΔV_T(2 µm) |
|---|---|---|---|
| No halo | −6.89 mV | −3.21 mV | −1.56 mV |
| Halo 5e17 | −2.84 mV | −1.39 mV | not run |
| Halo 2e18 | **+80.5 mV** | **+63.7 mV** | **+45.9 mV** |

The 2e18 row is reverse-short-channel behaviour. It does not scale as 1/L in this range: the halo regions' threshold shift is still developing.

- **H9-3 [PR] REJECTED.** I predicted that a 5e17 halo would make ΔL negative at 1.32 V. It stays positive (+18.6 nm), roughly halving the intrinsic value.
- **[I]** The near-threshold end correction is the sum of intrinsic charge sharing (positive) and halo-induced threshold elevation near the ends (negative). The balance is set by halo dose: weak halos leave a reduced intrinsic effect; strong halos reverse it by an order of magnitude.

**Not implemented** (limits stated, no ad-hoc substitute introduced):
- **Fermi–Dirac statistics.** It needs a generalised Scharfetter–Gummel flux and a generalised Einstein relation in the deck. [D] At the inversion levels where the artefacts are significant (N_s ≤ 5.5e11 cm⁻²), the surface density is ≲ 5e18 cm⁻³ ≪ N_c, so the Boltzmann error is ≲ 6 %.
- **Quantum confinement** (density gradient). Not available in the current DEVSIM deck. Expected to shift the inversion centroid by about 1 nm, which matters for t_ox = 2 nm. Its effect on lateral end electrostatics is untested.

## 4. What the closest prior literature establishes

| Work | Access level | What it establishes |
|---|---|---|
| Guo, Chung & Hsu, IEEE T-ED 1994, DOI 10.1109/16.324592 | **E3** (metadata only; abstract elided by the publisher; no open full text) | User-provided characterisation: connects 2-D charge sharing with L_eff, mobility and R_sd extraction |
| Liou, Ortiz-Conde & García-Sánchez, *Analysis and Design of MOSFETs*, ch. 4 (Springer 1998) | E2 (chapter page) | A review of L_eff extraction methods, including Terada–Muta, shift-and-ratio, Guo–Chung, capacitance methods, 2-D numerical analysis (Narayanan 1995) and gate-voltage-dependent L_eff/R_sd (LDD) |
| Romanjek et al., SSE 2005, DOI 10.1016/j.sse.2004.08.019 | E3 (abstract elided) | Multi-length split-C–V mobility correction |
| Romanjek et al., EDL 2004 | E3 (title only) | Multi-length split-C–V mobility correction |
| Thomas et al., SSE 2009 | E1 (via automated summary) | Low-inversion split-C–V anomalies partly unexplained |
| "Capacitance ratio method … halo implant" | E3 (title only) | L_eff extraction for halo devices is an established topic |
| Taur 2000; Narayanan 1995 | E2 | Electrical L_eff interpretation from 2-D models |
| BSIM3 (Huang et al. 1993) | E1 (via summary) | Exponential roll-off; no 1/L term |

**Conclusion.** Gate-dependent L_eff, extraction bias from 2-D charge sharing, spurious R_sd, multi-length cancellation and halo-induced reversal are all **known in substance**. The exact closed forms (R_sd,TM = 2ΔΔL/Δσ; ΔV_T,lin = (2σ²/(s²L))·dΔL/dV_G) were not found in what I could read. They are first-order algebra on the affine form, and I do not consider them a publishable contribution by themselves.

## 5. Claims that required correction, qualification or rejection

| Claim | Status |
|---|---|
| C-8 linear-V_T magnitude (−7.4 mV at 1 µm) | **Corrected** to −3.2 mV under standard maximum-g_m with field-dependent mobility. The form (1/L) is unchanged |
| C-8 significance for "modern-like" devices | **Qualified.** Halo-free or weak-halo only. Strong halos reverse the sign |
| C-8 split-C–V | **Qualified.** Model-terminal equivalence at the quasi-static level; not experimental equivalence |
| C-8 novelty of the extraction consequences | **Rejected** (§4) |
| H9-3 (weak halo reverses the sign) | **Rejected** by its own test |
| C-7 κ-law failures at α∞ = 0.2 | **Stand as recorded** (unchanged) |

## 6. What original contribution, if any, remains

| Remaining element | Evidence | Novelty risk |
|---|---|---|
| (a) Exact zero-bias formulation: conductance = equilibrium-weighted network (P4); Rayleigh two-sided bounds (P5); affine law verified to 1e-5 (C-6) and on direct solves to L/l₀ ≈ 38 (C-9) | [P], [V] | The mathematics is classical. Its use as an exact error decomposition is methodological |
| (b) First quantitative test of Ward 1988 §6.4 | [V], E1 basis. Refuted near threshold; the error budget locates the failure in the phenomenological surface-potential profile; Ward's corner-charge relation holds near threshold on unseen structures | Low risk, but narrow audience |
| (c) **Saturation of the end correction in subthreshold** at ΔL_sub ≈ 1.5–1.6 W_d, the same for two devices 100× apart in doping | **[PR] support** (C-8 E8-1), with an electron-free electrostatic mechanism | **Not found in any source read.** It is the most physics-like, falsifiable, non-obvious result. It concerns what weak-inversion L extraction measures, which Guo–Chung (strong-inversion L_eff/R_sd extraction, per user) may not address. That is **unverified** without their full text |

## 7. Most important unresolved limitations

- Guo–Chung and Romanjek full texts are unread. That is decisive for (c) and for any extraction paper.
- F0 physics: no Fermi–Dirac statistics, no quantum confinement, frozen holes. Lombardi mobility without distance damping, and only as post-processing (exact for zero-bias linear response, P4).
- Zero bias only; quasi-static split-C–V; 2-D.
- Halo modelled as Gaussian pockets in an abrupt-junction structure, not process-simulated.
- Mesh: ΔL changes by ≤ 0.2 nm across levels (C-6/C-7). The level-2 sequence was non-asymptotic (C-6), so the working uncertainty is ≈ 0.5 nm per end. It does not affect any conclusion here.

## 8. Is a publication-oriented campaign justified?

**Not for the extraction-consequences paper.** Its central claims are known in substance, and its practical sign is design-dependent. Not for IEEE T-ED in any current form.

**Conditionally** for one narrower paper:
- **Option 1, physics:** "Electrostatic saturation of the source/drain end correction in subthreshold" (§6(c)), together with (a) and (b) as the framework. Target: SSE or IEEE J-EDS. Justified **only if** a literature gate on subthreshold L_eff extraction (weak-inversion g_m/I_D methods; subthreshold L_eff interpretation) does not already contain it, and a short theory step shows how ΔL_sub/W_d depends on ε_Si t_ox/(ε_ox W_d) and on (ψ_bi − ψ_s)/V_T.
- **Option 2, methodology fallback:** a verification-methodology paper for *J. Computational Electronics*. Content:
  - the zero-bias T1 identity with adjoint-validated error estimates (C-4);
  - the exact network decoupling and bounds;
  - the affine law as a code-verification target.

  This is the most complete and best-evidenced material in the project, with a modest novelty claim.

## 9. Recommendation

1. A **literature-plus-theory gate on §6(c)**, with no simulation needed first:
   - search the subthreshold L_eff and g_m/I_D extraction literature;
   - derive the electron-free limit ΔL_sub/W_d = F(ε_Si t_ox/(ε_ox W_d), (ψ_bi − ψ_s)/V_T), at least asymptotically;
   - pre-register its prediction for a third structure with a different ε_Si t_ox/(ε_ox W_d), and test it (≈ 3 states).
2. If (c) survives, write the narrow physics paper.
3. If it does not, write the methodology paper (Option 2).
4. **Do not invest further in the extraction-artefact narrative.**

## 10. Artifacts and record

| Kind | Location |
|---|---|
| Code | `phase3/c9/probe9.py` (halo, terminal charges, short-L meshing fix), `phase3/c9/e9_analysis.py` |
| Registration | `phase3/specs/C9_Preregistration_2026-10-31.md` |
| Data | `results/phase3/C9/runs/` (121 state JSONs, logs, failed-attempt logs, FAILURES.txt, RUNS/RUNNER_END); `results/phase3/C9/dev/` (2 disclosed smoke states); `results/phase3/C9/analysis/e9_result.json` |

**Experiment record:**

| Hypothesis | Outcome |
|---|---|
| H9-1 | SUPPORTED |
| H9-2 | EQUIVALENT (model-terminal) |
| H9-3 | REJECTED |
| H9-4 | SUPPORTED, with the C-8 magnitude corrected |

**Failures:** attempt 1, 44 jobs (mesh assertion), fixed and re-run. No other failures.
