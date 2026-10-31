# C-9 Research Phase Report — Independent Audit of C-8, Halo Realism, and Originality Reassessment

**Date:** 2026-10-31. **Branch:** `phase3/C9-audit-realism-2026-10-31`, from C-8 `536e1ac`.

**Revision (2026-10-31b).** The first version of this report was committed as `1edc8ad`. Two additions were made after the delivery audit of that version:
- **Erratum E9-1 (§11).** The Lombardi impurity term was wrong for halo devices in probe9 v1. Every halo Lombardi number has been re-computed (probe9 v2, 55 states re-run).
- **Addendum A9 (§12).** Post-hoc mesh and affine-identity diagnostics.

Neither changes a pre-registered verdict. Changed text is marked *(rev. E9-1)* or *(A9)*. The original numbers are kept in §11 and in `results/phase3/C9/analysis/e9_result.json`.

**Labels:**
- **[P]** exact or proved;
- **[C]** conditional analytical identity;
- **[V]** numerical verification within the stated model;
- **[PR]** pre-registered (commit `b83892f`, before the runs);
- **[O]** observation;
- **[I]** interpretation;
- **[PH]** post-hoc, not pre-registered;
- **[E1–E3]** literature level: E1 full text (read directly or through an automated summary of the full text); E2 abstract or chapter page; E3 metadata or title only.

**Device:** N_A = 1e18 cm⁻³, t_ox = 2 nm, N_D = 1e20 cm⁻³, F0 physics, zero bias. **Every length was solved directly in 2-D.** No result in this phase uses affine extrapolation, except where a result is explicitly labelled as the affine identity [C].

**Runs:**

| Set | States | Outcome |
|---|---|---|
| Production, attempt 1 | 121 requested | 77 OK; **44 FAILED** (exit 1, AssertionError): every L = 0.5 µm job. The mesh generator required L > 0.7 µm. Logs preserved in `results/phase3/C9/runs/logs_failed_attempt1/` and recorded in `FAILURES.txt` |
| Fix | — | Shrink the graded lateral end segment to a = 0.35 L for L < 0.7 µm only. Meshes for L ≥ 1 µm are unchanged |
| Production, attempt 2 | 44 | 44/44 OK; the 77 completed jobs were skipped by checkpoint |
| Post-hoc (A9, E9-1) | 85 | All OK: 55 halo states re-run with probe9 v2; 29 mesh-diagnostic states (levels 0 and 1, both probe9 versions); 1 v2 no-halo identity check. Three diagnostic invocations failed before starting because `/usr/bin/time` was not installed; recorded in `diag/DIAG_RUNS.txt` |

About 50 CPU-minutes in total: 25 for production and 25 post-hoc.

**Figure:** `results/phase3/C9/analysis/e9_figure.png` (and `.pdf`), from `phase3/c9/e9_figures.py`. Panels (a)–(f) cover H9-1, H9-3, the strong halo, H9-4, the strong-halo threshold shift and H9-2.

## 1. Strongest finding of this phase

Inside the model:
- the C-8 mechanisms are numerically real;
- the affine law holds on direct 2-D solves down to L = 0.5 µm. This includes subthreshold, the strong-halo device, and field-dependent mobility *(A9)*;
- the equilibrium charge proxy matches the terminal-charge split-C–V within the registered tolerance.

But:
- **Their substance is prior art.** Gate-voltage-dependent L_eff and R_sd from 2-D charge sharing (Guo–Chung 1994, as characterised by the user); multi-length split-C–V correction (Romanjek et al.).
- **Their practical sign depends on halo design.** A strong halo (2e18 on a 1e18 channel) **reverses** the end correction to −550 to −740 nm near threshold. The extraction artefacts then become reverse-short-channel-type, which is well known.

The C-8 extraction-consequences paper **does not have a defensible novelty claim**. §6 gives what remains.

## 2. Audit of the C-8 claims

| C-8 claim | C-9 test | Outcome |
|---|---|---|
| **Affine length law** (ΔL independent of L) | H9-1 [PR]: direct 2-D at L = 0.5, 1, 2 µm, decks W and R | **SUPPORTED.** \|ΔL(0.5) − ΔL(2)\| ≤ 0.10 nm at every V_G ≥ 1.24 V (L/l₀ ≥ 38). [PH] At 1.08 V (subthreshold, outside the registered range): 0.12 nm (W), 0.155 nm (R). The Terada–Muta lines through the direct solves deviate from straight lines by ≤ 1e-4 without halo. With the 2e18 halo the deviation is ≤ 3e-4, for constant and Lombardi μ *(rev. E9-1; the original "≤ 2e-2" was the artefact of §11)*. **Scope:** a property of the 2-D solutions *in this range*, with end-layer overlap decaying like e^{−L/ℓ}; not a theorem for arbitrary L. The C-8 formulas built on it are **conditional identities [C]**: exact given the affine form, and the affine form is verified here |
| **Split-C–V proxy** (C-8 used the equilibrium electron excess) | H9-2 [PR]: terminal charge Q_SD = S/D contact displacement charge + electrons, integrated from V_lo = 1.08 V; L = 1 µm, V_G = 1.32 V | **EQUIVALENT** within the registered 20 %. Deck W: −33.0 % (terminal) vs −29.3 % (electrons); the S/D-boundary displacement charge is 5 % of the measured charge at 1.32 V. Deck R: identical, −33.7 %. **Qualifications:** quasi-static equilibrium charges only (low-frequency limit). No terminal partitioning beyond "electrons belong to S/D". No overlap, interface traps or frequency dispersion. It is a **model-terminal equivalence, not experimental equivalence** |
| **Split-C–V mobility error** (−29 % at L = 1 µm, 1.4e11 cm⁻²) | Direct 2-D at L = 1 and 2 µm, terminal charge | **Confirmed in the model:** −33.7 % (L = 1), −21.0 % (L = 2) at 1.32 V |
| **Two-length cancellation** | Direct deck R, L = 1 and 2 µm, terminal charge, constant μ | **Exact within numerics:** μ error 0.61 % at 1.24 V and ≤ 0.07 % for V_G ≥ 1.32 V. The residual is the charge missed below V_lo. **Not novel:** multi-length differential techniques are established prior art |
| **Linear-extrapolation V_T roll-off** (C-8: −7.43 mV at L = 1 µm, constant μ, extrapolation at the sweep edge 2.0 V) | H9-4 [PR]: maximum-g_m with Lombardi field-dependent mobility, giving an **interior** g_m peak (≈ 1.48 V); direct 2-D, deck W | **SUPPORTED but CORRECTED.** ΔV_T = −3.21 mV (L = 1 µm), −6.89 mV (0.5 µm), −1.56 mV (2 µm), scaling as 1/L (ratio 2.14). It does not change between upper sweep limits 1.8 and 2.0 V. **The C-8 magnitude (−7.4 mV) was an artefact of the extraction definition** (edge extrapolation with constant μ). The 1/L form survives; the magnitude is about 2.3× smaller under the standard procedure |
| **Terada–Muta spurious R_sd** | Regression on direct L = 0.5/1/2 µm | **Confirmed** in sign and order of magnitude. No halo: 23–126 Ω·µm (constant μ) and 87–364 Ω·µm (Lombardi). **2e18 halo: negative R_sd**, −66 to −3084 Ω·µm (constant μ) and −383 to −11290 Ω·µm (Lombardi) *(rev. E9-1; v1: −1340 to −30577)*. **Not novel in substance** (Guo–Chung 1994 per user; Liou et al. 1998 review, E2) |

## 3. Physical realism: halo (pocket) implants

Gaussian pockets: peak 20 nm from each metallurgical junction, σ = 15 nm laterally and vertically, at 15 nm depth. The pocket is truncated at the junction plane; inside the 1e20 box well it would change the doping by < 1 %.

ΔL per end in nm, constant μ, L = 1 µm. It is L-independent within 1.24 nm (0.5, 1 and 2 µm; 0.03 nm for the 5e17 halo). The level-0 → level-1 mesh change is ≤ 4.6 nm, ≤ 1.2 % *(A9)*.

| V_G (V) | No halo | Halo 5e17 | Halo 2e18 |
|---|---|---|---|
| 1.08 | 49.2 | 26.6 | −552.6 |
| 1.24 | 47.5 | 23.5 | −740.9 |
| 1.32 | 39.1 | 18.6 | −609.3 |
| 1.40 | 27.9 | 13.8 | −243.6 |
| 1.48 | 20.5 | 10.4 | −84.5 |

Effect on the maximum-g_m threshold shift (Lombardi μ). Halo rows are re-computed with probe9 v2 *(rev. E9-1)*:

| Device | ΔV_T(0.5 µm) | ΔV_T(1 µm) | ΔV_T(2 µm) |
|---|---|---|---|
| No halo, deck W (H9-4) | −6.89 mV | −3.21 mV | −1.56 mV |
| No halo, deck R | −6.95 mV | −3.25 mV | −1.57 mV |
| Halo 5e17 | −2.80 mV | −1.36 mV | not run (affine identity [C]: −0.67 mV) |
| Halo 2e18 | **+51.6 mV** | **+35.6 mV** | **+22.5 mV** |

The original v1 values were −2.84 / −1.39 mV for the 5e17 halo and +80.5 / +63.7 / +45.9 mV for the 2e18 halo. They are invalid (§11).

**The 2e18 row is reverse-short-channel behaviour and is sub-linear in 1/L** (ratios 1.45 and 1.58 instead of 2).
- *(rev. E9-1, replaces an unsupported [I])* **[C] + [V]** This is the exact affine identity G(L) = σ∞/(L − 2ΔL(V_G)) outside its first-order regime: near the g_m peak, |2ΔL_μ| is not small compared with L (0.20 µm at 1.48 V, i.e. 0.39 L at L = 0.5 µm).
- Evaluated from the L = 1 µm solves alone, the identity reproduces the direct ΔV_T at 0.5 and 2 µm within 0.07 mV *(A9)*.
- **H9-3 [PR] REJECTED.** I predicted that a 5e17 halo would make ΔL negative at 1.32 V. It stays positive (+18.6 nm; +18.9 nm at mesh level 1, A9), roughly halving the intrinsic value.
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
| Taur, IEEE T-ED 2000, DOI 10.1109/16.817582 | E2 (abstract, IBM page) | Shift-and-ratio L_eff; L interpreted through the current-injection points into the S/D; no subthreshold treatment in the abstract |
| Narayanan 1995 | E2 | Electrical L_eff interpretation from 2-D models |
| Cunha, Schneider, Galup-Montoro et al., Nanotech/WCM 2005, pp. 135–138 | E1 (full text via automated summary) *(A9)* | L_eff and W_eff from g_m/I_D in weak-to-moderate inversion. One ΔL per technology; no gate dependence; no S/D-depletion interpretation |
| US 6,275,972 B1 (charge-collection L extraction) | E1 (full text via automated summary) *(A9)* | Subthreshold-biased L extraction; treats charge sharing as an error to minimise; no gate-dependence analysis |
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
| *(rev. E9-1)* This report v1: halo Lombardi ΔV_T (+80.5 / +63.7 / +45.9 mV), Lombardi halo R_sd, "≤ 2e-2" line nonlinearity | **Corrected** to +51.6 / +35.6 / +22.5 mV, −383 to −11290 Ω·µm and ≤ 3e-4 (§11) |
| *(rev. E9-1)* This report v1: "the halo regions' threshold shift is still developing" (non-1/L scaling) | **Withdrawn.** Replaced by the affine-identity explanation (§3, §12) |

## 6. What original contribution, if any, remains

| Remaining element | Evidence | Novelty risk |
|---|---|---|
| (a) Exact zero-bias formulation: conductance = equilibrium-weighted network (P4); Rayleigh two-sided bounds (P5); affine law verified to 1e-5 (C-6) and on direct solves to L/l₀ ≈ 38 (C-9), including halo devices and field-dependent mobility (A9) | [P], [V] | The mathematics is classical. Its use as an exact error decomposition is methodological |
| (b) First quantitative test of Ward 1988 §6.4 | [V], E1 basis. Refuted near threshold; the error budget locates the failure in the phenomenological surface-potential profile; Ward's corner-charge relation holds near threshold on unseen structures | Low risk, but narrow audience |
| (c) **Saturation of the end correction in subthreshold** at ΔL_sub ≈ 1.5–1.6 W_d, the same for two devices 100× apart in doping | **[PR] support** (C-8 E8-1), with an electron-free electrostatic mechanism. [PH] L-independent within 0.16 nm down to L = 0.5 µm (A9) | **Not found in any source read**, including the two weak-inversion L-extraction sources added in A9, which extract a single ΔL and do not interpret it. It is the most physics-like, falsifiable, non-obvious result. It concerns what weak-inversion L extraction measures, which Guo–Chung (strong-inversion L_eff/R_sd extraction, per user) may not address. That is **unverified** without their full text, and the search was not exhaustive |

## 7. Most important unresolved limitations

- Guo–Chung and Romanjek full texts are unread. That is decisive for (c) and for any extraction paper.
- F0 physics: no Fermi–Dirac statistics, no quantum confinement, frozen holes.
- Mobility:
  - Lombardi without distance damping, and only as post-processing (exact for zero-bias linear response, P4).
  - Its impurity term was calibrated on uniformly doped channels. Inside a 15-nm pocket it is an extrapolation.
  - The halo Lombardi results depend on that term: E9-1 changed them by a factor of 1.6–3.6.
- Zero bias only; quasi-static split-C–V; 2-D.
- Halo modelled as Gaussian pockets in an abrupt-junction structure, not process-simulated.
- **Mesh:**
  - No-halo devices: ΔL changes by ≤ 0.2 nm across levels (C-6/C-7; 0.13 nm at 1.32 V here, A9). The level-2 sequence was non-asymptotic (C-6), so the working uncertainty is ≈ 0.5 nm per end.
  - Halo devices (A9): level 0 → 1 changes ΔL by ≤ 1.2 % (≤ 4.6 nm for the 2e18 halo near threshold) and ΔV_T by ≤ 0.2 mV.
  - No conclusion here depends on these.

## 8. Is a publication-oriented campaign justified?

**Not for the extraction-consequences paper.** Its central claims are known in substance, and its practical sign is design-dependent. Not for IEEE T-ED in any current form.

**Conditionally** for one narrower paper:
- **Option 1, physics:** "Electrostatic saturation of the source/drain end correction in subthreshold" (§6(c)), together with (a) and (b) as the framework. Target: SSE or IEEE J-EDS. Justified **only if** a literature gate on subthreshold L_eff extraction (weak-inversion g_m/I_D methods; subthreshold L_eff interpretation) does not already contain it, and a short theory step shows how ΔL_sub/W_d depends on ε_Si t_ox/(ε_ox W_d) and on (ψ_bi − ψ_s)/V_T.
- **Option 2, methodology fallback:** a verification-methodology paper for *J. Computational Electronics*. Content:
  - the zero-bias T1 identity with adjoint-validated error estimates (C-4);
  - the exact network decoupling and bounds;
  - the affine law as a code-verification target.

  This is the most complete and best-evidenced material in the project, with a modest novelty claim. E9-1 is itself an example of what such checks catch: an L-dependence that the affine law forbids exposed an implementation error.

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
| Code | `phase3/c9/probe9.py` (v2 since E9-1: halo, terminal charges, short-L meshing fix, Lombardi impurity fix), `phase3/c9/e9_analysis.py` (`--halo-runs`, `--tag`), `phase3/c9/e9_posthoc.py`, `phase3/c9/diag9_mesh.py`, `phase3/c9/e9_figures.py` |
| Registration | `phase3/specs/C9_Preregistration_2026-10-31.md` |
| Data | `results/phase3/C9/runs/`: 121 v1 state JSONs (kept unchanged), logs, failed-attempt logs, FAILURES.txt, RUNS/RUNNER_END |
| | `results/phase3/C9/runs_v2/`: 55 halo states, v2 |
| | `results/phase3/C9/diag/`: post-hoc mesh diagnostics; see its README |
| | `results/phase3/C9/dev/`: 2 disclosed smoke states |
| Analysis | `results/phase3/C9/analysis/`: `e9_result.json` (v1, original), `e9_result_v2.json` (corrected), `e9_posthoc.json`, `e9_figure.png/.pdf`, stdout records |

**Experiment record:**

| Hypothesis | Outcome |
|---|---|
| H9-1 | SUPPORTED |
| H9-2 | EQUIVALENT (model-terminal) |
| H9-3 | REJECTED |
| H9-4 | SUPPORTED, with the C-8 magnitude corrected |

All four verdicts are unchanged by E9-1: H9-1–H9-3 use constant-mobility or charge data; H9-4 is registered on deck W without halo.

**Failures:**
- Attempt 1: 44 jobs failed on the mesh assertion; fixed and re-run.
- Three diagnostic invocations failed on a missing tool before any simulation started.
- The implementation error of §11.

## 11. Erratum E9-1 — Lombardi impurity term in halo devices (found in the delivery audit)

**Error.**
- probe9 v1 implemented the halo by subtracting it from the donor function: donor9 = N_D − N_halo, so the net doping N_D − N_A − N_halo is correct.
- `endprobe.analyse` builds the Lombardi impurity concentration as N_A + donor(). In v1 that is N_A + N_D − N_halo.
- `endprobe.lombardi` clips it at 1 cm⁻³. Wherever N_halo > N_A + N_D, the C·N^λ term of μ_ac therefore used N = 1 cm⁻³ instead of N_A + N_halo ≈ 3e18 cm⁻³, a factor of about 200 in N^λ. Elsewhere in the pocket, N was under-counted.
- The v1 docstring described this as an "under-count", which misstated its size.

**How it was found.**
1. In the 2e18-halo device, the Lombardi ΔL at L = 0.5 µm differed from L = 1 and 2 µm by up to 169 nm (≈ 6 %). The constant-mobility ΔL agreed to 1.2 nm, and the affine law forbids such an L-dependence.
2. Diagnostic D1 solved L = 1 µm on the L = 0.5 µm end mesh and reproduced the L = 0.5 µm value (−2341.8 nm). So the deviation was a mesh effect, not a length effect.
3. Inspecting what drives a mesh-sensitive Lombardi integral exposed the clipped impurity term.

**Fix (v2, `probe9-2026-10-31b`).** After the equilibrium solve, the Lombardi term receives N_A + N_D + N_halo. The net doping, the equilibrium solution and all constant-mobility quantities are untouched.

**Verification.**
- All 55 halo states were re-run. Their constant-mobility outputs, terminal charges and 1COL references are identical to v1.
- A no-halo state re-run with v2 is identical to v1 (`diag/v2_nohalo_identity/`).
- The v2 halo Lombardi ΔL is L-independent within 1.8 nm (−697.7 / −698.1 / −698.1 nm at 1.32 V for L = 0.5 / 1 / 2 µm). The end-mesh dependence at 1.32 V fell from 5.5 % to 0.05 %.

**Corrected values** (v1 → v2):

| Quantity | v1 | v2 |
|---|---|---|
| 2e18 halo, Lombardi ΔL at 1.08 / 1.24 / 1.32 / 1.40 / 1.48 V (L = 1 µm) | −1964 / −2845 / −2479 / −1096 / −460 nm | −623 / −846 / −698 / −279 / −98 nm |
| 2e18 halo, max-g_m ΔV_T at L = 0.5 / 1 / 2 µm | +80.5 / +63.7 / +45.9 mV | **+51.6 / +35.6 / +22.5 mV** |
| 5e17 halo, max-g_m ΔV_T at L = 0.5 / 1 µm | −2.84 / −1.39 mV | −2.80 / −1.36 mV |
| 2e18 halo, Terada–Muta R_sd (Lombardi) | −1340 to −30577 Ω·µm | −383 to −11290 Ω·µm |
| 2e18 halo, Terada–Muta line nonlinearity (Lombardi) | ≈ 2e-2 | 2.2e-4 |

**Not affected:**
- H9-1–H9-4 verdicts;
- all constant-mobility results, including the halo ΔL table, H9-3 and the constant-μ R_sd;
- all no-halo Lombardi results (H9-4, deck W and R).

## 12. Addendum A9 — post-hoc diagnostics [PH]

Code: `phase3/c9/diag9_mesh.py` and `phase3/c9/e9_posthoc.py`. Output: `results/phase3/C9/analysis/e9_posthoc.json`. Figure panels (c) and (e).

**A9-1. Mesh sensitivity at V_G = 1.32 V** (ΔL per end in nm; Lombardi values from v2):

| Case | ΔL, constant μ | ΔL, Lombardi μ |
|---|---|---|
| 2e18, L = 0.5 µm, end segment a = 0.175 µm, level 0 | −609.35 | −697.68 |
| 2e18, L = 1 µm, a = 0.35 µm, level 0 | −609.29 | −698.06 |
| 2e18, L = 1 µm on the short-device end mesh (a = 0.175 µm), level 0 | −609.35 | −697.68 |
| 2e18, L = 0.5 µm, level 1 | −606.46 | −694.13 |
| 2e18, L = 1 µm, level 1 | −606.42 | −694.30 |
| 5e17, L = 1 µm, level 0 → level 1 | 18.63 → 18.91 | 19.42 → 19.69 |
| No halo, L = 1 µm, level 0 → level 1 | 39.08 → 39.21 | 40.44 → 40.56 |

Over the full 2e18 sweep at L = 1 µm, level 0 → 1 changes ΔL by at most 1.1 % (constant μ) and 1.2 % (Lombardi, v2). Under v1 it was 2.9 %. The short-device end mesh changes ΔL by only 0.06 nm (constant μ) and 0.4 nm (Lombardi) here. The small L = 0.5 µm deviations of §2 are therefore at the mesh level and are not a length effect.

**A9-2. Affine law in subthreshold.** At V_G = 1.08 V, |ΔL(0.5) − ΔL(2)| = 0.12 nm (W) and 0.155 nm (R). This is outside the registered H9-1 range and is reported as an observation. It supports using the subthreshold saturation value (§6(c)) as a length-independent quantity down to 0.5 µm.

**A9-3. Affine identity for the threshold shift [C] + [V].**
- **Method:** V_T(L) is computed from the L = 1 µm solves alone, through G(L) = σ∞(V_G)/(L − 2ΔL(V_G)).
- **Comparison with the direct maximum-g_m ΔV_T:**

| Device | L (µm) | Direct | From the identity |
|---|---|---|---|
| Deck W, no halo | 0.5 | −6.89 mV | −6.88 mV |
| Deck W, no halo | 2 | −1.56 mV | −1.56 mV |
| Deck R, no halo | 0.5 | −6.95 mV | −6.94 mV |
| 5e17 halo | 0.5 | −2.80 mV | −2.80 mV |
| 2e18 halo | 0.5 | +51.64 mV | +51.57 mV |
| 2e18 halo | 2 | +22.51 mV | +22.51 mV |

- **Meaning:** every reported ΔV_T(L) is fixed by one ΔL(V_G) curve. The sub-linear 1/L behaviour of the strong halo is the nonlinearity of the identity itself, not a separate short-channel mechanism.
- **Mesh sensitivity of the strong-halo threshold shift:** the level-1 identity gives +51.6 / +35.8 / +22.7 mV, against +51.6 / +35.6 / +22.5 mV at level 0.

**A9-4. What these diagnostics do not establish.** They test mesh and length consistency inside the F0 + post-processed-Lombardi model. They do not test the physical adequacy of Lombardi inside a pocket, Fermi–Dirac statistics or quantum confinement (§7).
