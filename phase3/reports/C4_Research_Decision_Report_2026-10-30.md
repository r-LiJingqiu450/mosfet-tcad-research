# C-4 Research Decision Report — Numerical Reference Closure and Novelty Decision

**Date:** 2026-10-31 (UTC). **Branch:** `phase3/C4-2026-10-30`, from verified C-3 head `e4e4528`.

Reference validity and novelty are decided **separately**.

| Decision | Status |
|---|---|
| **Reference status** | **`REFERENCE_READY`** for the zero-bias identity and the boundary/deck verification of R-I(a), R-I(b) (corrected) and R-I(c) (f = 0.1, 0.3, 1). Scope and stated limits in §1 |
| **Novelty status** | **`NOVELTY_UNPROVEN`**. Several sub-claims are rejected or downgraded, and D1 as originally framed is withdrawn (§3) |
| **Limited finite-bias stop/go (48-state FB-1)** | **Not recommended in its registered form** (§5). A reformulated, smaller test is proposed for separate authorization after independent review |
| Package E, parameter sweep, manuscript | Not started |

---

## 1. Why `REFERENCE_READY`, and its exact scope

**What is established:**
- **T1 v2** (registered before execution): 192/192 PASS on both independent paths.
  - The paths are equilibrated SuperLU and DEVSIM AC with an L2-native estimate.
  - Six deck variants, L = 1, 3, 10 µm, levels 0–3, V_G = 0.8 / 1.25 / 2.0 V. No FAIL.
- **The C-3 class I point is explained and resolved:**
  - Cause: a non-backward-stable unscaled solve (η = 9.4e-3).
  - v2: r₀ = −1.75e-8 with E = 8.4e-8 (L2v2), and r₀ = −2.4e-7 with E = 6.0e-7 (L1v2).
  - The historical v1 verdict is preserved.
- **The error estimates are validated, not assumed.** The equilibrated-SuperLU − UMFPACK difference lies within the sum of their estimated errors at 192/192 points (max ratio 0.93), with the correct sign at 192/192. Backward error ≤ 1.6e-15 for every solver used in a verdict.
- **R-I(b) bottom corners are corrected.**
  - Both R04 §4.2 conditions are imposed on the corner nodes, with no node shared by two contacts.
  - The lateral ψ spread falls from 1.7e-11 V to 1e-16 V.
  - The effect on G is below resolution (≤ 0.59 of the combined linear-solve estimate).
- **R-I(c) is implemented, verified and T1-tested** (84/84).

**Stated limits:**
1. **Conservative verdict.** Under the worst-case rounding term (E_cond = κ_G·u), all 18 points at L = 10 µm, level 3 (every deck) are NON-RESOLVING. The primary rounding term (RSS) was designated before execution and is the one validated against cross-solver differences. **A reader who adopts the worst-case term should read the status as `REFERENCE_PARTIAL` for L = 10 µm at level 3.** The worst case is reported, not hidden.
2. **B5v2 (≤ 1e-8 linear-solve part) is not attainable** for L ≥ 3 µm on finer meshes. κ_G reaches 1.2e10, so the error is conditioning-limited, not solver-limited. T1's E ≤ 1e-6 is met.
3. **T1 is a zero-bias identity.** Finite-bias correctness rests on:
   - conservation (δ_KCL, δ_cont ≤ 3e-14);
   - asymptotic mesh convergence (ℓ_res diagnostic: 5 of 6 references asymptotic);
   - agreement of the same-mesh and continuum metrics (≤ 0.2 %).
4. **What READY does and does not mean.** It means the decks are correct implementations of their **specified** boundary-value problems. It does **not** mean any of them is a physically correct reference for real devices (§2).

## 2. Task F — scientific interpretation of the opposite-sign discrepancies

| Question | Answer | Evidence |
|---|---|---|
| Numerically resolved? | **Yes**, at L = 1 µm, V_G = 1.25 V, V_D = 50 mV | ε^a = −5.126e-4 ± 3.5e-6; ε^b = +5.703e-4 ± 2.9e-6; asymptotic; same-mesh and continuum metrics agree ≤ 0.2 % |
| Sensitive to the reservoir definition? | **Yes, dominantly.** Admissible ideal reservoirs give ε from −5.1e-4 to +5.7e-4 at the same point | Within the gated-strip family ε falls from +5.70e-4 (f → 0) through +2.42e-4 and +6.8e-5 to an unresolved value at f = 1 |
| Explained by known asymptotic theory? | **Structure partly explained by the project's own perturbation analysis; Ward unknown.** | (i) O(V_D) exact for any φ_n-Dirichlet reservoir (T1). (ii) O(V_D²) in current exact for mirror-symmetric reservoirs, hence ε ∝ V_D² (Prop. 2). (iii) The leading coefficient has same-order (L⁻²) interior and reservoir-dependent parts (Prop. 3, a scaling argument). Ward 1988 §6.4 and 1990 could not be read |
| Robust under R-I(c)? | **No.** The R-I(a)/R-I(b) values are not robust. They are end-condition truncation effects of the reference | ε^c(f = 1) is ~100× smaller than ε^a or ε^b. Long-strip convergence (f = 3) is non-asymptotic, so open |
| Reservoir-symmetric quantity | **Not a reservoir-independent reduction error.** σ_ab = +2.9e-5 ± 3.2e-6 differs from ε^c(f = 1) = +1.2e-6 ± 3.6e-6 | Prop. 3; ℓ_res data. Three or four tested reservoirs never establish independence, and these show dependence |
| Scientifically useful beyond the literature? | **Undetermined.** Possibly as a verification/methodology statement (§3) | Literature gate not passed |

**Mechanism test TC-1** (registered verdict INCONCLUSIVE):
- Mirror-asymmetric reference pairs produce an O(V_D) discrepancy, c₁ = ±0.0336 V⁻¹ at L = 1 µm, scaling as L^−2.00, with c₁(ab) = −c₁(ba).
- Symmetric pairs do not, to 5–7×10⁻⁶ of that. But they sit at 3.1–4.0σ of a truncation-dominated fit, outside the registered 3σ bound.
- The mechanism of Proposition 2 is therefore strongly indicated **post hoc** but **not confirmed by the registered criterion**.

## 3. Task E and the novelty decision

**Ward full-text audit: not performed.**
- The 1990 SIAM paper is behind a signed-URL redirect the proxy rejects.
- The thesis text extraction stops at p. 51 of ~150 (`C4_Ward_FullText_Audit_2026-10-30.md`).
- No content of those works is asserted.

**Claim decisions** (`C4_Novelty_Assessment_Update_2026-10-30.md`):

| Claim | Decision |
|---|---|
| K3a (δ² law) | **Rejected as novel**: an elementary symmetry consequence |
| K3b (α ∝ L⁻²) | **Rejected as novel by itself** |
| K2 (corner/contact treatment) | **Rejected as a scientific contribution** |
| K1 (verification test with validated DWR estimate) | **Downgraded** to methodology |
| K4 (same-mesh metric) | **Downgraded** to methodology |
| K3c (end-condition dependence, quantified) | **Unproven**: Ward §6.4 is the direct competitor |
| K5, K6 | **Unproven, premises weakened** |
| K7 (new: mirror-symmetry selection rule) | **Unproven** |
| **D1 as originally framed** | **Withdrawn**: its central quantity is not defined independently of an arbitrary reservoir specification at leading small-bias order |

**Overall: `NOVELTY_UNPROVEN`.**
- Not `REJECTED`: no publication has been shown to establish K3c or K7.
- Not `POTENTIAL`: there is no literature comparison.

## 4. Proposed scientific hypothesis and falsification criteria (Deliverable 8)

**H-R (reservoir truncation).**
- **Setting:** idealized 2-D references with perfectly conducting (φ_n-pinned) gated reservoirs, F0 physics, small V_D.
- **Statement:** the small-bias GCA discrepancy coefficient c^c(ℓ) = ε^c/V_D² converges, as ℓ_res/W_dep → ∞, to a limit c_∞(V_G, L) with |c_∞| ≪ |c^a|, |c^b|.
- **Corollary:** the R-I(a)/R-I(b) discrepancies are end-condition truncation effects of the references, scaling as L⁻². They are not a reduction error of the GCA.

H-R is **falsified** if, in a newly registered test with a fixed strip cell size and f ∈ {1, 2, 3, 5}, levels 0–3, L ∈ {1, 2, 4} µm, V_G ∈ {1.0, 1.25, 2.0} V, V_D ∈ {8, 50} mV, any of the following holds:
- **F1:** the successive differences ε^c(f_{k+1}) − ε^c(f_k) do not decrease monotonically in f (beyond 2·U_disc); or the f = 3 → 5 change exceeds 2·U_disc.
- **F2:** the extrapolated |ε_∞| ≥ 0.1·min(|ε^a|, |ε^b|) at any resolved (L, V_G), i.e. the long-reservoir discrepancy is not small compared with the truncated references.
- **F3:** |ε^a − ε_∞| (or |ε^b − ε_∞|) does not scale as L^−2 ± 0.3 across the three lengths, which would contradict Prop. 3's end-layer scaling.

**Why H-R would matter, if it survives and the literature does not already contain it:**
- Reduced-model (GCA/compact) validation against idealized 2-D references can be dominated by the reference's reservoir truncation.
- A symmetry rule (TC-1) and a convergence criterion in ℓ_res give a practical check.
- This is a **verification-methodology** statement, not a new device model. Its suitability for T-ED, as opposed to a numerical-methods or TCAD-verification venue, is doubtful and must be judged after the literature gate.

**Companion mechanism hypothesis H-S** (symmetry rule; re-test of TC-1 under a truncation-robust rule):
- Statement: |c₁,sym| ≤ 10⁻³·|c₁,asym| at every level and L, with a quartic fit and δ down to 0.25 mV.
- c₁,asym ∝ L^−2 ± 0.1 over L ∈ {1, 2, 3, 4} µm.
- c₁(ab) + c₁(ba) = 0 within 2σ.
- Falsified by any one of these failing.

## 5. Recommendation on the limited finite-bias stop/go experiment (Deliverable 10)

**Do not authorize the registered 48-state FB-1 stop/go.**
- Its decision quantity (H-C: a reservoir-independent σ) is shown to be ill-posed:
  - Prop. 3;
  - σ_ab ≠ ε^c(f = 1);
  - the 100× spread across admissible reservoirs.
- Running it would produce a well-converged number with no defensible interpretation.

**Instead, in this order:**
1. **Literature gate (required).** Ward 1990 and thesis ch. 5–6 attached and audited equation by equation. Then:
   - if §6.4 already contains the end-condition dependence or the symmetry rule, **reject K3c/K7** and redirect;
   - if not, proceed.
2. **Independent review** of C-4 (R05 rule 6): T1 v2 registration and results, Props. 1–3, TC-1, the ℓ_res diagnostic.
3. **Only then, and with separate approval:** a newly registered, minimal test of H-R and H-S (§4).
   - About 40–60 states at levels 0–3.
   - Fixed strip cell size.
   - Estimated 3–5 CPU-hours on this machine, based on measured C-4 timings: level-3 strip states take ~2.5–3 min per bias point.
4. **Strategic alternative (needs Package E authorization):** if the goal remains a T-ED-level physical contribution, move the reference to a physically structured device (R-II, n⁺ source/drain). There the "GCA error" is defined by the device rather than by an idealized reservoir. C-4 is the main argument for doing so.

## 6. Deliverables index

| # | Deliverable | Location |
|---|---|---|
| 1 | C-4 numerical reference report | `phase3/reports/C4_Numerical_Reference_Report_2026-10-30.md` |
| 2 | Complete revised T1 matrix and raw error data | `results/phase3/C4/t1v2_{B,C,D,Dprime}/` (64 JSON + logs + RUNS.txt); `results/phase3/C4/analysis_t1v2/t1v2_table.csv` (192 rows, every error component), `t1v2_summary.json`, `v1_vs_v2.csv`, `T1v2_overview.png` |
| 3 | Independent L1 and L2 error assessments | Numerical report §2.3–2.4; per-point `eval` blocks (η, residual, DWR, E_round, E_cond, κ_G, sections) in every T1 v2 JSON |
| 4 | Corrected R-I(b) deck and evidence | `phase3/c4/deck4.py` (model `b`); numerical report §3; `t1v2_C/` |
| 5 | R-I(c) specification and T1 results | `phase3/specs/R-Ic_Reservoir_Strip_Model_Specification_2026-10-30.md`; numerical report §4; `t1v2_D/`, `t1v2_Dprime/`; ℓ_res: `C4_lres_and_TC1_Results_2026-10-30.md` |
| 6 | Full-text Ward analysis | `C4_Ward_FullText_Audit_2026-10-30.md` — **not performed (access blocked); attempts documented** |
| 7 | Updated claim-vs-prior-art assessment | `C4_Novelty_Assessment_Update_2026-10-30.md` |
| 8 | Hypothesis with falsification criteria | §4 above; derivation in `C4_Small_Bias_Perturbation_Analysis_2026-10-30.md` |
| 9 | Source ZIP, Git bundle, SHA-256 manifests | Archive files (final HEAD SHA in `bundle_verification.txt`; a commit cannot contain its own SHA) |
| 10 | Stop/go recommendation | §5 above |

**Pre-registrations:**
- `phase3/specs/C4_T1v2_Preregistration_2026-10-30.md` (`409e521`);
- TC-1 rules in the perturbation analysis (`eb8309b`).

**Run infrastructure:**
- `phase3/c4/jobrunner.py`: checkpointing, bounded, failure-reporting. It replaced the open-ended wait chains after the mid-run diagnosis.
- No simulation was restarted or discarded; the stopped processes were wait loops only.
