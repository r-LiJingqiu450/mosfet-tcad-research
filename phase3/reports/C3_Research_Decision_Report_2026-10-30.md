# C-3 Research Decision Report — Scientific Reference Validation and Novelty Gate

**Date:** 2026-10-30. **Branch:** `phase3/C3-2026-10-30` (base `c909bf11a0c6b4a41dcef613be7a4bb0d4a29b61`, C-2, provenance verified at start).

| Decision | Status |
|---|---|
| **Reference classification** | **`REFERENCE_PARTIAL`** |
| **Novelty status** | **`NOVELTY_UNPROVEN`** |
| **Limited finite-bias discrepancy experiment** | **Not yet scientifically ready.** The design (FB-1) is complete and falsifiable. Its 48-state stop/go stage becomes runnable once three preconditions are met (§5). Not started |
| Package E, large MOSFET sweeps, manuscript | Not started |

---

## 1. Workstream outcomes

| WS | Task | Outcome | Evidence |
|---|---|---|---|
| WS1 | Preservation | Verified git bundle, source ZIP and SHA-256 manifests (file-level and archive-level). **No push.** GitHub write access was denied in an earlier phase and was not retried. `main` (`c64c463`) and all remote branches are untouched | Archive files `bundle_verification.txt`, `SHA256SUMS.txt`, `MANIFEST_SHA256_files.txt` |
| WS2 | C-1(b) final acceptance | B1–B6 decided. The original 1e-10 figure is withdrawn and recorded as **not met** (1.4–1.9e-10). **Addendum A:** on the T1 matrix B5 is met on 41/72 points, and B3 is not met for the T1 right-hand side of the L1 path (2 points examined). C-1(b) is ACCEPTED for identification and Newton solves, and ACCEPTED WITH LIMITATION for L1 linear response | `phase3/reports/C1b_Final_Acceptance_Review_2026-10-30.md` |
| WS3 | Formal T1, R-I(a) | 35/36 PASS for L1 and for L2. **1 NON-RESOLVING (class I):** L = 10 µm, level 3, V_G = 2.0 V. 0 FAIL. Post-hoc diagnostic: caused by the unscaled SuperLU solve of the L1 path (equilibrated backward error 9.4e-3). DEVSIM's own UMFPACK step agrees with the identity to 1.7e-8. The verdict is not reclassified | `phase3/reports/C3_T1_RIa_RIb_Report_2026-10-30.md` §2–§5 |
| WS4 | R-I(b) boundary verification and T1 | 36/36 PASS. Ownership, interface rows, continuity (0 V), displacement flux (Gauss closure ≤ 4.3e-14), conservation, refinement and the absence of artificial constraints were checked from the exported matrix. **Open item R-I(b)-1:** the two bottom corner nodes do not carry the body ψ condition. Estimated impact on G ≤ 7.5e-12 | same report, §6–§7 |
| WS5 | Specification reconciliation | Research_08 v2.2: T0/T2 factor of two, D-1 grading, T1 meaning, B1–B6, triple-point rule, 1COL admissibility, S_res wording. Old and new wording side by side; historical documents unchanged. **Erratum E1:** observed order with non-constant ratios; the check value 1.991 is corrected to 1.981 | `phase3/specs/Research_08_Specification_Corrections_v2.2.md` |
| WS6 | Finite-bias experiment design | FB-1: definitions; falsifiable hypotheses H-A to H-D; minimal matrix (576 states, about 12 CPU-hours); stop/go stage (48 states, about 1 CPU-hour); inclusion, exclusion and failure rules. Diagnostic previews included | `phase3/reports/Finite_Bias_Discrepancy_Experiment_Design_2026-10-30.md` |
| WS7 | Novelty challenge | Claim-vs-prior-art matrix (K1–K6 against P1–P11, evidence levels E1/E2/E3 kept separate). Ward, Odeh & Cohen 1990 and Ward 1988 §6.4 overlap the most interesting observation (K3). Their full texts were not obtainable through this environment | `phase3/reports/Novelty_Challenge_Matrix_2026-10-30.md` (with the C-3 close addendum) |

**Also produced in C-3:**
- Final regression: `pytest semiconductor_sim phase3/tests` gives **146 passed, 4 skipped, exit 0** in both environments (pinned 3.12.3 / NumPy 2.3.5 / SciPy 1.17.0; current 3.13.18 / 2.5.3 / 1.18.1).
- Golden `stage3.json` SHA-256 `17ca3274…7828` is unchanged. The production diff vs `7073f20` is 0 lines (`results/phase3/C3/regression_final/`).
- New unit tests for the convergence classification. They found an analysis error, which was corrected and versioned (v1 kept, v2 cited).

## 2. Why `REFERENCE_PARTIAL`

**Supports reference quality (zero drain bias):**
- The discrete T1 identity holds within the pre-registered resolving power (E ≤ 10⁻⁶) on 71/72 points per path:
  - R-I(a) and R-I(b);
  - L = 1, 3, 10 µm;
  - four mesh levels;
  - subthreshold, moderate and strong inversion.
- There is no FAIL. The finite-δ path L3 is consistent with L2 on 18/18 points.
- The triple-point repair holds at all levels: continuity jump 0, Gauss closure at round-off, no artificial oxide constraints.
- G_h converges to the independent continuum value M0 at second order. The classification is asymptotic, and the difference is ≤ 0.8·U_disc.
- The R-I(b) boundary equations were read from the assembled matrix, not assumed.

**Prevents `REFERENCE_READY`:**
1. **One class I point** (R-I(a), L = 10 µm, level 3, V_G = 2.0 V). R-I(a) is not demonstrated at its finest level for the longest channel in strong inversion under the pre-registered L1 definition.
2. **The L1 path is weaker than specified.**
   - B5 is not met on 31/72 points.
   - B3 does not transfer to the T1 right-hand side.
   - Two decoupling diagnostics fail at noise level.
   - Fixing this needs a versioned, pre-registered L1 definition (equilibrated factorization) re-run on all points.
3. **Open item R-I(b)-1:** bottom-corner ownership, a deviation from the stated bottom condition on two nodes. The estimated impact is negligible, but it is not yet removed.
4. **R-I(c) is not implemented.** FB-1 requires three reservoir models, so that the reservoir dependence is not judged from two arbitrary choices.
5. **T1 is a zero-bias identity.** At finite V_D, reference correctness rests on:
   - conservation (δ_KCL ≤ 2.0e-14, δ_cont ≤ 2.9e-14 over all previews);
   - refinement previews on 2–3 levels;
   - no identity-level test.
   - This is inherent and is stated as a limitation, not a defect.

**Why not `REFERENCE_BLOCKED`:**
- No verdict is FAIL. Every limitation above is localized, quantified and has a defined remedy.
- None of them affects the 71 resolved points.

## 3. Why `NOVELTY_UNPROVEN`

- **What C-3 established is verification infrastructure:**
  - a deck-level discrete identity test (K1);
  - a TCAD triple-point treatment (K2);
  - a same-mesh comparison device (K4).
  - None of these is a T-ED-level scientific contribution by itself. **Passing T1 says nothing about D1 originality.**
- **The most interesting new observation (K3)** is the small-V_D structure of the reduction discrepancy: δ² law, reservoir-antisymmetric part ∝ L⁻², peak near threshold. It lies in the territory of Ward, Odeh & Cohen 1990 (abstract: interface/boundary-layer effects on the averaged channel conductivity) and Ward 1988 §6.4 (linear-regime I–V including end effects). Those full texts have **not** been read (E2/E3 only).
- **The C-3 lateral-refinement previews sharpen the D1 framing risk** (diagnostic, not results):
  - At L = 1 µm and V_D = 50 mV the two reservoir models give same-mesh discrepancies of **opposite sign** (−5.3e-4 and +5.5e-4).
  - The reservoir-antisymmetric part α is mesh-converged and about **20×** the reservoir-symmetric part σ.
  - σ is mostly lateral discretization error, with a small extrapolated limit: ≈ 2.8e-5 at L = 1 µm and ≈ 1.5e-7 at L = 10 µm (V_G = 1.25 V).
  - If this persists, "the GCA reduction error" measured against idealized 2-D references is mainly a property of the reference's reservoir model at short L.
- **Moving to `NOVELTY_POTENTIAL` requires both:**
  - (i) the Ward full texts read, and K3/K4 shown to go beyond them (coefficients, end-condition dependence, verified uncertainty);
  - (ii) the FB-1 stop/go stage showing a resolved, reservoir-independent σ (H-C not falsified), with R-I(c) included.

## 4. Claims ledger (what may and may not be said now)

| Statement | Status |
|---|---|
| "The R-I(a) and R-I(b) DEVSIM decks satisfy the discrete zero-bias identity within the stated resolving power on 71/72 tested configurations (per path)" | **Supported**, with the class I point named |
| "The 2-D conductance converges at second order to the independent M0 reference" | **Supported** (asymptotic, R04 §7.4, estimates not bounds) |
| "The L1 linear-solve error is ≤ 1e-8 for T1" | **Not supported** in 31/72 points |
| "R-I(b) imposes ψ = ψ_b0 on the whole bottom boundary" | **Not exactly.** Two corner nodes deviate (open item) |
| "S_res bounds the reservoir uncertainty" | **Not permitted.** It is a sensitivity between two specified models (Research_08 §7) |
| "The GCA reduction error at L = 1 µm is X" | **Not supported.** It is reservoir-model-dependent in the previews |
| Any novelty claim (K1–K6) | **Not supported** (`NOVELTY_UNPROVEN`) |

## 5. Readiness for a limited finite-bias discrepancy experiment

**Answer: not yet.** The project is technically able to run the 48-state FB-1 stop/go stage, but it is not scientifically ready. Running it before the preconditions below risks computing a quantity that prior art may already give analytically (H-A/H-B), or one whose meaning depends on an untested reservoir choice.

| # | Precondition (design §7 and this report) | Status | Who |
|---|---|---|---|
| P1 | Independent review of Packages A–D (R05 rule 6) | **Not done** | Project lead / external reviewer |
| P2 | R-I(c) implemented (reservoir strips, ℓ_res ∈ {0.1, 0.3, 1}·W_dep) and passing T1 with the Research_08 §5 checks | **Not done** (needs authorization; about 1 day of implementation plus the T1 matrix) | Next package |
| P3 | Full texts of Ward, Odeh & Cohen 1990 and Ward 1988 ch. 5–6 read, so that H-A/H-B are checked against existing derivations | **Not done.** Both are open-access Caltech PDFs, but the fetch proxy here refused them | **User:** download and attach the two PDFs (links in the novelty matrix) |
| P4 | T1 coverage of every deck and level used | Stop/go states (L = 1, 10 µm; V_G = 1.25 V; levels 0–3; models a, b): **all T1 PASS**. Full FB-1 additionally needs T1 for L = 2, 4 µm, and a resolution of the L = 10 µm, level 3, V_G = 2.0 V class I point (versioned L1) or its exclusion as class I | Next package |
| P5 | R-I(b)-1 bottom-corner fix | Recommended, not blocking (≤ 7.5e-12) | Next package |

## 6. Recommended next package (C-4) — for authorization; nothing started

In order:
1. **L1 v2** (pre-registered): equilibrated factorization (and optionally iterative refinement in extended precision) for the L1 path. Re-run T1 on all 72 points plus L = 2, 4 µm. Report v1 and v2 side by side. About 2 CPU-hours.
2. **R-I(b)-1 fix:** end the lateral contacts one node above the bottom. Re-verify from the exported matrix and re-run R-I(b) T1. Under 1 CPU-hour.
3. **R-I(c)** implementation, boundary verification and T1 (same checks as WS4).
4. **Literature gate:** read Ward 1990 and thesis §5–6 (P3). Update the novelty matrix to E1 for P1/P2. Decide whether H-A/H-B are already derived.
5. **Only then,** and after independent review: the FB-1 stop/go stage (48 states). Its outcome decides between the two branches in design §6:
   - (i) a verified reservoir-independent coefficient c_σ(V_G, L), a possible D1 core;
   - (ii) a methodological result: GCA-error estimates against idealized references are dominated by the reference's reservoir model. Modest novelty; likely not sufficient alone for T-ED.

## 7. Honest assessment of the publication path

- After C-3 the project has a **well-verified zero-bias reference capability** and a sharply posed, falsifiable question. It does **not** yet have a scientific contribution.
- The question's most promising answer (K3) is under direct threat from 1988–1990 asymptotic work that has not been read.
- If H-C fails, the remaining contribution is methodological and probably fits a verification-focused venue or a T-ED brief better than a full T-ED paper. That judgment should be revisited after P3 and the stop/go stage, not before.

## 8. Compliance with the C-3 restrictions

| Restriction | Status |
|---|---|
| `main` and golden baselines | Untouched (`main` = `c64c463`; golden SHA unchanged) |
| Production S1/S2/S3 physics | Unchanged (diff vs `7073f20`: 0 lines) |
| Failed and terminated runs | Kept: the class I point, the decoupling failures, `l1_diag/run1_normwise` (including the terminated job) |
| Reference parameters | Not modified to improve agreement |
| Uncertainty estimates | Presented as estimates, never as bounds |
| Novelty from T1 | Not claimed |
| Package E, large sweeps, manuscript | None launched |
| No force push; no overwrite of remote branches | Complied (no push at all) |
| Credentials or environments | None committed |
| Failed authentication method | Not retried |

## 9. Deliverables index

| # | Deliverable | Location |
|---|---|---|
| 1 | Preservation: bundle, source ZIP, SHA-256 manifests | Archive `mosfet-tcad-phase3-B-C1-C2-C3-2026-10-30.bundle`, `…C3…_source_<sha>.zip`, `SHA256SUMS.txt`, `MANIFEST_SHA256_files.txt`, `bundle_verification.txt` |
| 2 | R-I(a) formal T1 report | `phase3/reports/C3_T1_RIa_RIb_Report_2026-10-30.md` §1–§5 |
| 3 | R-I(b) boundary verification and T1 comparison | same, §6–§7 |
| 4 | C-1(b) final acceptance (+ Addendum A) | `phase3/reports/C1b_Final_Acceptance_Review_2026-10-30.md` |
| 5 | Specification corrections v2.2 (+ Erratum E1) | `phase3/specs/Research_08_Specification_Corrections_v2.2.md` |
| 6 | Finite-bias experiment design (FB-1) | `phase3/reports/Finite_Bias_Discrepancy_Experiment_Design_2026-10-30.md` |
| 7 | Novelty challenge matrix | `phase3/reports/Novelty_Challenge_Matrix_2026-10-30.md` |
| 8 | Raw data, analysis, figures, regression | `results/phase3/C3/` (`t1/`, `analysis_v2/`, `analysis/` v1, `t1_rib_bodycorners/`, `rib_inspect/`, `l1_diag/`, `fb/`, `regression_final/`) |
| 9 | Code | `phase3/c3/` (`t1.py`, `t1_analysis.py`, `rib_inspect.py`, `fb_preview.py`, `l1_refine_diag.py`, `run_t1_matrix.sh`); `phase3/tests/test_c3_analysis.py`; `phase3/c1/{deck,onecol}.py` (extended, defaults unchanged) |
| 10 | This decision report | `phase3/reports/C3_Research_Decision_Report_2026-10-30.md` |

The final HEAD SHA of `phase3/C3-2026-10-30` cannot be recorded inside its own commit. It is given in `bundle_verification.txt` and `SHA256SUMS.txt` of the archive.
