# C-1(b) — Final Numerical Acceptance Review of the DEVSIM Newton-System Interface

**Date:** 2026-10-30. **Scope:** final decision on the proposed criteria B1–B6 (`C1b_Numerical_Acceptance_Assessment_2026-10-30.md`).
**Evidence used (unchanged; no new run was needed):**
- `results/phase3/C1/repro_f4ca0bf/c1b/c1b.json`
- `results/phase3/C2/c1b_acceptance/run3/` (bitwise reproduced in `repro_fe44d68`)
- the C-3 T1 runs (`results/phase3/C3/t1/`), which exercise the same interface on 36 configurations

**Retained observation.** The external update (SciPy SuperLU on the exported A, b) differs from DEVSIM's own update (UMFPACK 5.1) by **1.4–1.9e-10** relative per unknown block, against the original R05 figure of ≤ 1e-10. This number stays in the record. It is neither re-labeled as a pass nor "fixed".

## 1. What C-1(b) is actually asking — six distinct questions

| # | Question | Kind of evidence | Can it be settled exactly? |
|---|---|---|---|
| Q1 | Algebraic identification: is the exported (A, b) the system the Newton iteration factors at the same state? | Bitwise comparison with the matrix and RHS captured inside the solver callback | Yes: exact equality or not |
| Q2 | Jacobian verification: is A = ∂F/∂x of the residual F = −b at the current state? | Finite-difference columns and directional derivatives | No: limited by FD truncation and round-off; an *agreement level* is reported |
| Q3 | Backward stability: does each linear solve return the exact solution of a nearby system? | Normwise backward error of the equilibrated system | Measurable per solve |
| Q4 | Cross-solver forward sensitivity: how much do two backward-stable solvers' solutions differ? | Direct difference | Measurable; its *expected* size is set by conditioning |
| Q5 | Conditioning: how sensitive is the solution to perturbations of this size? | Condition-number **estimates** (Hager–Higham 1-norm; Skeel) | Estimates only: Hager–Higham gives a lower bound on ‖A⁻¹‖₁; Skeel is computed through the same estimator |
| Q6 | Damped vs undamped updates: do DEVSIM's damped steps differ from A⁻¹b only through its documented update formulas? | Reconstruction of the damped step from the undamped one | Yes, up to round-off |

## 2. Decision per criterion

| ID | Proposed criterion (2026-10-30 assessment) | Measured | **Decision** | Required change to wording |
|---|---|---|---|---|
| **B1** | Static matrix and RHS bitwise equal to the factored Newton system at the same state; export does not change the state | ΔA = 0, Δb = 0, Δstate = 0 (2 states; reproduced at f4ca0bf and fe44d68) | **ACCEPT** | None. This is the *identification* result and the only criterion that answers "is the export the Newton system" |
| **B2** | Central-difference directional derivatives agree with A·v to ≤ 1e-7 at the optimal step of a ladder (≥ 3 random scaled directions plus column checks at contact, interface, gate and circuit unknowns) | ≤ 8.9e-11 (4 directions), ≤ 2.6e-10 (23 columns) | **ACCEPT, as verification of consistency** | Add: "agreement to finite-difference accuracy; this verifies A = ∂F/∂x to that accuracy, not exactly." The 1e-7 threshold is a detection level, chosen to expose missing or mis-signed terms (which produce O(1) errors). It is **not** a statement about the precision of A |
| **B3** | Normwise backward error of the row- and column-equilibrated system ≤ 1e-13 for every solver used | UMFPACK ≤ 1.7e-14; SuperLU ≤ 6.4e-16 | **ACCEPT** | Add: "computed after the documented row/column max-equilibration; unscaled values are dominated by unit disparity (V vs cm⁻³) and are not reported as evidence." Keep the note that componentwise backward error is degenerate on exactly-decoupled contact rows |
| **B4** | Cross-solver forward difference per block ≤ 10·κ_Skeel·u (≈ 3.4e-8 here) | 1.4–1.9e-10 | **ACCEPT AS A CONSISTENCY CHECK, NOT AS AN ERROR BOUND** | Rewrite as: "The cross-solver difference is reported. It is compared with the first-order perturbation scale κ_Skeel·u, where κ_Skeel is an **estimate** (Hager–Higham lower-bound estimator). The comparison checks consistency; it does not bound the true error. A difference far above that scale would indicate a solver or assembly fault." The factor 10 is a heuristic margin and is labeled so |
| **B5** | For T1: the linear-solve part of E_lin ≤ 0.01·E_req = 1e-8 | C-1(a) ≤ 1.1e-10; corner study ≤ 1.9e-9; C-3 T1 (all 36 configurations): see the T1 report | **ACCEPT, with the definition made explicit** | Specify that this term is \|G_L1 − G_NS\|/G plus the contact-vs-section difference of the same linear solution, as pre-registered in `phase3/c3/t1.py`. It is part of E (an estimate), not a bound |
| **B6** | Damped and undamped steps compared only after applying DEVSIM's update formulas | reconstruction ≤ 9.1e-16 | **ACCEPT** | Add the formula reference (DEVSIM `Equation.cc`: `log_damp` with its hard-coded 0.0259 V; `positive` limiting at 0.001·n_old). The T1 driver uses `variable_update='default'` for every Newton-step comparison |

**The original R05 figure.** "One Newton step computed externally must reproduce DEVSIM's own update to ≤ 1e-10 relative" is **withdrawn** as an acceptance criterion. It is recorded as **not met** (1.4–1.9e-10).

Why it is withdrawn:
- It compares two different factorizations.
- For two backward-stable solvers the expected difference scales with the conditioning of the system, estimated here at κ_Skeel·u ≈ 3e-9.
- A fixed 1e-10 equality figure has no derivation from that analysis. Whether it is met depends on the solver pair, not on whether the export is correct.

The question it was meant to answer, "is the export the Newton Jacobian?", is answered exactly by B1 and verified by B2.

## 3. Language rules adopted for all later reports

1. "Identical" is used **only** for bitwise equality (B1).
2. "Verified to finite-difference accuracy" is used for B2.
3. "Backward stable (estimated backward error X)" is used for B3.
4. Condition numbers are always called **estimates**. They are never written as "error ≤ κ·u", and never presented as guaranteed bounds.
5. Cross-solver differences are reported as **measured sensitivities**.

## 4. Status

**C-1(b): ACCEPTED** for its stated purpose: the DEVSIM matrix/RHS export can be used as the Newton linear system for L1, with B1–B6 as above. The ≤ 1e-10 R05 figure is withdrawn and recorded as not met. The superseding wording is in `Research_08_Specification_Corrections_v2.2` §4.
