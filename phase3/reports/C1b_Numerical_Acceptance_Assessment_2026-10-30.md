# C-1(b) — Revised Numerical Acceptance Assessment

**Date:** 2026-10-30. **Script:** `phase3/c1/c1b_acceptance.py`. **Results:** `results/phase3/C2/c1b_acceptance/run3/`, reproduced bit-for-bit in `repro_fe44d68/c1b_acceptance`. Run 1 and run 2 are kept: their Skeel estimator failed (`RUNS.txt`).
**System:** the same legacy C-1 deck and states as C-1(b) (V_G = 1.25 V; V1 stepped 0 → 1 mV and 0 → 50 mV). 9 851 unknowns, nnz 77 217. Both stored C-1(b) systems are reproduced **bitwise** (ΔA = 0, Δrhs = 0, Δ(DEVSIM update) = 0).

## 1. Original observation (retained)

The external update (SciPy SuperLU on the exported A, b) differs from DEVSIM's own update (UMFPACK 5.1) by **1.4–1.9e-10** relative per unknown block. The circuit `drain_bias` block differs by 6.5e-16. R05's acceptance figure was ≤ 1e-10.

## 2. Measurements

| Quantity | dV1mV | dV50mV | Comment |
|---|---|---|---|
| Residual ‖Ax − b‖∞: UMFPACK / SuperLU | 4.1e-13 / 2.4e-15 | 5.2e-13 / 9.0e-14 | raw units (mixed V and cm⁻³ rows) |
| Normwise backward error, raw | 7e-33 / 4e-35 | 2e-34 / 3e-35 | meaningless: dominated by \|\|A\|\| × \|\|x\|\| of the cm⁻³ rows |
| **Normwise backward error, row+column equilibrated** | **1.7e-14 / 6.4e-16** | **3.3e-15 / 5.3e-16** | both solvers backward-stable (≤ ~150 u) |
| Componentwise ω (Oettli–Prager), all rows | 9.2e-2 / 1.0 | 1.4e-3 / 1.0 | not meaningful, see below |
| Componentwise ω, decoupled rows excluded | 9.2e-2 / 1.1e-4 | 1.4e-3 / 5.7e-4 | dominated by deep-bulk ψ rows whose exact update is ≈ 0 (\|x\| ~ 1e-16 V; denominator 1e-27 of its maximum) |
| Structure | 497 single-entry rows, 496 of them decoupled (b = 0) | same | contact Dirichlet rows |
| 1-norm condition estimate (Hager–Higham, a *lower-bound* estimate): raw / row-equilibrated / row+col-equilibrated | 5.5e35 / 1.2e28 / **9.8e10** | same | raw values only reflect unit disparity |
| **Skeel componentwise condition** cond(A, x) per block (estimate) | **3.1–3.6e7** (ψ, n, circuit V1.I); drain_bias 1.0 | same | relevant for the observed component errors |
| κ_Skeel · u (u = 1.1e-16) | ≈ 3.4e-9 | ≈ 3.4e-9 | first-order forward-error scale; **≈ 20× the observed 1.7–1.9e-10** |
| Reference solution by mixed-precision refinement (80-bit residuals) | relative corrections stagnate at 2.6–6.7e-11 | 1e-11 – 1.7e-10 | reference reliable only to ≈ 1e-10; residual cancellation in the long double accumulation limits it |
| Forward error vs that reference: UMFPACK / SuperLU / SuperLU + 3 double refinements | 1.15–1.34e-10 / 3.4–3.7e-11 / 3.3–3.8e-11 | 2.3–2.7e-10 / 6.4–7.5e-11 / 3.0–3.4e-11 | ranking UMFPACK < SuperLU is **within the reference's own uncertainty**; not claimed |
| Directional derivative (central differences of F = −rhs along 4 random, variable-scaled directions; h = 1e-1 … 1e-8) | best 3.2–8.9e-11 at h = 1e-6; V-shape (truncation ∝ h² down to 1e-6, round-off above) | identical | the exported matrix is the Jacobian of F at the current state |
| Column checks (C-1(b), 23 sampled unknowns incl. contact, gate and circuit) | ≤ 2.6e-10 | ≤ 2.6e-10 | from C-1 |
| Damped vs undamped Newton step, V1 = 50 mV | ψ: no component damped (\|Δψ\| ≤ 0.022 V < 0.0259 V); n: 991 components altered by `positive` (undamped n would be ≤ 0); max change 1.3e-3 of max\|Δn\| | — | damping is per component, not a line search |
| Damped vs undamped Newton step, V1 = 200 mV | ψ: 2 270 Si and 206 oxide components damped (`log_damp`), up to 56 % of max\|Δψ\|; n: 3 578 altered | — | — |
| Damped step reproduced from the undamped step with DEVSIM's formulas (`Equation.cc`: log_damp with the hard-coded 0.0259 V; positive: 0.001·n_old) | **≤ 9.1e-16** | — | damped and undamped updates may only be compared through these formulas |

## 3. Interpretation

1. The exported static matrix and RHS are **exactly** the Newton system at the current state:
   - bitwise identity with the factored system (C-1(b));
   - Jacobian of F = −rhs confirmed by column and directional derivatives to the finite-difference optimum.
2. Both linear solvers are **backward stable** on the equilibrated system. The cross-solver difference of 1.4–1.9e-10 is **ordinary forward error**, about 20× below the first-order scale κ_Skeel·u ≈ 3.4e-9.
3. **The 1e-10 equality figure cannot be derived from standard error analysis.** For two backward-stable solvers on this system, any bound is ≥ κ_Skeel·u ~ 3e-9; normwise bounds are looser still (κ_eq·u ~ 1e-5). Requiring ≤ 1e-10 would be a property of particular solver pairs, not of the matrix export.
4. The componentwise backward error is not a usable acceptance metric for this system. Its maximum sits in rows whose exact solution component is zero to round-off.
5. An undamped DEVSIM step equals A⁻¹·rhs. Damped steps differ by design (up to 56 % in ψ at 200 mV), and they are reproduced exactly from the undamped step with DEVSIM's own component-wise formulas. Any comparison must therefore use `variable_update='default'`, or apply those formulas, as done here.

## 4. Proposed revised acceptance criteria for C-1(b)

These replace "≤ 1e-10 vs DEVSIM's own update". They are a proposal; reviewer approval is required.

| ID | Criterion | Justification | Measured |
|---|---|---|---|
| B1 | `get_matrix_and_rhs` static matrix and RHS bitwise equal to the system factored by the Newton iteration at the same state (pattern and values); the export does not change the state | Definition of "is the Newton Jacobian" | 0.0 / 0.0 / 0.0 — **met** |
| B2 | Central-difference directional derivatives of F = −rhs agree with A·v to ≤ 1e-7 at the optimal step of a ladder (≥ 3 random scaled directions plus column checks at contact, interface, gate and circuit unknowns) | FD optimum for double precision is ~u^(2/3) ≈ 2e-11 times derivative-scale factors; 1e-7 leaves margin while still detecting a missing or mis-signed term (errors of O(1)) | ≤ 8.9e-11 (directions), ≤ 2.6e-10 (columns) — **met** |
| B3 | Normwise backward error of the equilibrated system ≤ 1e-13 for every linear solver used | ≈ 1000 u; detects a broken factorization | ≤ 1.7e-14 — **met** |
| B4 | Cross-solver forward difference per block ≤ 10 · κ_Skeel · u (≈ 3.4e-8 here), with κ_Skeel estimated at the state | First-order forward-error scale with a factor-10 margin for estimator bias (Hager–Higham estimates are lower bounds) | ≤ 1.9e-10 — **met** |
| B5 | For T1: the linear-solve part of E_lin (relative difference of G between two solvers, or between the L1 and Newton-step routes) ≤ 0.01·E_req = 1e-8 | Keeps the linear algebra two orders below the T1 resolving target | ≤ 1.1e-10 (C-1(a)), ≤ 1.9e-9 (corner study L2) — **met** |
| B6 | Damped and undamped steps are compared only after applying DEVSIM's update formulas | Avoids attributing damping to the Jacobian | ≤ 9.1e-16 — **met** |

The original discrepancy of 1.4–1.9e-10 is retained in the record; it is not "fixed".
