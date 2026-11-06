# WP8 report: why the coupled small-signal solve fails, and grid convergence of the A-R references

**Date:** 2026-11-06 · **Author:** Claude · **Branch:** `benchmark/wp8-mechanism-gci-2026-11-06`
**Pre-registration:** `benchmark/specs/WP8_Preregistration_SmallSignal_Mechanism_and_Grid_Convergence_2026-11-06.md`
- predictions committed at `6337336`, frozen at `5d5a9af`;
- amendment 1 (exploratory analysis X) at `c34c3db`; amendment 2 (check X4) at `154521b`.

This is step A of the pre-submission plan agreed with the user on 2026-11-06.

## Summary

- **Q1, mechanism (registered).** DEVSIM's coupled small-signal output is wrong because of its data, not because of its solver.
  - The double-precision system of each of the 11 states was rebuilt bit for bit and solved exactly, by iterative refinement with double-double residuals.
  - The exact solution is still wrong in all 11 states: by 4.2e-4 (AR-6, the largest current) up to 1.6e3 relative (AR-1).
  - Rounding every Jacobian entry by one ulp moves G by 20% to 93% of the double-precision error.
  - **P-M1 fails** (1 of 10 gate states) and **P-M2 fails** (0 of 10). By the interpretation fixed in advance: *rounding the Jacobian entries alone moves G by at least 1% of the error. No double-precision solver can repair the output; it needs a higher-precision assembly and solve, or a cancellation-free formulation.*
- **Q1, where and why (exploratory, amendments 1 and 2).**
  - **Linear.** The response to the one-ulp perturbations is linear. The first-order adjoint prediction reproduces all 33 measured shifts to 1.1e-11 relative or better; the second-order terms are at most 1.1e-11 of the shifts.
  - **Where.** At least 99.9% of the first-order sensitivity lies in the electron-continuity rows of the drain-side n-type region, and at least 99.98% in the potential columns. There:
    - the potential response is uniform (δψ = 1 V/V);
    - the drain collects any injected current (λ = 1);
    - the potential derivatives of the edge fluxes reach 3.4e6 S/cm on the thin surface cells.
  - **Why.** The potential entries of each such row sum to zero to within one rounding, because the Scharfetter–Gummel flux depends only on potential differences. The large terms cancel, and the rounding of each entry leaves a current that the drain collects.
  - **The floor.** This leaves an absolute data floor u‖c‖₂ that does not depend on the gate bias: 3.2225e-9 S/cm on dev1 and on dev1-L2, 8.6136e-9 S/cm on dev2, the same to within 3.4e-6 relative in every state of a device.
  - **Estimating it.** The unrefined double-precision solutions give the same floor to 2.5e-12, so one adjoint solve with the same factorization estimates it without a reference.
  - **The actual errors.** The exact-solve errors were 0.015 to 1.7 times the floor.
- **Q2, grid convergence (registered). P-G1, P-G2 and P-G3 all hold.**
  - Convergence is monotone in 12 of 12 pairs, with observed order 1.98–2.22 and asymptotic ratio 0.998–1.000.
  - The level-2 references lie 1.2e-4 to 2.0e-3 from the Richardson-extrapolated values.
  - The registered tolerances contain those differences with a margin of 2.0–2.3, and contain GCI_fine in all 12 pairs.
  - The independent finite-volume value (P4) is closer to the extrapolated value than the finest network value in 6 of 6 states.
- **A correction to a WP2 diagnostic.** WP2's κ_G divided by |G_E|, the wrong double-precision output, instead of the true G; its range was 4.9e14 to 1.2e18. With the true G and exact solutions, the condition number of the output with respect to the entries is 4.8e14 to 8.7e20 (section 3). The manuscript now reports this value.

## 1. What was run

| Step | States | Instrument | Records |
|---|---|---|---|
| Smoke | dev1, V_G = 0.55 V, level 0 (not a study state) | `ssmech.py` | `results/smoke/` (not evaluated) |
| Q1 study, registered | 11 A-R and HO states, level 1, body B1 | `ssmech.py` via `run_wp8.py` | `results/production/SS_*.json` |
| Q2, registered | 6 A-R states × {G_DD, G_SD}, levels 0–2 (case files) | `gci_wp8.py` | `results/gci_wp8.json`, `results/gci_table.md` |
| Evaluation | — | `evaluate_wp8.py` | `results/evaluation_wp8.json`, `results/report_tables_wp8.md` |
| X, exploratory | the same 11 states | `xsens.py` via `run_xsens.py` | `results/exploratory/XS_*.json`, `summary_xsens.json` |
| X4, exploratory | AR-1 (dev1), AR-4 (dev2) | `xsens_rowsum.py` | `results/exploratory/X4_*.json` |

- Every run is a fresh DEVSIM 2.11.2 process: Python 3.12.3, NumPy 2.3.5, SciPy 1.17.0, reference BLAS/LAPACK, one thread.
- Every run has a ledger line with the record's SHA-256 (`results/ledger.jsonl`).
- All study runs exited 0. Wall time per state:
  - dev1 states: 22–27 s;
  - dev2 states: 149–153 s, almost all of it building the equilibrium.
- The exact solve (M1) took at most 1.6 s, and the three perturbed solves (M2) at most 4.8 s.

Devices:
- **dev1:** L = 1 µm, N_A = 1e18 cm⁻³, t_ox = 2 nm. States AR-1, AR-2, AR-3, HO-1, HO-2.
- **dev2:** L = 2 µm, N_A = 1e16 cm⁻³, t_ox = 20 nm. States AR-4, AR-5, AR-6, HO-3, HO-4.
- **dev1-L2:** dev1 with L = 2 µm. State HO-5.

## 2. Q1, registered results

**P-M0 (instrument): holds.** All 11 states reproduce WP2's G_E, G_M and G_net bit for bit (`same_*` true), so none is excluded.

**Exact solve (M1).**
- Converged in all 11 states, in 3–4 iterations, to a final double-double componentwise backward error of 6e-32 to 2.4e-26.
- The starting solution, equilibrated SuperLU, has a componentwise backward error of up to 0.71 (AR-1) in double-double, although its normwise backward error is 1.6e-17.
- M1b (UMFPACK, factorized once) agrees with M1 to at most 6.7e-15 relative in every state.

| State | Device | G_net (S/cm) | e_E | e_M | e_dbl | e_IR | e_IR/e_dbl | S | S/e_dbl | G_IR − G_net (S/cm) | it | ω_dd | M1b vs M1 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| AR-1 | dev1 | 3.5326e-12 | 6.28e+03 | 1.74e+03 | 1.74e+03 | 1.56e+03 | 0.896 | 1.62e+03 | 0.929 | −5.50e-09 | 4 | 2.8e-28 | 0 |
| AR-2 | dev1 | 1.6192e-10 | 137 | 41.7 | 41.7 | 27.2 | 0.652 | 35.2 | 0.846 | −4.40e-09 | 4 | 3.4e-30 | 0 |
| AR-3 | dev1 | 2.1147e-08 | 1.06 | 0.313 | 0.313 | 0.236 | 0.755 | 0.27 | 0.863 | −5.00e-09 | 4 | 1.7e-30 | 0 |
| AR-4 | dev2 | 2.5505e-09 | 22.6 | 25.1 | 22.6 | 1.31 | 0.0579 | 6.64 | 0.294 | +3.34e-09 | 4 | 6.4e-32 | −3.3e-16 |
| AR-5 | dev2 | 1.2035e-07 | 0.461 | 0.553 | 0.461 | 0.0624 | 0.135 | 0.141 | 0.305 | +7.51e-09 | 4 | 6.6e-32 | 0 |
| AR-6 | dev2 | 1.5976e-05 | 0.00358 | 0.00458 | 0.00358 | 0.000422 | 0.118 | 0.00106 | 0.295 | +6.74e-09 | 3 | 6.7e-32 | 0 |
| HO-1 | dev1 | 6.1776e-11 | 363 | 112 | 112 | 74 | 0.661 | 92.4 | 0.825 | −4.57e-09 | 4 | 1.9e-29 | 0 |
| HO-2 | dev1 | 1.1211e-09 | 19.5 | 6.18 | 6.18 | 4.15 | 0.671 | 5.09 | 0.823 | −4.65e-09 | 4 | 2.9e-31 | +2.2e-16 |
| HO-3 | dev2 | 4.5692e-08 | 1.4 | 1.59 | 1.4 | 0.0204 | 0.0145 | 0.37 | 0.264 | −9.32e-10 | 4 | 6.1e-32 | 0 |
| HO-4 | dev2 | 8.4012e-07 | 0.0753 | 0.0809 | 0.0753 | 0.00523 | 0.0694 | 0.0201 | 0.267 | +4.39e-09 | 4 | 9.2e-32 | 0 |
| HO-5 | dev1-L2 | 1.7171e-12 | 7.0e+03 | 7.79e+03 | 7.0e+03 | 28.4 | 0.00406 | 1.41e+03 | 0.201 | +4.88e-11 | 4 | 2.4e-26 | −6.7e-15 |

Column definitions:
- e_E, e_M: |G/G_net − 1| of the two double-precision solves; e_dbl = min(e_E, e_M).
- e_IR: the same for the exact solve.
- S: the largest one-ulp shift |G_k − G_IR| / |G_net| over three seeds.

**Gate.** The gate states (e_dbl > 1e-2) are the 10 states other than AR-6, as anticipated in the pre-registration.

**P-M1 (the solve): fails.**
- e_IR ≤ 1e-2 · e_dbl holds in 1 of 10 gate states (HO-5).
- The exact solve leaves:
  - 65–90% of the error on dev1;
  - 1.5–13.5% on dev2;
  - 0.4% on HO-5.

**P-M2 (the data): fails.** S ≤ 1e-2 · e_dbl holds in 0 of 10 gate states; S/e_dbl = 0.20–0.93.

**Registered interpretation (row "either / fails").** Rounding the Jacobian entries alone moves G by at least 1% of the error. No double-precision solver can repair the output; it needs a higher-precision assembly and solve, or a cancellation-free formulation.

**Descriptive.**
- No state is resolved to 1e-6 by the exact solve. Even AR-6, at 1.6e-5 S/cm, keeps e_IR = 4.2e-4.
- Within a device, the absolute error of the exact solve does not depend on the state. On dev1 it is −4.4e-9 to −5.5e-9 S/cm; on dev2 it is −0.9e-9 to +7.5e-9 S/cm. So the relative error scales as 1/G.
- The M2 shift of each seed is the same to three digits in every state of a device:
  - dev1: −2.23e-9, −5.71e-9, −5.07e-9 S/cm;
  - dev2: −2.59e-9, −1.69e-8, −1.46e-8 S/cm.
- On dev2 and dev1-L2 the factorization error is most of the double-precision error: the exact solve removes 86–99.6% of it. What remains, the data floor, is still 20–30% of e_dbl (S), so a better solver alone does not make the output usable.

## 3. Exploratory analysis X (amendments 1 and 2; no predictions)

Specified after the first records were seen and before it was run. Reported separately from the registered results.

Records: `results/exploratory/XS_*.json`, summarized in `summary_xsens.json` and `summary_tables_xsens.md`. The tables below are copied from `summary_tables_xsens.md`.

**X1, linearity.** For every seed of every state the first-order prediction reproduces the measured shift to at most 1.1e-11 relative. With the second-order term added it matches to 5.3e-13; the second-order term is at most 1.1e-11 of the shift. Re-solving the perturbed systems reproduced the M2 records bit for bit (33 of 33). The response is linear, so first-order sensitivity analysis is exact for this question.

**X2, size.**
- **Condition number.** κ_A = Σ|c_ij| / |G_net| ranges from 4.8e14 (AR-6) to 8.7e20 (HO-5).
- **The floor.** The absolute floor u‖c‖₂, with u = 2⁻⁵³, is 3.2225e-9 S/cm on dev1 and dev1-L2 and 8.6136e-9 S/cm on dev2. It is the same to within 3.4e-6 relative in every state of a device, so it does not depend on the gate bias.
- **Relative error.** The relative floor therefore scales as 1/G.
- **Unrefined solutions.** Computed from the unrefined double-precision x and λ (the WP2 computation), the floor agrees with the refined one to 2.5e-12. The adjoint and the solution are accurate where they matter (λ = 1, δψ = 1 in the drain); their errors sit in the small components that carry G.
- **WP2's κ_G.** It equals the unrefined κ_A times |G_net/G_E|, so it was misnormalized, not miscomputed.

**X3, location.**
- In every state at least 99.90% of Σ|c_ij| lies in the electron-continuity rows of the drain-side n-type region, and at least 99.98% in the potential columns.
- The circuit rows carry 0.1% on dev1; the source side, the p-type body, the oxide and the contact rows carry none to the printed precision.
- The largest rows are at the surface next to the drain contact.

#### X2: first-order sensitivity to the Jacobian entries (refined x and lambda)

| State | G_net (S/cm) | kappa_A | u (sum c^2)^1/2 (S/cm) | u sum|c| (S/cm) | same, unrefined | WP2 kappa_G | |G_IR - G_net| (S/cm) | x vs lambda |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| AR-1 | 3.533e-12 | 4.22e+20 | 3.22e-09 | 1.65e-07 | 3.22e-09 | 6.7e+16 | 5.50e-09 | +0.0e+00 |
| AR-2 | 1.619e-10 | 9.20e+18 | 3.22e-09 | 1.65e-07 | 3.22e-09 | 6.8e+16 | 4.40e-09 | +0.0e+00 |
| AR-3 | 2.115e-08 | 7.05e+16 | 3.22e-09 | 1.65e-07 | 3.22e-09 | 1.2e+18 | 5.00e-09 | +0.0e+00 |
| AR-4 | 2.550e-09 | 3.03e+18 | 8.61e-09 | 8.58e-07 | 8.61e-09 | 1.4e+17 | 3.34e-09 | -3.3e-16 |
| AR-5 | 1.204e-07 | 6.42e+16 | 8.61e-09 | 8.58e-07 | 8.61e-09 | 1.2e+17 | 7.51e-09 | +0.0e+00 |
| AR-6 | 1.598e-05 | 4.84e+14 | 8.61e-09 | 8.58e-07 | 8.61e-09 | 4.9e+14 | 6.74e-09 | +0.0e+00 |
| HO-1 | 6.178e-11 | 2.41e+19 | 3.22e-09 | 1.65e-07 | 3.22e-09 | 6.7e+16 | 4.57e-09 | +0.0e+00 |
| HO-2 | 1.121e-09 | 1.33e+18 | 3.22e-09 | 1.65e-07 | 3.22e-09 | 7.2e+16 | 4.65e-09 | +0.0e+00 |
| HO-3 | 4.569e-08 | 1.69e+17 | 8.61e-09 | 8.58e-07 | 8.61e-09 | 4.2e+17 | 9.32e-10 | +0.0e+00 |
| HO-4 | 8.401e-07 | 9.20e+15 | 8.61e-09 | 8.58e-07 | 8.61e-09 | 9.9e+15 | 4.39e-09 | +0.0e+00 |
| HO-5 | 1.717e-12 | 8.68e+20 | 3.22e-09 | 1.65e-07 | 3.22e-09 | 1.2e+17 | 4.88e-11 | +6.7e-15 |

#### X1: response to the M2 perturbations (S/cm)

| State | seed | measured G_k - G_IR | first order dG1 | second order dG2 | dG1/measured - 1 | (dG1+dG2)/measured - 1 | rho rms | mean (ulp) | = M2 record |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| AR-1 | 0 | -2.2384e-09 | -2.2384e-09 | -1.84e-21 | -8.22e-13 | +0.00e+00 | 1.568e-16 | +0.0087 | True |
| AR-1 | 1 | -5.7086e-09 | -5.7086e-09 | -1.53e-21 | -2.68e-13 | +0.00e+00 | 1.568e-16 | +0.0101 | True |
| AR-1 | 2 | -5.0693e-09 | -5.0693e-09 | -1.38e-21 | -2.72e-13 | +2.22e-16 | 1.568e-16 | +0.0094 | True |
| AR-2 | 0 | -2.2350e-09 | -2.2350e-09 | -1.84e-21 | -8.23e-13 | +0.00e+00 | 1.568e-16 | +0.0087 | True |
| AR-2 | 1 | -5.7070e-09 | -5.7070e-09 | -1.53e-21 | -2.68e-13 | +0.00e+00 | 1.568e-16 | +0.0101 | True |
| AR-2 | 2 | -5.0687e-09 | -5.0687e-09 | -1.38e-21 | -2.72e-13 | +2.22e-16 | 1.568e-16 | +0.0096 | True |
| AR-3 | 0 | -2.2345e-09 | -2.2345e-09 | -1.84e-21 | -8.23e-13 | -5.55e-16 | 1.568e-16 | +0.0087 | True |
| AR-3 | 1 | -5.7103e-09 | -5.7103e-09 | -1.53e-21 | -2.68e-13 | -1.11e-16 | 1.568e-16 | +0.0102 | True |
| AR-3 | 2 | -5.0654e-09 | -5.0654e-09 | -1.38e-21 | -2.73e-13 | -1.11e-16 | 1.568e-16 | +0.0094 | True |
| AR-4 | 0 | -2.6020e-09 | -2.6020e-09 | -6.84e-21 | -2.63e-12 | -3.33e-16 | 1.565e-16 | +0.0089 | True |
| AR-4 | 1 | -1.6926e-08 | -1.6926e-08 | -1.10e-20 | -6.48e-13 | +0.00e+00 | 1.565e-16 | +0.0088 | True |
| AR-4 | 2 | -1.4648e-08 | -1.4648e-08 | -8.92e-21 | -6.09e-13 | +2.22e-16 | 1.565e-16 | +0.0091 | True |
| AR-5 | 0 | -2.5926e-09 | -2.5926e-09 | -6.85e-21 | -2.65e-12 | -6.55e-15 | 1.565e-16 | +0.0090 | True |
| AR-5 | 1 | -1.6921e-08 | -1.6921e-08 | -1.10e-20 | -6.49e-13 | -7.77e-16 | 1.565e-16 | +0.0091 | True |
| AR-5 | 2 | -1.4644e-08 | -1.4644e-08 | -8.92e-21 | -6.10e-13 | -4.44e-16 | 1.564e-16 | +0.0093 | True |
| AR-6 | 0 | -2.5881e-09 | -2.5881e-09 | -6.85e-21 | -2.12e-12 | +5.29e-13 | 1.565e-16 | +0.0093 | True |
| AR-6 | 1 | -1.6915e-08 | -1.6915e-08 | -1.10e-20 | -6.63e-13 | -1.49e-14 | 1.565e-16 | +0.0093 | True |
| AR-6 | 2 | -1.4643e-08 | -1.4643e-08 | -8.92e-21 | -5.56e-13 | +5.37e-14 | 1.565e-16 | +0.0095 | True |
| HO-1 | 0 | -2.2339e-09 | -2.2339e-09 | -1.84e-21 | -8.23e-13 | +2.22e-16 | 1.568e-16 | +0.0085 | True |
| HO-1 | 1 | -5.7075e-09 | -5.7075e-09 | -1.53e-21 | -2.68e-13 | +2.22e-16 | 1.567e-16 | +0.0102 | True |
| HO-1 | 2 | -5.0702e-09 | -5.0702e-09 | -1.38e-21 | -2.72e-13 | +2.22e-16 | 1.568e-16 | +0.0095 | True |
| HO-2 | 0 | -2.2349e-09 | -2.2349e-09 | -1.84e-21 | -8.23e-13 | -3.33e-16 | 1.568e-16 | +0.0088 | True |
| HO-2 | 1 | -5.7074e-09 | -5.7074e-09 | -1.53e-21 | -2.68e-13 | +0.00e+00 | 1.568e-16 | +0.0102 | True |
| HO-2 | 2 | -5.0671e-09 | -5.0671e-09 | -1.38e-21 | -2.72e-13 | -1.11e-16 | 1.568e-16 | +0.0096 | True |
| HO-3 | 0 | -2.5923e-09 | -2.5923e-09 | -6.85e-21 | -2.64e-12 | -4.44e-16 | 1.565e-16 | +0.0090 | True |
| HO-3 | 1 | -1.6922e-08 | -1.6922e-08 | -1.10e-20 | -6.48e-13 | +0.00e+00 | 1.565e-16 | +0.0090 | True |
| HO-3 | 2 | -1.4646e-08 | -1.4646e-08 | -8.92e-21 | -6.09e-13 | +2.22e-16 | 1.564e-16 | +0.0093 | True |
| HO-4 | 0 | -2.5914e-09 | -2.5914e-09 | -6.85e-21 | -2.64e-12 | +4.44e-16 | 1.565e-16 | +0.0092 | True |
| HO-4 | 1 | -1.6920e-08 | -1.6920e-08 | -1.10e-20 | -6.47e-13 | +6.66e-16 | 1.565e-16 | +0.0092 | True |
| HO-4 | 2 | -1.4646e-08 | -1.4646e-08 | -8.92e-21 | -6.04e-13 | +4.66e-15 | 1.565e-16 | +0.0093 | True |
| HO-5 | 0 | +2.4171e-09 | +2.4171e-09 | -1.64e-21 | +6.78e-13 | -5.55e-16 | 1.601e-16 | +0.0121 | True |
| HO-5 | 1 | -1.4478e-09 | -1.4478e-09 | -5.62e-22 | -3.88e-13 | -1.11e-16 | 1.600e-16 | +0.0139 | True |
| HO-5 | 2 | -1.4799e-10 | -1.4799e-10 | -1.61e-21 | -1.09e-11 | +1.55e-15 | 1.600e-16 | +0.0129 | True |

#### X3: share of sum|c_ij| by the class of the row node

| State | circuit | contact_body | contact_drain | contact_source | n_drain_side | n_source_side | oxide | p_type |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| AR-1 | 0.001 | 0.000 | 0.000 | 0.000 | 0.999 | 0.000 | 0.000 | 0.000 |
| AR-2 | 0.001 | 0.000 | 0.000 | 0.000 | 0.999 | 0.000 | 0.000 | 0.000 |
| AR-3 | 0.001 | 0.000 | 0.000 | 0.000 | 0.999 | 0.000 | 0.000 | 0.000 |
| AR-4 | 0.000 | 0.000 | 0.000 | 0.000 | 1.000 | 0.000 | 0.000 | 0.000 |
| AR-5 | 0.000 | 0.000 | 0.000 | 0.000 | 1.000 | 0.000 | 0.000 | 0.000 |
| AR-6 | 0.000 | 0.000 | 0.000 | 0.000 | 1.000 | 0.000 | 0.000 | 0.000 |
| HO-1 | 0.001 | 0.000 | 0.000 | 0.000 | 0.999 | 0.000 | 0.000 | 0.000 |
| HO-2 | 0.001 | 0.000 | 0.000 | 0.000 | 0.999 | 0.000 | 0.000 | 0.000 |
| HO-3 | 0.000 | 0.000 | 0.000 | 0.000 | 1.000 | 0.000 | 0.000 | 0.000 |
| HO-4 | 0.000 | 0.000 | 0.000 | 0.000 | 1.000 | 0.000 | 0.000 | 0.000 |
| HO-5 | 0.001 | 0.000 | 0.000 | 0.000 | 0.999 | 0.000 | 0.000 | 0.000 |

#### X3: share by equation block (rows / columns)

| State | si_potential rows | si_electrons rows | ox_potential rows | circuit rows | si_potential cols | si_electrons cols | ox_potential cols | circuit cols |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| AR-1 | 0.000 | 0.999 | 0.000 | 0.001 | 1.000 | 0.000 | 0.000 | 0.000 |
| AR-2 | 0.000 | 0.999 | 0.000 | 0.001 | 1.000 | 0.000 | 0.000 | 0.000 |
| AR-3 | 0.000 | 0.999 | 0.000 | 0.001 | 1.000 | 0.000 | 0.000 | 0.000 |
| AR-4 | 0.000 | 1.000 | 0.000 | 0.000 | 1.000 | 0.000 | 0.000 | 0.000 |
| AR-5 | 0.000 | 1.000 | 0.000 | 0.000 | 1.000 | 0.000 | 0.000 | 0.000 |
| AR-6 | 0.000 | 1.000 | 0.000 | 0.000 | 1.000 | 0.000 | 0.000 | 0.000 |
| HO-1 | 0.000 | 0.999 | 0.000 | 0.001 | 1.000 | 0.000 | 0.000 | 0.000 |
| HO-2 | 0.000 | 0.999 | 0.000 | 0.001 | 1.000 | 0.000 | 0.000 | 0.000 |
| HO-3 | 0.000 | 1.000 | 0.000 | 0.000 | 1.000 | 0.000 | 0.000 | 0.000 |
| HO-4 | 0.000 | 1.000 | 0.000 | 0.000 | 1.000 | 0.000 | 0.000 | 0.000 |
| HO-5 | 0.000 | 0.999 | 0.000 | 0.001 | 1.000 | 0.000 | 0.000 | 0.000 |


**X4, the cancellation** (`X4_AR-1_B1_lev1.json`, `X4_AR-4_B1_lev1.json`). Over the electron-continuity rows of the drain-side n-type region, contact nodes excluded (4935 rows on AR-1, 19719 on AR-4):
- **Row sums.** The stored potential-column entries of each row sum (exactly, `math.fsum`) to at most 1.2e-16 of their absolute sum. The median is 2.4e-17 (AR-1) and 2.3e-17 (AR-4). The full row behaves the same, because the electron-column entries are negligible in these units.
- **Adjoint.** λ = 1 to within 8e-11 in the whole region.
- **Potential response.** δψ = 1 to within 1e-11 at the median. It falls to 0.41 (AR-1) and 0.73 (AR-4) only in the depleted part of the n-side junction, where n is as low as 2e11 and 1e15 cm⁻³ and the sensitivity is negligible.
- **Largest rows.** The ten largest-s_i rows lie at the surface (first two node rows below the oxide, y ≤ 0.7 nm) next to the drain contact. There n = 1e20 cm⁻³, the largest entries are 2.5e6–3.4e6 S/cm, and λ = δψ = 1.

**Interpretation (exploratory).**
- **Shift invariance and its rounding.** The Scharfetter–Gummel flux of an edge depends on the potentials only through their difference (`vdiff = (Potential@n0 − Potential@n1)/V_t` in the deck). The potential columns of an electron-continuity row therefore sum to zero in exact arithmetic. The recombination term's potential derivative is negligible where holes are few.
- **Why the drain dominates.** In the neutral drain the whole region follows the drain voltage, δψ = 1, so each row applies entries up to 3e6 S/cm to a constant vector, and the result is their sum. Rounded entries do not sum to zero: the diagonal is the rounded sum of several edge derivatives. The leftover, of order u × 3e6 S/cm per row, is a spurious current that the drain collects completely (λ = 1).
- **The floor.** Summed over the region it gives the floor u‖c‖₂. The floor is fixed by the drain region, the mesh and the doping, not by the channel, which explains the state-independent absolute numbers.
- **How the actual error compares.** The actual error of DEVSIM's assembly plays the role of one draw of such a perturbation. The exact-solve errors were 0.015 to 1.7 times the floor (X2 table). HO-5 shows the lower end: it has the dev1 floor but an exact-solve error 66 times smaller. So the floor is the uncertainty of G given the double-precision data, not a bound on the error of one particular assembly.
- **Contrast with the nonlinear residual (untested).** The residual evaluates each flux from the potential difference, which a uniform shift leaves unchanged. That may be why bias-step differences of the circuit current, which never form this Jacobian product, were accurate to 1.2e-9 (WP2). This is a hypothesis; WP8 did not test it.

**The sign of the M2 shifts.**
- 32 of the 33 shifts are negative. There are 9 independent draws (3 seeds × 3 devices), and 8 of them are negative; only HO-5's first seed is positive.
- The response is linear (X1), so a common sign would come from the realized perturbations, not from second-order effects.
- The realized relative perturbation is slightly biased: mean −1.2e-17 to −1.3e-17, about 8% of its RMS 1.57e-16, because A_ij(1 + 2⁻⁵² u) is rounded to double. A mantissa-dependent bias of this kind can give the seeds of one matrix a common first-order component.
- Not quantified further; it concerns the M2 instrument, not the device.

## 4. Q2, grid convergence (registered)

`gci_wp8.py` applies the Celik et al. (2008) procedure:
- h ∝ N^−1/2, with r21 and r32 between 1.410 and 1.415;
- safety factor 1.25;
- input: the P1 values at levels 0, 1 and 2 from the released case files.

| State | quantity | p | f(L2)/f_ext − 1 | GCI_fine | asympt. ratio | tolerance | within | P4/f_ext − 1 |
|---|---|---:|---:|---:|---:|---:|---|---:|
| AR-1 | G_DD | 2.21 | +1.26e-03 | 1.57e-03 | 0.999 | 2.85e-03 | yes | −1.59e-04 |
| AR-1 | G_SD | 2.22 | +1.25e-03 | 1.56e-03 | 0.999 | 2.85e-03 | yes | −1.73e-04 |
| AR-2 | G_DD | 2.22 | +1.60e-03 | 1.99e-03 | 0.998 | 3.65e-03 | yes | −2.30e-04 |
| AR-2 | G_SD | 2.22 | +1.60e-03 | 1.99e-03 | 0.998 | 3.65e-03 | yes | −2.31e-04 |
| AR-3 | G_DD | 2.21 | +2.04e-03 | 2.54e-03 | 0.998 | 4.66e-03 | yes | −2.94e-04 |
| AR-3 | G_SD | 2.21 | +2.04e-03 | 2.54e-03 | 0.998 | 4.66e-03 | yes | −2.94e-04 |
| AR-4 | G_DD | 2.02 | +1.15e-04 | 1.44e-04 | 1.000 | 2.60e-04 | yes | −1.47e-05 |
| AR-4 | G_SD | 2.02 | +1.15e-04 | 1.44e-04 | 1.000 | 2.61e-04 | yes | −1.49e-05 |
| AR-5 | G_DD | 2.00 | +1.19e-04 | 1.49e-04 | 1.000 | 2.50e-04 | yes | −5.98e-06 |
| AR-5 | G_SD | 2.00 | +1.19e-04 | 1.49e-04 | 1.000 | 2.50e-04 | yes | −5.98e-06 |
| AR-6 | G_DD | 1.98 | +1.19e-04 | 1.49e-04 | 1.000 | 2.36e-04 | yes | +1.93e-05 |
| AR-6 | G_SD | 1.98 | +1.19e-04 | 1.49e-04 | 1.000 | 2.36e-04 | yes | +1.93e-05 |

- **P-G1 holds:** monotone convergence in 12 of 12 pairs (required ≥ 10).
- **P-G2 holds:** 1 ≤ p ≤ 3 in every monotone pair. The range is 1.98–2.22, consistent with the second-order scheme.
- **P-G3 holds:** |f(L2)/f_ext − 1| ≤ tolerance in 12 of 12 pairs.

Descriptive:
- The tolerance is 1.98–2.28 times the distance to f_ext and contains GCI_fine in 12 of 12 pairs.
- The asymptotic ratio is 0.998–1.000, so the three meshes are in the asymptotic range.
- In 6 of 6 states, P4 (on its own refined mesh) lies closer to f_ext than f(L2) does: 6e-6 to 2.9e-4 from f_ext. This independent implementation therefore supports f_ext as the continuum limit.
- The released tolerances are unchanged.

## 5. Consequences

- **Manuscript, section 3.3.** A grid-convergence paragraph and a table (`tab_gci`) were added.
- **Manuscript, section 5.3.**
  - The condition-number sentence now uses the correctly normalized κ_J.
  - New paragraph "The data, not the solver": the registered results, with the interpretation as fixed.
  - New paragraph "Where the sensitivity lies (exploratory)", which states what it is, and the floor lines in Fig. 3(b).
- **Manuscript, other sections.** The abstract, introduction, recommendations and conclusions follow. Refinement does not rescue the output, and the floor can be estimated from one adjoint solve.
- **Benchmark README, atlas entry.** κ_G ≈ 1e17 is replaced by the WP8 findings.
- **Asset script.** The figure PDFs are now written without a creation date. Two runs of `make_assets.py` give byte-identical figures and generated files.
- **Release.** The release notes and the snapshot import test cover WP8; the snapshot still has to be rebuilt once the WP7 author items are in.
- **Unchanged:** the benchmark cases, tolerances and scorer, and every WP1–WP7 record.

## 6. Deviations and disclosures

1. **Smoke run (section 7 of the pre-registration, before any study run).**
   - The registered M1b inner solver, plain SuperLU, diverged, and was replaced by UMFPACK factorized once.
   - A convergence requirement ω_dd ≤ 1e-20 was added.
2. **Grid instrument.** `gci_wp8.tolerance()` first read a non-existent key and returned None. It was corrected to read the registered tolerance `uncertainty.tolerance.<q>_rel` before the evaluation was run. The observed orders and extrapolated values had been computed at that point; no threshold or tolerance changed.
3. **Exploratory analyses.** Amendments 1 and 2 were added after the first study records were seen. They are specified in the pre-registration before being run, carry no predictions, and are reported as exploratory.
4. **Scoping disclosure (pre-registration section 2).** One grid sequence (AR-1 G_DD) was seen before the predictions were written.
5. **WP2's κ_G.** It was normalized by the double-precision output G_E. This report gives the correctly normalized value. The WP2 records are unchanged (append-only).
6. **Author.** The same model (Claude) designed, ran and evaluated the study; no independent analyst checked it.

## 7. Files

- **Instruments:**
  - `benchmark/wp8/ddlinalg.py`, `ssmech.py`, `run_wp8.py`, `gci_wp8.py`, `evaluate_wp8.py`;
  - exploratory: `xsens.py`, `run_xsens.py`, `xsens_summary.py`, `xsens_rowsum.py`.
- **Tests:** `benchmark/tests/test_wp8.py`, 6 tests: exact transforms, the double-double residual, refinement, vector and adjoint refinement against rational arithmetic, the three-grid procedure, and non-convergence detection.
- **Records:**
  - `benchmark/wp8/results/{smoke,production,exploratory}/`, `ledger.jsonl`;
  - `evaluation_wp8.json`, `gci_wp8.json`, `report_tables_wp8.md`, `gci_table.md`;
  - `exploratory/summary_xsens.json`, `exploratory/summary_tables_xsens.md`.
- **Manuscript:** `paper/main.tex`, `paper/main.pdf`, regenerated by `paper/scripts/make_assets.py`.
