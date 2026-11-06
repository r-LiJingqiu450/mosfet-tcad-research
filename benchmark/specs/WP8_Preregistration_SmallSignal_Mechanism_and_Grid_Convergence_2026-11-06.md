# WP8 pre-registration: why the coupled small-signal solve fails, and grid convergence of the A-R references

**Date:** 2026-11-06 (MST) · **Author:** Claude · **Branch:** `benchmark/wp8-mechanism-gci-2026-11-06`, from the WP7 head `75f989d`.

**Status: FROZEN** (2026-11-06, after the smoke run and before any study run).
- The predictions of sections 3 and 4 were committed at `6337336`, before any WP8 computation on a study state and before the smoke run. They are unchanged.
- The smoke run required two procedure changes, listed in section 7. Later changes are dated amendments appended at the end.

This is step A of the pre-submission plan agreed with the user on 2026-11-06.

## 1. Questions

- **Q1. Mechanism.**
  - DEVSIM's coupled small-signal output misses the network conductance on the same mesh by up to 7.8e3 relative (WP2 P-SS), although its normwise backward error is at most 3e-17. DEVSIM's AC analysis returns the same solve (WP4 B2).
  - Is the error caused by the rounding of the double-precision factorization (the solve)?
  - Or is the double-precision system itself, with its rounded Jacobian entries, unable to determine G?
- **Q2. Grid convergence.** Do the A-R reference conductances converge under the published refinement with an observed order consistent with the discretization? Do the registered A-R tolerances contain the Richardson-extrapolated values?

## 2. What is known before WP8

- **WP2 P3 (level 1, body B1, 11 states).** The small-signal system was solved twice in double precision:
  - G_E: SuperLU after row and column equilibration;
  - G_M: UMFPACK.

  Relative errors against the network G_net (`benchmark/wp2/results/production/P1_*_B1_lev1.json`):

  | State | G_E | G_M | State | G_E | G_M |
  |---|---:|---:|---|---:|---:|
  | AR-1 | −6.28e3 | +1.74e3 | HO-1 | −3.63e2 | +1.12e2 |
  | AR-2 | −1.37e2 | +4.17e1 | HO-2 | −1.95e1 | +6.18 |
  | AR-3 | −1.06 | +3.13e-1 | HO-3 | −1.40 | +1.59 |
  | AR-4 | −2.26e1 | +2.51e1 | HO-4 | −7.53e-2 | +8.09e-2 |
  | AR-5 | −4.61e-1 | +5.53e-1 | HO-5 | −7.00e3 | +7.79e3 |
  | AR-6 | −3.58e-3 | +4.58e-3 | | | |

  The two solvers disagree with each other by the same order as their error, with opposite signs.
- **WP4 B2.** DEVSIM's AC analysis reproduced G_M bit for bit in 5 of 6 A-R states. It fails in extended precision too: in DEVSIM 2.11.2 the AC system is converted to double precision and solved once.
- **Grid data.** The P1 values at levels 0, 1 and 2 of every A-R state are in the released case files (`same_mesh_references`).
  - While scoping I looked at one sequence, AR-1 G_DD: 3.5434e-12, 3.5326e-12, 3.5276e-12. By hand that gives an observed order of about 2.2.
  - No other WP8 quantity has been computed.

## 3. Q1 design

**States.** The 11 three-terminal states of A-R and HO, on their level-1 meshes, body B1. These are the states and mesh of WP2 P-SS and WP4 B2. Each runs in a fresh DEVSIM process.

**Procedure per state** (`benchmark/wp8/ssmech.py`):
1. Build the WP2 deck and its zero-bias equilibrium. Compute the network reference G_net on the same mesh.
2. Extract DEVSIM's coupled Jacobian A at equilibrium and rhs = (b(V1 = 1 mV) − b(V1 = 0)) / 1 mV, exactly as WP2 P3 did. The output is G = −x[V1.I].
3. **M0, reproduction.** Recompute G_E and G_M, and compare them and G_net with the WP2 record.
4. **M1, exact solve of the same double-precision system.** Run mixed-precision iterative refinement from the G_E solution:
   - residuals r = rhs − A x are computed in double-double (`benchmark/wp8/ddlinalg.py`);
   - the inner solver is the same equilibrated SuperLU factorization;
   - the run is converged when |G_k − G_{k−1}| ≤ 1e-14 |G_k| at two consecutive iterations, with at most 40 iterations.
5. **M1b.** The same refinement with plain SuperLU (no equilibration) as the inner solver, as a check that M1 does not depend on the inner solver.
6. **M2, sensitivity to rounding of the Jacobian.** Build three systems in which every nonzero of A is multiplied by (1 + 2⁻⁵² u), with u ~ U[−1, 1] independently and seeds 20261106 + k. Solve each as in M1.

**Instruments.** The double-double primitives, the residual and the refinement are tested against exact rational arithmetic (`benchmark/tests/test_wp8.py`):
- the error-free transformations are exact;
- the residual is exact to 1e-29 of its scale;
- refinement recovers a tiny solution component that a double LU misses, and reports non-convergence for a useless inner solver.

**Quantities per state.**
- e_E = |G_E/G_net − 1| and e_M = |G_M/G_net − 1|; e_dbl = min(e_E, e_M), the better of the two double-precision solves.
- e_IR = |G_IR/G_net − 1| (M1).
- S = max over the three M2 samples of |G_k − G_IR| / |G_net|.
- The gate states are those with e_dbl > 1e-2. By section 2 these are expected to be the 10 states other than AR-6.

**Predictions.**
- **P-M0 (instrument).** M0 reproduces WP2's G_E, G_M and G_net to ≤ 1e-12 relative. A state that does not is excluded; if more than 2 states are excluded, Q1 is reported as undetermined.
- **P-M1 (the solve).** In at least all but one of the gate states, M1 converges and e_IR ≤ 1e-2 · e_dbl. That is, an exact solve of the same double-precision system removes at least 99% of the error.
- **P-M2 (the data).** In at least all but one of the gate states, S ≤ 1e-2 · e_dbl. That is, one-ulp rounding of the Jacobian entries moves G by less than 1% of the double-precision error.

Reason for the prediction: two backward-stable solvers produce errors of the same size with opposite signs. That points to the solve and not to the data, since the data are the same for both.

**Interpretation, fixed now.**

| P-M1 | P-M2 | What the manuscript will state |
|---|---|---|
| holds | holds | The error is the rounding of the double-precision factorization, amplified by the ill-conditioning of the output. The double-precision system determines G to e_IR, and refinement with double-double residuals recovers it. Remedy: residual-based refinement in higher precision. |
| fails | holds | The exact solution of the double-precision system is itself wrong, while one-ulp perturbations do not move it. The exported Jacobian differs systematically from the exact linearization, which is not a rounding effect. The cause is reported as not identified. |
| either | fails | Rounding the Jacobian entries alone moves G by at least 1% of the error. No double-precision solver can repair the output; it needs a higher-precision assembly and solve, or a cancellation-free formulation. |

Descriptive only, with no prediction:
- the number of states with e_IR ≤ 1e-6;
- S per state;
- the agreement of M1b with M1;
- the iteration histories.

## 4. Q2 design

**Data.** For each of the 6 A-R states and each of G_DD and G_SD (12 pairs), take the P1 values f3, f2, f1 at levels 0, 1, 2 from the case files.

**Procedure** (`benchmark/wp8/gci_wp8.py`): the three-grid procedure of Celik et al. (2008), which extends Roache's grid convergence index (1994) to non-constant refinement ratios.
- h ∝ N^−1/2, with N = nx·ny of the published mesh; r21 = √(N_L2/N_L1) and r32 = √(N_L1/N_L0) (about 1.41).
- Observed order p by fixed-point iteration; extrapolated value f_ext = (r21^p f1 − f2)/(r21^p − 1).
- GCI_fine = 1.25 |(f1 − f2)/f1| / (r21^p − 1), and the asymptotic ratio GCI_32/(r21^p GCI_21).
- Oscillatory convergence (opposite signs of the two differences) is reported as non-monotone.
- The procedure is tested on a manufactured sequence with known order and limit.

**Predictions.**
- **P-G1.** Monotone convergence in at least 10 of the 12 pairs.
- **P-G2.** Observed order 1 ≤ p ≤ 3 in every monotone pair.
- **P-G3.** The registered A-R tolerance contains the extrapolated value, |f1/f_ext − 1| ≤ tolerance, in all 12 pairs.

**Consequence.** The released tolerances are not changed (append-only). If P-G3 fails for any pair, the manuscript reports the GCI-based uncertainty beside the registered tolerance and names the pairs. A revision of the tolerance rule would need its own pre-registration.

Descriptive only: GCI_fine against the tolerance, the asymptotic ratio, and the independent finite-volume value (P4, its refined mesh) against f_ext.

## 5. Execution

- **Smoke.** One non-study state (dev1, V_G = 0.55 V, level 0, body B1; `run_wp8.py --smoke`) checks the pipeline: memory, time, convergence of the instruments. Its outcome is not evaluated or reported as a result.
- **Study.** `run_wp8.py` runs the 11 states one after another, each in a fresh process. A ledger line is written per run (`benchmark/wp8/results/ledger.jsonl`), with the record's SHA-256.
- **Evaluation.** `evaluate_wp8.py` applies sections 3 and 4 and writes `benchmark/wp8/results/evaluation_wp8.json`.
- **Afterwards.** A report (`benchmark/reports/WP8_*.md`), then the manuscript, updated from the new records by the asset script.
- **Untouched.** Every WP1–WP7 file, case file, result and the WP6 instruments.

## 6. Threats

- **M2 models rounding as one independent ulp per entry.** Real assembly errors can be several ulps and correlated. S is therefore an estimate of the sensitivity to rounding at one ulp, not a bound. If P-M2 holds only narrowly (S within a factor of 10 of the 1% gate), that is reported.
- **M1 is exact only to double-double accuracy** (a residual relative to its scale near 1e-32). With an output condition number up to about 1e21 relative to the true G, the attainable accuracy of G_IR is about 1e-11 relative. That is far below the 1e-6 level that matters here.
- **G_net and A come from the same double-precision equilibrium.** The network forms agree to ≤ 1e-9, so G_net is an adequate reference for errors above 1e-6.
- **This is the system DEVSIM's AC analysis solves.** The Jacobian is exported in double precision. In extended mode DEVSIM also converts the AC system to double, so this double-precision system is the one its AC analysis solves.

## 7. Procedure changes required by the smoke run (before any study run)

The smoke state (dev1, V_G = 0.55 V, level 0) is not a study state; its outcome is not evaluated. Its records are kept in
`benchmark/wp8/results/smoke/` (the first run under the name `*.run1_plain_superlu_M1b.json`).

1. **M1b inner solver.**
   - Plain SuperLU without equilibration, the registered M1b inner solver, diverged in the smoke: G grew by a factor of about 20 per iteration over 40 iterations. As an inner solver for this matrix it is useless, so it says nothing about M1.
   - M1b now starts from the UMFPACK solution and uses UMFPACK, factorized once through DEVSIM's interface, as its inner solver.
   - In the smoke this converged to the same G as M1 (relative difference 0).
2. **Convergence of M1 and M2.**
   - Besides the stability of G, a converged refinement must reach a final componentwise backward error in double-double of at most 1e-20.
   - The smoke reached 1.8e-27 (M1). This criterion is stricter than the registered one and was fixed before any study run.

Neither change touches a prediction, a threshold of section 3 or 4, or the interpretation table.

## 8. Amendment 1 (2026-11-06): exploratory sensitivity analysis X, no predictions

**When.** Added after the study records of AR-1 to AR-6 and HO-1 to HO-3 were seen (HO-4 and HO-5 still running), and before the registered evaluation was run. Sections 1 to 7 are unchanged, and nothing below changes a prediction, a threshold or the interpretation table.

**What was seen that motivates it.**
- In every record so far, the exact solve of the double-precision system (M1) misses G_net by an amount of the same size as the one-ulp spread S (M2).
- Within a device the absolute numbers do not depend on the state:
  - the M2 shift G_k − G_IR of each seed is the same to three digits in all states of the device (dev1: −2.23e-9, −5.71e-9, −5.07e-9 S/cm; dev2: −2.59e-9, −1.69e-8, −1.46e-8 S/cm);
  - the M1 error G_IR − G_net is within a factor of 2 per device (dev1: −4.4e-9 to −5.5e-9 S/cm; dev2: −0.9e-9 to +7.5e-9 S/cm).
- All M2 shifts so far have the same sign (6 independent draws: 3 seeds × 2 devices). A response linear in a symmetric perturbation gives that with probability 1/32. The realized perturbation is not exactly symmetric: rounding A_ij (1 + 2⁻⁵² u) to double biases it slightly.
- WP2 reported a componentwise condition number κ_G computed with the double-precision solution and adjoint and normalized by |G_E|. The M1 records show that this solution has a componentwise backward error up to 0.7, so that κ_G may be inaccurate.

**Questions** (descriptive):
- **X1.** Is the response of G to the M2 perturbations linear?
- **X2.** How large is the first-order sensitivity of G to the Jacobian entries, measured with accurate solutions? What absolute and relative spread does it predict?
- **X3.** Where is that sensitivity: in which equation block and which part of the device?

**Procedure** (`benchmark/wp8/xsens.py`, one fresh DEVSIM process per state; all 11 states, level 1, body B1):
1. **Build.** Build the state, G_net, A and rhs exactly as `ssmech.py`. Compute x (the M1 refinement) and the adjoint λ = A⁻ᵀ e_V1.I by the same refinement on Aᵀ, with the transposed equilibrated SuperLU as inner solver. Both are refined until the double-double componentwise backward error is at most 1e-24 (at most 12 iterations). Check: λ[V1.I] = x[V1.I], since both equal (A⁻¹)_{V1.I,V1.I}.
2. **First-order coefficients.** c_ij = λ_i A_ij x_j over the nonzeros, so that a perturbation δA changes G by δG₁ = Σ c_ij δA_ij / A_ij. Report:
   - κ_A = Σ|c_ij| / |G_net|, the componentwise condition number of G with respect to the entries;
   - the absolute scales u·Σ|c_ij| and u·(Σ c_ij²)^½ in S/cm, with u = 2⁻⁵³;
   - the same κ computed from the unrefined double-precision x and λ, normalized as in WP2.
3. **Linearity (X1).** For each M2 seed:
   - rebuild the perturbed matrix exactly as `ssmech.py` and take δA = A_k − A, which is exact in floating point;
   - compute δG₁ = λᵀ δA x and the second-order term δG₂ = −λᵀ δA z, where z = A⁻¹(δA x) is refined as in step 1;
   - re-solve the perturbed system as M2 did, record its backward error, and compare G_k − G_IR with δG₁ and δG₁ + δG₂.
4. **Location (X3).**
   - Per row i, s_i = Σ_j |c_ij|. Aggregate the s_i by equation block (silicon potential, silicon electron continuity, oxide potential, circuit) and by the class of the row's node: contact node (source, drain, body), n-type node on the drain side or the source side (net doping > 0, x above or below L/2), p-type node, oxide node.
   - Report the fractions of Σ s_i, the 10 largest rows with their position and doping, and a map of s binned on 64 × 32 cells of the silicon region.

**Output.** `benchmark/wp8/results/exploratory/XS_<id>_B1_lev1.json`, and a ledger line per run in the same ledger. This is reported as exploratory, separately from the registered results, in the WP8 report and the manuscript.

## 9. Amendment 2 (2026-11-06): check X4 of the cancellation, no predictions

**When.** Added after the X records of AR-1 to AR-3 were seen. Those records put 99.9% of Σ|c_ij| in the electron-continuity rows of the drain-side n-type region, with λ_i ≈ 1 there and the large c_ij in the potential columns. Sections 1 to 8 are unchanged.

**Proposed explanation, to be checked.**
- The Scharfetter–Gummel flux of an edge depends on the two potentials only through their difference.
- So in each electron-continuity row the potential-column entries sum to zero, up to the recombination term, whose potential derivative is small where holes are few.
- In the drain-side n-type region the potential response is nearly uniform: δψ ≈ 1 V/V and δn ≈ 0. The large c_ij of a row then cancel, and rounding the entries leaves a current that the drain collects (λ ≈ 1).

**X4** (`benchmark/wp8/xsens_rowsum.py`; AR-1 and AR-4, one state per device of the A-R tier, level 1, B1). Over the electron-continuity rows of the drain-side n-type region, report:
- the distribution (median, 99th percentile, maximum) of |Σ_{j∈ψ} A_ij| / Σ_{j∈ψ} |A_ij|, the relative row sum over the potential columns;
- the same quantity over the full row;
- the range of the refined δψ_j (x on the potential unknowns) and of δn_j / n_j in that region;
- the equilibrium n and the largest |A_ij| of the rows with the 10 largest s_i.

Descriptive; reported with X.
