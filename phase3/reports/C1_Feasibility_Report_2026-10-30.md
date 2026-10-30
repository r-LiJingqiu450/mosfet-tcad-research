# Package C-1 — DEVSIM Feasibility Report: (a) AC, (b) Jacobian export, (c) current extraction

**Date:** 2026-10-30. **Specification:** Research_05 v2.1, C-1 (Research_07 Part III).
**Branch / commit:** `phase3/C1-2026-10-30`. Code at `f4ca0bf`. Reproduction runs on a clean tree are in `results/phase3/C1/repro_f4ca0bf/`.
**Scope:** software and numerical feasibility on one mesh level. **No T0, T1 or T2 test was executed**, and none is claimed. The "identity preview" numbers below are diagnostics on a C-1 deck that does not meet the R04/R05 mesh or uncertainty requirements.

**Final status: `C1_PARTIAL`.**
- All three feasibility questions are answered **yes**.
- Two items keep C-1 short of READY:
  1. The literal C-1(b) acceptance figure (≤ 1e-10 against DEVSIM's own update) is missed by about 2× because of cross-solver forward error. The figure is not numerically justified for this system (§3).
  2. A DEVSIM mesh-ownership defect at the Si/SiO₂/lateral-contact triple points leaves a 1.5–5.2e-6 identity residual. That is above E_req = 1e-6, and it must be closed before C-3 (§5).

## 1. Environment and solver parameters

| Item | Value |
|---|---|
| DEVSIM | 2.11.2 (wheel SHA-256 `fef14e3f…c0a8f793`). `info`: extended_precision True; direct_solver `custom` = UMFPACK 5.1 through DEVSIM's stock Python callback; BLAS suffix `_` |
| Math libraries | `DEVSIM_MATH_LIBS=/lib/x86_64-linux-gnu/libblas.so.3:/lib/x86_64-linux-gnu/liblapack.so.3` |
| Python environment | venv 3.12.3: NumPy 2.3.5, SciPy 1.17.0 (pinned) |
| Units | cm. Contact currents are A/cm per unit depth; I′ [A/µm] = I × 1e-4; G′ [S/µm] = G × 1e-4 |
| Physics (PM checklist, printed into every JSON) | F0: N_A = 1e16 cm⁻³, t_ox = 20 nm, x_max = 1.2 µm, L = 1 µm. T = 300 K, V_T = kT/q, n_i = 1e10 cm⁻³, ε 11.7/3.9, Boltzmann statistics, μ_n = 1350 cm²/Vs, no G–R, holes frozen p = n_i e^(−ψ/V_T), ψ intrinsic-referenced, ψ_gate = V_G + ψ_b0 (ψ_b0 = −0.3571586 V) |
| Mesh | DEVSIM built-in 2-D tensor mesher. 21 x-columns (dx = 50 nm). 228 Si y-lines (0.5 nm at the interface, 20 nm at the bottom). 13 oxide lines. Thin inert padding regions are needed for edge contacts. 9851 unknowns |
| Contacts | gate (oxide top, ψ). body (Si bottom, corners excluded, ψ only, so electron flux is natural = 0). source/drain (Si lateral, full depth): n Dirichlet φ_n = bias, plus R-I(a) ψ = PsiLat + s_couple·bias. Drain integrated into circuit node `drain_bias` with `V1` (acreal = 1 V) |
| **Newton stopping** | R05 D-3 absolute rule: ‖Δψ‖∞ ≤ 1e-10 V and V_T·max\|Δn\|/n ≤ 1e-10 V, then 2 further iterations. Attained: \|Δψ\| ≈ 1e-16 V and \|Δφ_n\| ≈ 1e-16 V |

**Why DEVSIM's own stopping test was replaced.** DEVSIM computes rel = \|Δx\|/(\|x\| + 1e-10) (source `Equation.cc`; minError is not exposed to Python). The intrinsic-referenced ψ crosses 0 in the depletion layer, so a fully converged Newton can be reported as "not converged". This happened above about 1.0 V in attempts 2 and 3 of C-1(b), with rel ≈ 1e-6 at the zero-crossing node. DEVSIM's device-level absolute error also mixes V and cm⁻³ in one maximum.

## 2. C-1(a) — small-signal AC with a circuit voltage source on a custom contact: **YES**

Script `phase3/c1/c1a_ac.py`. Final run `repro_f4ca0bf/c1a` (exit 0, 8.0 s); attempts 1–3 are kept.

**API behavior.**
- `solve(type='ac', frequency=f)` needs the contact equations to be declared with `circuit_node`, including the node-model derivative w.r.t. that node.
- The AC source is `circuit_element('V1', acreal=1)`.
- Results are read from `get_circuit_node_value(solution='ssac_real'|'ssac_imag', node='V1.I')`.
- AC node solutions appear as `Potential_Real/_Imag` and `Electrons_Real/_Imag`.
- Sign: Y = −I_V1 / 1 V.

| Quantity (V_G = 0.8 / 1.25 / 2.0 V) | Result |
|---|---|
| Re Y(f) vs f | Flat for f ≤ 100 Hz. E_freq(1 Hz, 0.1 Hz) = 9.0e-15 / 9.2e-15 / 4.4e-14 |
| Im Y | Exactly ∝ f (C·ω) from 1e-3 Hz to 1 MHz |
| Re Y(0.1 Hz) vs L1 (direct solve with the exported matrix) | −4.1e-11 / 1.1e-10 / 9.8e-11 |
| L1 vs DEVSIM's own first Newton step / δ | ≤ 1.1e-10 (δ = 1 mV and 1 µV); RHS exactly linear in δ |
| Decoupling diagnostic (s_couple 0 → 1, lateral ψ data follow V_D) | ΔG_L1/G = −1.6e-11 / 1.4e-11 / 1.3e-11. Im Y at 1 MHz changes (9.9e-9 → 9.8e-5 S/cm at 0.8 V), as expected |
| G′_2D (S/µm) | 3.2036e-9 / 3.4761e-5 / 1.9401e-4 |

The earlier smoke-test discrepancy (golden AC conductance off by 1.75e-6 in `ssac_diode.py`) was **not re-run**. On this deck, AC equals L1 to ≤ 1.1e-10 once DC is converged to 1e-16 V. That supports, but does not prove, the DC-stopping-point hypothesis.

## 3. C-1(b) — `get_matrix_and_rhs(format='csc')` vs the Newton system: **YES** (literal 1e-10 figure missed by ~2×; reported, not forced)

Script `phase3/c1/c1b_jacobian.py`. Final run `repro_f4ca0bf/c1b` (exit 0, 14.5 s).
- Attempt 1: exit 1 (JSON tuple keys).
- Attempt 2: invalid; the equilibrium had not converged.
- Attempt 3: exit 1 (ill-posed DEVSIM stopping test).
- Attempt 4: valid, but its backward-error metric was degenerate.

All attempts are kept with notes.

**Method.**
- A *recording* solver callback captures the exact (A, b) that Newton factors. DEVSIM's default direct solver is already a Python callback (UMFPACK), so this changes no physics.
- Variable updates are re-declared `default`, i.e. undamped: x ← x + Δ (`Equation.cc::DefaultSolutionUpdate`). No damped update is compared.
- Two states are tested at V_G = 1.25 V, starting from the converged V_D = 0 state: V1 stepped to 1 mV and to 50 mV.

| Property | Result (both states) |
|---|---|
| Returned keys | `static`, `dynamic` (charge matrix for AC/transient), `permutation` (`row`, `copy`), `format` |
| Static matrix vs the matrix Newton factors | Same sparsity pattern (nnz 77 217); **max \|ΔA\| = 0.0** |
| Static RHS vs Newton RHS | **max \|Δb\| = 0.0** |
| Does the export change the state? | No: 0.0 in every unknown block, including circuit nodes |
| RHS sign | b = −F: source sets `rhssign = −1` (`Newton.hh`). Confirmed: undamped DEVSIM update = **+A⁻¹b** (relative difference 6e-14 / 3e-16); the opposite sign gives a relative difference of 2.0 |
| Linearization state | A = ∂F/∂x at the *current* state. Central-difference columns of −b agree to ≤ **2.56e-10** at 23 sampled unknowns (interior, interface, drain/source/body contact, gate, oxide, circuit `drain_bias`). Flipped sign: 2.0 |
| Unknown ordering | Oxide Potential 0–272, Si Potential 273–5060, Si Electrons 5061–9848, circuit `V1.I` 9849, `drain_bias` 9850. Covers 0…n−1 exactly. Valid only during a solve (captured in the callback) |
| Contact rows | Bulk rows at contact nodes are replaced. Permutation: 971 moved, 952 dropped (= 2×228 + 2×228 + 19 + 21 contact rows), 0 copied. Drain electron rows have 3 nonzeros (n, ψ, `drain_bias`); diagonal = NodeVolume × ∂(node model)/∂n, i.e. contact node models are integrated over NodeVolume |
| External update (SciPy SuperLU on the exported A, b) vs DEVSIM's own update (UMFPACK) | **1.4–1.9e-10** relative per unknown block. R05 target: ≤ 1e-10 |
| Normwise backward error (row/column-equilibrated system) | UMFPACK 1.7e-14 / 3.3e-15; SuperLU 6.4e-16 / 5.3e-16 |
| 1-norm condition estimate | 5.5e35 raw (unit disparity between V and cm⁻³); **9.8e10 after equilibration** |

**Assessment.**
- The question C-1(b) poses is whether the static matrix is the Newton Jacobian including contact rows. The answer is **exactly yes**.
- The ≤ 1e-10 agreement figure compares two *different* backward-stable direct solvers. For κ ≈ 1e11 the forward difference can reach κ·ε ≈ 2e-5. The observed 1.9e-10 is well inside that, so the figure is **not a numerically justified acceptance threshold**. When the update is computed through the same factorization path (the callback), it reproduces to 1e-13.
- Proposed restatement, for a reviewer to decide: bitwise identity of matrix and RHS with the Newton system, plus backward error ≤ 1e-13, plus the cross-solver difference reported as part of E_lin.

## 4. C-1(c) — contact currents vs SG edge-flux sums: **YES**

Script `phase3/c1/c1c_currents.py`. Final run `repro_f4ca0bf/c1c` (exit 0, 2.3 s). V_G = 1.25 V. Cross-sections are mid-cell lines at x = 0.025, 0.275, 0.525, 0.725 and 0.975 µm.

| Quantity | V_D = 1 mV | V_D = 50 mV |
|---|---|---|
| I_drain (`get_contact_current`, electrons) | +3.4607e-4 A/cm = 3.4607e-8 A/µm | +1.33492e-2 A/cm = 1.33492e-6 A/µm |
| δ_cont (5 sections vs drain) | **1.7e-13** | **7.8e-15** |
| δ_KCL (source + drain + displacement terms) | **6.4e-14** | **6.5e-16** |
| \|V1.I + I_drain\| / \|I_drain\| | 6.0e-15 | 2.6e-16 |
| Diagonal edges crossing a section | 227 per section; EdgeCouple ≤ 1e-36, i.e. zero on right-triangle tensor meshes | same |

**Sign map (V_D > 0).**
- Drain contact current is +I_D and source contact current is −I_D.
- V1.I = −I_D.
- The edge model F_e (`ElectronCurrent`) is conventional current density from n0 to n1 (A/cm²).
- Σ s_e F_e·EdgeCouple_e, with s_e = +1 for edges oriented in +x, equals **−I_D**: current flows drain → source.

**AC (V_D = 0).** The linearized edge-flux sum is Σ_v (∂F/∂v@n0·δv0 + ∂F/∂v@n1·δv1)·EdgeCouple, using DEVSIM's own derivative edge models. Its real part matches Re Y to ≤ 3.0e-11 at 0.1 Hz and ≤ 4.1e-9 at 1 MHz. The imaginary parts vary along x (distributed jωQ) and are not expected to equal Im Y.

Caveat: the finite-V_D states keep the V_D = 0 lateral ψ data (s_couple = 0). They are current-extraction tests, **not** R-I(a) physics states.

## 5. Deck defect found: Si/SiO₂ interface lost at the triple points

**What happens.**
- DEVSIM's built-in 2-D mesher assigns a mesh point to only one contact or interface, and the **first declaration wins**. Four orderings were tested in `order.py` (scratch).
- The Si nodes at (0, 0) and (L, 0) belong to the source/drain contacts. Their interface pairs are therefore dropped: the mesh log shows "19 node pairs" for 21 columns.
- The oxide corner node then has no continuity link to Si.

**Effect on the T1 identity preview** (r̂₀ = G_2D/G_GCA,h − 1; L1, L2 and the Newton step agree on it to ≤ 1e-10):

| Lateral ψ data | Oxide treatment | Lateral ψ spread at the interface row | r̂₀ preview (0.8 / 1.25 / 2.0 V) | Run |
|---|---|---|---|---|
| Copied from a 2-D Neumann run (itself corrupted at the corners) | none | 0.18–0.23 V | −0.90 / −0.92 / −0.94 | `c1a/attempt1` |
| Independent 1COL FV (`onecol.py`, identical y-lines) | none | 2.4–4.4 mV | 2.1e-2 / 9.6e-3 / 5.8e-3 | `c1a/attempt2` |
| Independent 1COL FV | oxide node next to each triple point pinned to 1COL | 0.6–1.1 µV | **5.2e-6 / 2.5e-6 / 1.5e-6** | `c1a/attempt3`, `repro` |

**Notes on the pin.**
- The "pin" is a one-segment oxide contact. It captured only the oxide node at y = −0.243 nm (deviation from 1COL exactly 0.0); the corner point itself stays owned by the Si contact and floats at +5.38 mV at V_G = 1.25 V (`repro_f4ca0bf/oxide_pin_check_VG1p25.txt`).
- `get_element_node_list` returns no element nodes for it.
- The floating corner oxide node remains, and it accounts for the residual.
- In the interior column, 2-D and 1COL ψ agree to 1.4e-9 V (Si) and 1e-10 V (oxide).

**Consequence.** With the residual at 1.5–5.2e-6 > E_req = 1e-6, a C-3 T1 run on this construction would fail as **class P** (boundary-condition implementation). Closure options, in order of preference:
1. DEVSIM `custom_equation`, which assembles rows by hand (the shipped `examples/diode/diode_1d_custom.py` uses it), to restore the corner continuity.
2. An explicit mesh (`create_gmsh_mesh`) to test whether corner ownership can be shared.
3. A documented deviation (e.g. R-I(c) reservoir strips).

R-I(b) (Neumann lateral ψ) has the same corner issue and was not previewed.

## 6. T1 method capability matrix

| T1 route | Feasible? | Evidence (this report) | Measured error components on the C-1 deck | Remaining blocker for a T1 *pass* |
|---|---|---|---|---|
| **L1** direct linear response (`get_matrix_and_rhs`) | **YES** | §3: export ≡ Newton system; G_L1 = Newton-step G to ≤ 1.1e-10 | E_freq = 0. E_lin (cross-solver / section sum) ~1e-10. Backward error ≤ 2e-14 | Triple-point defect (§5) |
| **L2** small-signal AC on the custom circuit contact | **YES** | §2 | E_freq ≤ 4.4e-14 at f₁ = 1 Hz. Re Y vs L1 ≤ 1.1e-10. Decoupling ≤ 2e-11 | Triple-point defect (§5) |
| **L3** finite-V_D ratio extrapolation | **Not exercised** | Prerequisites met: DC Newton reaches \|Δψ\| ≈ 1e-16 V; δ_cont ≤ 1.7e-13 at 1 mV | — | Needs the R-I(a) ψ data at finite V_D (ψ₁(x; V_G, V_D) from 1COL) and the triple-point fix |

Because (a) and (b) are both feasible, R05's L3-only fallback does **not** apply.

## 7. Recommendation on Package C-2

**C-2 can begin**, under these conditions.
1. **Rebuild D-1** (`phase3/meshmap.py`, x-map per R04 §7.3, r = √2, K ≥ 3) with its tests first. It was lost with the 2026-10-29 files.
2. **Avoid the triple-point defect in the 1COL deck.** Use a one-cell-wide 2-D mesh with **no lateral contacts**: Poisson only, with n = n_i e^((ψ−φ_n)/V_T) at a prescribed **uniform** φ_n.
   - With equal lateral φ_n data and zero-flux top and bottom, continuity is satisfied identically. Lateral φ_n Dirichlet is therefore redundant.
   - The interface keeps both corner pairs.
   - This departs from the literal R05 C-2 boundary-condition wording and needs reviewer acknowledgment.
   - DEVSIM node positions must be checked to coincide with D-1 coordinates (explicit mesh lines).
3. Inputs: M0 ψ_s (intrinsic-referenced) and N_s from `phase3/m0.py` / `m0_F0.csv`. Stopping: `deck.newton` (R05 D-3). Threshold via `vt_cc_tabulated` on a ≤ 0.01 V grid.

**C-3/C-4 are not ready.** They need the triple-point closure (§5) and a decision on the C-1(b) acceptance wording (§3).

## 8. Known limitations and open issues

1. GitHub push was not available: the request to elevate to push access was denied by the session's safety classifier and was not retried. Preservation is by a verified git bundle and a source ZIP.
2. The C-1 deck is one mesh level (not D-1), with no uncertainty analysis. The previews are not T1.
3. The `ssac_diode.py` 1.75e-6 golden AC difference was not re-run in this session.
4. The C-1(b) acceptance figure needs a reviewer decision (§3).
5. Package A remains open (Phase 2D / G-0 not on the remote). The D_pre calibration was not rebuilt.
6. The test suite writes `_test_common_npz.npz` into `results/golden/` (`test_common.py:430`). This is documented and isolated: the file is untracked, excluded from commits, the bundle and the ZIP, and `stage3.json` is unchanged. See `results/phase3/C1/regression/scratch_file_isolation.txt`.

## 9. Regression and golden integrity

`pytest semiconductor_sim phase3/tests` gives **106 passed, 4 skipped, 4 warnings, exit 0** in both environments (94 production + 12 Package B). Golden `stage3.json` SHA-256 `17ca32745b69c34355113347eb0f51f95772cf8a76f76b1decdd85925ef67828` is the same before and after. `git diff 7073f20 -- semiconductor_sim` is empty.
