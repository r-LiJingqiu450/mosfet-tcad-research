# Package C-2 — One-Column Reference (1COL), Tests T0/T2, and the D-1 Mesh Family

**Date:** 2026-10-30. **Specification:** Research_05 v2.1 C-2 and D-1; criteria from R04 v2.1 §7.4 and §8.1.
**Branch / commit:** `phase3/C2-2026-10-30`. Code and results at `fe44d68`; bit-identical reproduction is in `results/phase3/C2/repro_fe44d68/`.
**Nature of the evidence:** this is **solution verification of a code against a semi-analytic reference of the same continuum problem** (M0). It is not model validation, says nothing about real devices, and makes no novelty claim.

**Final status: `C2_READY`** (§8).

## 1. D-1 mesh-family generator (restored by independent replacement)

The 2026-10-29 `phase3/meshmap.py` was **not recoverable**. The search covered the container filesystem, every commit in Git history, and the Project documents. `phase3/meshmap.py` is a new implementation with 36 tests (`phase3/tests/test_meshmap.py`).

**Construction.**
- Each axis is a sequence of segments between **fixed breakpoints**.
- Each segment is a fixed monotone map X(ξ): uniform, or geometric X = a + (b − a)(Q^ξ − 1)/(Q − 1).
- Level k samples ξ uniformly with N_s,k = round(N_s,0 · √2^k) cells per segment.

**Properties that hold by construction and are tested:**
- breakpoints are nodes on every level;
- the same map is used on all levels;
- monotone coordinates and correct cell counts;
- level-0 specification values.

**Properties that do not hold and are not required:**
- **Adjacent-cell grading is level-dependent** (Q^(1/N_k): 1.0220, 1.0155, 1.0109, 1.0077 for the geometric segment, tested). R05's D-1 test text ("grading ratios identical across levels") contradicts this; per the 2026-10-30 instruction it is not imposed. *Specification conflict recorded for the reviewer.*
- The effective refinement ratio per segment is N_k+1/N_k = 1.400–1.429, not exactly √2 (rounding). It is recorded per level and segment.
- No convergence order is assumed.

**F0 maps.**

| Axis | Breakpoints | Segments |
|---|---|---|
| Si depth | 0, 3·W_dep = 0.9118 µm, 1.2 µm | geometric, N₀ = 202, Q = 81.27, first cell 0.25 nm; then uniform, N₀ = 15, 19.2 nm |
| Oxide | −20 nm, 0 | uniform, N₀ = 5 |
| Lateral | 0, L | uniform; N₀ = 50 per R05 C-3, and N₀ = 20 used in the corner study |

**Per-level quality** (`results/phase3/C2/D1/mesh_quality.csv`; all coordinates in `mesh_family_F0.json`, SHA-256 `56e9ab76…22c6f4`):

| Level | Si nodes | h_min (nm) | h_max (nm) | Max adjacent ratio | Ratio across 3·W_dep joint | Oxide cells |
|---|---|---|---|---|---|---|
| 0 | 218 | 0.250 | 19.88 | 1.0347 | 0.966 | 5 |
| 1 | 308 | 0.176 | 14.09 | 1.0264 | 0.974 | 7 |
| 2 | 435 | 0.124 | 9.99 | 1.0403 | 0.961 | 10 |
| 3 | 614 | 0.088 | 7.08 | 1.0321 | 0.969 | 14 |

At the interface the oxide spacing is about 16× the Si spacing on every level (16.0, 16.2, 16.1, 16.3). This is harmless for the charge-free 1-D oxide, whose linear potential is represented exactly.

## 2. 1COL deck and physics configuration (`phase3/c2/onecol_deck.py`)

| Item | Setting |
|---|---|
| Physics | PM checklist of R05 §0.2, F0. Constants shared with the C-1 deck: V_T = kT/q (constants.py), n_i = 1e10 cm⁻³, ε 11.7/3.9, μ_n = 1350 cm²/Vs, N_A = 1e16 cm⁻³, t_ox = 20 nm, x_max = 1.2 µm, holes frozen p = n_i e^(−ψ/V_T), ψ intrinsic-referenced, ψ_b0 = −0.35715858 V |
| Mesh | Explicit (`create_gmsh_mesh`) at **exact D-1 coordinates** (no inserted nodes); width 1 nm |
| Gate | Contact on the oxide top: ψ = V_G + ψ_b0 (V_FB = 0). Checked: error 0.0 V on all states |
| Bottom | Contact on the Si bottom (all bottom nodes): ψ = ψ_b0. Electrons: no bottom contact, so zero normal flux. Checked: error 0.0 V |
| Interface | `continuous` ψ. Every pair is active: no lateral contact solves PotentialEquation, so the C-1 corner rule does not apply |
| **Formulation P (primary)** | Unknown ψ only; n = n_i exp((ψ − φ_n)/V_T) with a **spatially uniform** φ_n parameter; one cell wide; laterals homogeneous Neumann |
| **Formulation F (equivalence check)** | Unknowns ψ and n; SG electron continuity (flux identical to the C-1 deck and to DEVSIM `simple_dd.py`); **two cells wide**, so the middle column is a genuine interior where continuity is solved; lateral electron contacts n = n_i exp((ψ − φ_n)/V_T) (φ_n Dirichlet, equal on both sides); ψ Neumann |
| Newton stopping | R05 D-3 absolute rule: ‖Δψ‖∞ ≤ τ and V_T·max\|Δn\|/n ≤ τ, τ = 1e-10 V, plus 2 further iterations; U_iter from τ = 1e-12 V. Attained \|Δψ\| ≈ 1e-17 to 1e-16 V |
| N_s | N_s,FV = Σ_i n_i·w_i over one Si column, with w_i the 1-D dual widths (trapezoid weights incl. half cells at x = 0 and x_max) — R04 condition D4. Units cm⁻² |
| Q_s | Q_s,FV = q Σ_i (p_i − n_i − N_A) w_i (C/cm²); a diagnostic |
| Runs | V_G = 0 → 2 V in 0.05 V steps at φ_n = 0; φ_n continuation (0.05 V steps) to 0.05 and 0.5 V at V_G ∈ {0.8, 1.25, 2.0} V; V_G = V_T,long = 0.99635 V (constants.VT_LONG) for T0. Levels 0–3; zero failed states |

## 3. Equivalence of the reduced formulation (P) to the zero-current, frozen-hole formulation (F)

**Scope of the claim.** The equivalence holds under these boundary conditions: equal φ_n on both lateral edges, zero normal electron flux at the interface and bottom, no generation/recombination. Under them the discrete continuity problem has the unique solution φ_n ≡ const, and the SG flux of a Boltzmann equilibrium vanishes identically. It is **not** a statement about the 2-D drift–diffusion problem at V_D ≠ 0.

Numerically, over 47 common states × 4 levels:

| Check | Max over all states and levels |
|---|---|
| \|ψ_s,P − ψ_s,F\| | **1.1e-16 V** |
| \|N_s,P/N_s,F − 1\| | **6.7e-16** |
| F: electron-continuity residual per node / sum of \|flux terms\| | 6.0e-15 |
| F: SG edge flux / its largest term | 1.1e-14 (P evaluated the same way: 9.8e-15) |
| F: quasi-Fermi relation \|n/(n_i e^((ψ−φ_n)/V_T)) − 1\| | 8.2e-15 |
| P and F: hole relation \|p/(n_i e^(−ψ/V_T)) − 1\| | ≤ 1e-14 |
| F: lateral ψ spread | 1.2e-13 V |
| Gate / bottom BC errors | 0.0 / 0.0 V |
| Gauss \|D_ox + Q_s,FV\|/\|Q_s\| (V_G > 0) | 2.4e-11 – 3.4e-10 |

The Gauss imbalance is the displacement flux leaving through the bottom contact: the depletion tail at x_max, of order e^(−(x_max − W_dep)/L_D) ≈ e^(−21). It quantifies the finite-depth vs semi-infinite specification difference relative to M0, and it is below every uncertainty in §4.

## 4. T0 and T2 results (`results/phase3/C2/t0t2/t0t2.csv`, `c2_t0t2.json`)

M0 comes from `phase3/m0.py` (Package B). ψ_s,M0 = φ_s + ψ_b0. U_GCA combines the |default − tightened ÷ 100| difference with the `quad` estimate (N_s) or the brentq xtol (ψ_s). At V_G = 0 M0 refuses V_G < 0.05 V by design; flat band is exact there (φ_s = 0 is the root of V_FB + φ_s − Q_s/C_ox = 0 with Q_s(0) = 0).

U_disc follows R04 §7.4 on levels 3 (finest) … 0 with no assumed order; U_iter = \|f(1e-10) − f(1e-12)\|.

### T0 — surface potential ψ_s (V, intrinsic-referenced)

| V_G (V) | M0 | 1COL L3 | Signed diff | Class (p₁₂₃ / p₂₃₄) | U_disc | U_iter | U_GCA | R05 limit (U+1e-6) | R04 limit (2U+1e-6) | Richardson − M0 (diagnostic) | Verdict R05 / R04 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | −0.357158576 | −0.357158576 | 0.0 | converged to round-off | 0 | 0 | 0 (analytic) | 1e-6 | 1e-6 | 0 | PASS / PASS |
| 0.99635 | 0.357158576 | 0.357157776 | −8.00e-7 | asymptotic (1.99 / 2.02) | 1.00e-6 | 0 | 2e-14 | 2.0e-6 | 3.0e-6 | +3.6e-9 | PASS / PASS |
| 2.0 | 0.512750248 | 0.512738061 | −1.22e-5 | asymptotic (1.98 / 2.01) | 1.53e-5 | 0 | 2e-14 | 1.63e-5 | 3.16e-5 | +8.9e-8 | PASS / PASS |

### T2 — electron sheet density N_s (cm⁻²; includes the neutral background n_b·x_max as M0 does)

| V_G (V) | φ_n (V) | M0 | Rel. diff L3 | Class (p₁₂₃/p₂₃₄) | U_disc,rel | U_GCA,rel | U_iter (abs) | R05 limit | R04 limit | Richardson rel. − M0 | Verdict R05 / R04 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.8 | 0 | 1.48026e7 | +1.49e-5 | asym. (1.99/2.02) | 1.87e-5 | 2.6e-14 | 1.9e-9 | 1.87e-5 | 3.73e-5 | −6.7e-8 | PASS / PASS |
| 0.8 | 0.05 | 2.14061e6 | +1.49e-5 | asym. | 1.87e-5 | 9.5e-14 | 0 | 1.87e-5 | 3.73e-5 | −6.7e-8 | PASS / PASS |
| 0.8 | 0.5 | 5.90061e-2 | +1.49e-5 | asym. | 1.87e-5 | 1.1e-13 | 0 | 1.87e-5 | 3.73e-5 | −6.7e-8 | PASS / PASS |
| 1.25 | 0 | 1.60661e11 | +9.83e-6 | asym. | 1.24e-5 | 1.1e-14 | 0 | 1.24e-5 | 2.47e-5 | −4.7e-8 | PASS / PASS |
| 1.25 | 0.05 | 1.09921e11 | +1.06e-5 | asym. | 1.34e-5 | 1.1e-14 | 0 | 1.34e-5 | 2.67e-5 | −5.1e-8 | PASS / PASS |
| 1.25 | 0.5 | 9.06325e4 | +2.43e-5 | asym. | 3.05e-5 | 1.1e-14 | 0 | 3.05e-5 | 6.09e-5 | −1.1e-7 | PASS / PASS |
| 2.0 | 0 | 8.96605e11 | +1.39e-5 | asym. (1.98/2.01) | 1.75e-5 | 2.4e-14 | 0 | 1.75e-5 | 3.50e-5 | −1.05e-7 | PASS / PASS |
| 2.0 | 0.05 | 8.35431e11 | +1.37e-5 | asym. | 1.72e-5 | 2.3e-14 | 0 | 1.72e-5 | 3.44e-5 | −9.9e-8 | PASS / PASS |
| 2.0 | 0.5 | 3.08790e11 | +1.25e-5 | asym. | 1.57e-5 | 1.6e-10 | 0 | 1.57e-5 | 3.13e-5 | −6.7e-8 | PASS / PASS |

**Reading the results.**
- Every sequence is in the asymptotic range, with observed order 1.98–2.02. No order was assumed.
- The finest-level differences are ≈ 0.8·U_disc. The passes therefore rest on the R04 safety factor 1.25 (R05 form) or on 2×1.25 (R04 form).
- The sharper evidence is the Richardson diagnostic, which is not a reported value: the extrapolated 1COL results agree with M0 to 3.6e-9 V / 8.9e-8 V (ψ_s) and ≤ 1.1e-7 relative (N_s). That is 100–200× smaller than the finest-level discretization error, so the discrete solutions converge onto the M0 solution of the same continuum problem.
- The residual ~1e-7 is consistent with higher-order discretization terms. It is far above the finite-depth effect (~1e-10) and the M0 quadrature uncertainty (≤ 1.6e-10).
- U_iter is 0 or ≤ 1.9e-9 cm⁻² (relative 1e-16), so T5 (U_iter ≤ 0.1·U_disc) holds.

**Specification inconsistency, recorded.** R05 C-2 and R04 §8.1 state T0/T2 tolerances that differ by a factor 2 (U vs 2U). Both are evaluated; every point passes both.

**Full sweep, diagnostic** (`sweep_L3_vs_M0.csv`, figure `C2_sweep_L3_vs_M0.png`): over V_G = 0.05–2.0 V at φ_n = 0 on L3, max \|Δψ_s\| = 1.22e-5 V and max \|ΔN_s/N_s\| = 1.74e-5. Convergence plots are in `C2_convergence.png`.

## 5. Execution provenance (all exit 0 unless noted)

| Step | Command | Exit |
|---|---|---|
| 1COL P levels 0–3, τ = 1e-10 and level 3, τ = 1e-12 | `python phase3/c2/onecol_deck.py --level k --form P --tau τ --out results/phase3/C2/onecol` | 0 (×5) |
| 1COL F levels 0–3, τ = 1e-10 and level 3, τ = 1e-12 | same with `--form F` | 0 (×5) |
| V_T,long points, P levels 0–3 and L3 at τ = 1e-12 | `… --points [[0.9963499970604937,0.0]] --tag vtlong` | 0 (×5). First attempt **exit 1** (×5): invalid JSON (`np.float64` repr); logs kept as `*_attempt1_FAILED.log` |
| Analysis | `python phase3/c2/c2_analysis.py` | 0 |
| Reproduction at `fe44d68` | all of the above plus corner study and C-1(b) study (24 runs) | 0 (×24); every numerical value bitwise identical (`repro_fe44d68/COMPARISON.txt`) |

`RUNS.txt` files list exact commands with start and end times. DEVSIM 2.11.2; venv Python 3.12.3, NumPy 2.3.5, SciPy 1.17.0.

## 6. Regression and protected files

`pytest semiconductor_sim phase3/tests`: **142 passed, 4 skipped, 4 warnings, exit 0** in both the pinned and current environments (94 production + 12 Package B + 36 D-1).
- Golden `stage3.json` SHA-256 `17ca3274…ef67828` is unchanged before and after.
- `git diff 7073f20 -- semiconductor_sim` is empty.
- The test-suite scratch file `results/golden/_test_common_npz.npz` is still untracked and excluded from commits and archives.

## 7. Limitations

1. T0/T2 verify the DEVSIM 1COL discretization against M0 for the same physics. They do not validate the physics, and they say nothing about 2-D effects.
2. The P ↔ F equivalence is shown only for uniform φ_n with zero-flux top and bottom. It is not a claim about finite-V_D 2-D states.
3. Only F0 was run. F1/F2 are future work per R05.
4. The finest level (h = 0.088 nm at the interface) is a numerical device. No physical (quantum) significance is implied (Boltzmann, no quantum correction, by PM design).
5. U values are estimates, not bounds (R04 §7.1). The Richardson values are diagnostics only.
6. The D-1 replacement deliberately deviates from R05's "identical grading across levels" test wording (§1).

## 8. Status

**`C2_READY`.** The one-column physics comparison is validated against M0 for the specified points:
- T0 at 3 gate biases and T2 at 9 (V_G, φ_n) points pass both published criteria, with asymptotic classification;
- iteration and quadrature uncertainties are negligible;
- the equivalence of the reduced formulation is demonstrated under the stated boundary conditions;
- the D-1 generator is restored with tests.

The C-6 convention audit and the C-3/C-4/C-5 tasks are **not** part of this package and were not started.
