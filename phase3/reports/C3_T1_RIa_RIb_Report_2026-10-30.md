# C-3 — Formal T1 for R-I(a), R-I(b) Boundary Verification and T1 Comparison

**Date:** 2026-10-30. **Branch:** `phase3/C3-2026-10-30` (base `c909bf1`, C-2).
**Scope:** WS3 and WS4 of the C-3 package. T1 is the zero-drain-bias consistency test (R04 v2.1 §5.2/§8.1, R05 v2.1 C-3/C-4, Research_08 §3 and §5).

**What T1 can and cannot show:**
- It verifies the **discrete** identity G_2D,h = (qμ_n/L)·Σᵢ n₀,ᵢwᵢ on each mesh level. This checks the deck's boundary conditions and its linear-response implementation.
- It does **not** show mesh convergence (reported separately, §4.1).
- It does not show the GCA at finite V_D.
- It says nothing about novelty.

All "E" values are sums of **estimates** (Research_08 §3), never bounds.

---

## 1. Configuration and provenance

| Item | Value |
|---|---|
| Lengths | L ∈ {1, 3, 10} µm |
| Mesh levels | k = 0–3 (D-1 maps, r = √2; actual whole-depth-axis h_rep ratios 1.4124–1.4147, used in §4.1) |
| Gate biases | V_G ∈ {0.8, 1.25, 2.0} V (subthreshold control, moderate inversion, strong inversion; V_T,cc = 1.043 V) |
| Models | R-I(a): lateral ψ = 1COL data + φ_n Dirichlet; R-I(b): lateral ψ Neumann + φ_n Dirichlet |
| Triple-point treatment | Explicit tensor mesh with one-node oxide contacts at the triple points (Research_08 §5) |
| Paths | L1 (exported Newton matrix, SciPy SuperLU, δ = 1 mV); L2 (DEVSIM small-signal AC, Re Y at 1 Hz and 0.1 Hz); L3 (finite-δ Neville extrapolation, δ = 0.25–8 mV, level 0 only, cross-check) |
| Runs | 12 runs (3 L × 4 levels) × 6 points each = **72 T1 points** (36 per model); L3 on 18 points |
| Pre-registered definitions | `phase3/c3/t1.py` module docstring, committed at `9dc3307` before any T1 result |

**Mesh levels (identical for all L except the lateral spacing):**

| Level | Lateral cells | Si depth nodes | First Si cell | Oxide nodes | Unknowns | Lateral/first-cell aspect (L = 1 / 3 / 10 µm) |
|---|---|---|---|---|---|---|
| 0 | 50 | 218 | 0.250 nm | 6 | 22,544 | 80 / 240 / 800 |
| 1 | 71 | 308 | 0.176 nm | 8 | 44,930 | 80 / 240 / 800 |
| 2 | 100 | 435 | 0.124 nm | 11 | 88,983 | 80 / 240 / 800 |
| 3 | 141 | 614 | 0.0878 nm | 15 | 176,508 | 81 / 242 / 808 |

**Provenance:**
- DEVSIM 2.11.2; Python 3.12.3; NumPy 2.3.5; SciPy 1.17.0.
- All 12 runs exited 0 (`results/phase3/C3/t1/RUNS.txt`). No point was dropped or re-run.
- L = 1 and 10 µm ran at commit `9dc3307` (clean tree).
- L = 3 µm ran at `752f024`, `e38146f`, `9b906bc` and `061c166`. These are flagged `git_dirty`. The per-commit diffs show that the only modified tracked files at those times were result logs of concurrently running jobs and one report. The solver code (`t1.py`, `phase3/c1/*`, `meshmap.py`, `m0.py`) differs from `9dc3307` only by the `body_corners` option, whose default (off) reproduces the original behavior.

## 2. Verdicts (pre-registered rule, unchanged)

Rule: **PASS** if |r₀| ≤ 3E and E ≤ 10⁻⁶. **NON-RESOLVING (class I)** if E > 10⁻⁶. **FAIL** otherwise.

| Path | PASS | NON-RESOLVING (class I) | FAIL | max \|r₀\| | max E |
|---|---|---|---|---|---|
| L1 | 71 | 1 | 0 | 2.68e-6 | 9.34e-6 |
| L2 | 71 | 1 | 0 | 9.07e-8 | 2.71e-6 |

**The single non-resolving point** is the same for both paths: L = 10 µm, level 3, R-I(a), V_G = 2.0 V (§3). **R-I(b) passes on all 36 points with both paths.**

### 2.1 Error components per configuration

Maximum over V_G of each quantity, as computed. E = E_eq + E_lin (L1) or E_freq + E_eq + E_lin (L2). For L1, E_freq = 0 by construction.

| L (µm) | Level | Model | L1: max \|r₀\| | L1: max E | L1: E_eq / E_lin | L2: max \|r₀\| | L2: max E | L2: E_freq / E_eq / E_lin |
|---|---|---|---|---|---|---|---|---|
| 1 | 0 | a | 1.5e-10 | 3.7e-10 | 4.9e-11 / 3.6e-10 | 7.5e-11 | 3.6e-10 | 9.0e-14 / 9.5e-11 / 2.7e-10 |
| 1 | 0 | b | 5.1e-10 | 1.9e-9 | 5.0e-10 / 1.4e-9 | 1.2e-10 | 6.0e-10 | 0 / 4.8e-11 / 5.5e-10 |
| 1 | 1 | a | 3.0e-10 | 5.8e-10 | 2.1e-10 / 5.1e-10 | 3.3e-10 | 6.9e-10 | 1.2e-13 / 5.7e-11 / 6.3e-10 |
| 1 | 1 | b | 2.3e-10 | 5.0e-10 | 1.8e-10 / 4.5e-10 | 2.9e-10 | 7.3e-10 | 0 / 1.6e-10 / 7.2e-10 |
| 1 | 2 | a | 5.3e-10 | 1.9e-9 | 4.0e-10 / 1.5e-9 | 3.1e-10 | 1.5e-9 | 9.9e-14 / 2.1e-10 / 1.3e-9 |
| 1 | 2 | b | 4.5e-10 | 3.9e-9 | 2.1e-9 / 1.8e-9 | 2.8e-10 | 1.3e-9 | 0 / 2.5e-10 / 9.9e-10 |
| 1 | 3 | a | 5.1e-9 | 1.2e-8 | 5.8e-9 / 1.1e-8 | 9.0e-10 | 6.3e-9 | 3.8e-13 / 2.4e-10 / 6.0e-9 |
| 1 | 3 | b | 9.4e-9 | 3.9e-8 | 9.4e-9 / 2.9e-8 | 9.7e-10 | 1.2e-8 | 0 / 3.4e-10 / 1.2e-8 |
| 3 | 0 | a | 7.5e-10 | 4.3e-9 | 5.8e-10 / 3.7e-9 | 2.7e-9 | 6.2e-9 | 5.4e-14 / 3.7e-10 / 5.9e-9 |
| 3 | 0 | b | 2.4e-9 | 2.2e-8 | 2.1e-8 / 3.4e-9 | 2.4e-9 | 3.6e-9 | 0 / 3.4e-10 / 3.6e-9 |
| 3 | 1 | a | 2.6e-9 | 9.1e-9 | 4.2e-9 / 4.9e-9 | 2.0e-9 | 4.7e-9 | 1.1e-13 / 5.2e-10 / 4.2e-9 |
| 3 | 1 | b | 4.1e-9 | 2.3e-8 | 4.3e-9 / 1.9e-8 | 2.4e-9 | 9.7e-9 | 0 / 3.8e-10 / 9.3e-9 |
| 3 | 2 | a | 2.2e-9 | 1.7e-8 | 5.9e-9 / 1.1e-8 | 5.2e-9 | 1.7e-8 | 1.3e-13 / 2.0e-9 / 1.5e-8 |
| 3 | 2 | b | 9.9e-9 | 4.6e-8 | 3.3e-8 / 1.3e-8 | 7.7e-9 | 1.5e-8 | 0 / 2.3e-9 / 1.3e-8 |
| 3 | 3 | a | 1.6e-8 | 4.7e-8 | 2.6e-8 / 3.5e-8 | 1.3e-8 | 2.7e-8 | 1.3e-13 / 7.9e-9 / 1.9e-8 |
| 3 | 3 | b | 4.5e-8 | 1.3e-7 | 4.9e-9 / 1.3e-7 | 6.8e-9 | 6.4e-8 | 0 / 2.6e-9 / 6.1e-8 |
| 10 | 0 | a | 2.8e-8 | 6.3e-8 | 4.9e-8 / 5.7e-8 | 9.7e-9 | 3.4e-8 | 4.0e-14 / 3.1e-9 / 3.1e-8 |
| 10 | 0 | b | 1.0e-8 | 3.2e-8 | 1.3e-8 / 1.9e-8 | 1.7e-8 | 4.2e-8 | 0 / 4.6e-9 / 3.7e-8 |
| 10 | 1 | a | 2.4e-8 | 5.2e-8 | 1.1e-8 / 4.3e-8 | 2.7e-8 | 5.6e-8 | 1.2e-13 / 1.6e-8 / 4.5e-8 |
| 10 | 1 | b | 1.4e-8 | 2.4e-8 | 2.8e-9 / 2.1e-8 | 1.2e-8 | 3.3e-8 | 0 / 4.1e-9 / 3.0e-8 |
| 10 | 2 | a | 8.1e-8 | 2.2e-7 | 7.3e-8 / 1.4e-7 | 6.1e-8 | 1.3e-7 | 1.5e-13 / 3.4e-8 / 1.1e-7 |
| 10 | 2 | b | 2.4e-7 | 7.1e-7 | 2.0e-7 / 5.1e-7 | 3.1e-8 | 2.6e-7 | 0 / 2.1e-8 / 2.5e-7 |
| 10 | 3 | a | **2.7e-6** | **9.3e-6** | 3.8e-6 / 5.6e-6 | 6.9e-8 | **2.7e-6** | 1.8e-13 / 2.5e-8 / 2.7e-6 |
| 10 | 3 | b | 6.8e-8 | 4.6e-7 | 2.7e-7 / 1.9e-7 | 9.1e-8 | 2.1e-7 | 0 / 1.8e-8 / 2.0e-7 |

Full per-point data: `results/phase3/C3/analysis_v2/t1_table.csv` (bitwise identical to v1).

**Observations:**
1. **Resolving power degrades with L and with level.** E grows by one to two orders of magnitude from L = 1 to 10 µm at fixed level, and by one to two orders from level 0 to 3 at fixed L. The lateral/first-cell aspect ratio scales with L (80 → 800) and the system size with level (22.5k → 176.5k unknowns). The linear-solve terms E_lin dominate, not E_freq or E_eq.
2. **E_freq is negligible** (≤ 4e-13). In R-I(b) it is exactly 0: Re Y(1 Hz) and Re Y(0.1 Hz) are bitwise equal. A plausible, unverified reason is that the drain contact carries only φ_n there, so the capacitive (ω²) part of Re Y is below round-off at 1 Hz.
3. **L2 is more accurate than L1 for large systems.** L2 is solved inside DEVSIM with UMFPACK. Its own r₀ stays ≤ 9.1e-8 everywhere.
   - But E_lin,L2 contains the pre-registered cross term |Re Y − G_L1|/G. It therefore inherits the L1 linear-solve error (§3).

## 3. The non-resolving point (L = 10 µm, level 3, R-I(a), V_G = 2.0 V)

| Quantity | Value (relative to G_GCA,h = 0.193932900 S/cm) |
|---|---|
| r₀ from L1 (SuperLU, exported matrix) | −2.68e-6 |
| r₀ from G_NS (DEVSIM UMFPACK Newton step) | −1.75e-8 |
| r₀ from L2 (Re Y at 0.1 Hz) | −1.75e-8 |
| L1 section sums (L/4, L/2, 3L/4) vs L1 contact current | +5.2e-7 to +5.6e-7 S/cm, i.e. up to 2.9e-6 relative |
| L1 repeat at τ_eq/100 (E_eq,L1) | 3.8e-6 (the state changes by < 1e-16 V; the L1 solve itself is not reproducible to better than ~1e-6) |
| E_L1 / E_L2 | 9.3e-6 / 2.7e-6 → both > 10⁻⁶ → **class I** |

**The verdict stands as computed. It is not re-evaluated** (Research_08 §3.5).

### 3.1 Post-hoc diagnostic of the cause (not used for any verdict)

`phase3/c3/l1_refine_diag.py`, `results/phase3/C3/l1_diag/run3/`. The diagnostic re-solves the same exported system.
- Step 0 reproduces the T1 value of G_L1 bit for bit.
- The metric is the normwise backward error of the row+column-equilibrated system, i.e. the B3 metric of the C-1(b) review.
- The two earlier diagnostic runs are kept:
  - `run1_normwise`: raw normwise metric, dominated by the largest rows (~1e-37), uninformative. Its L = 10 µm job was terminated before completion and is recorded as such.
  - `run2`: adds the componentwise metric, degenerate = 1 on decoupled contact rows.

| Solve of the T1 right-hand side b(δ) − b(0) | L = 1 µm, level 0 (control) | L = 10 µm, level 3 (class I point) |
|---|---|---|
| Plain SuperLU, as used by L1: equilibrated backward error | 9.1e-4 | 9.4e-3 |
| Plain SuperLU, as used by L1: r₀ | −1.5e-10 | −2.7e-6 |
| + 1/2/3 steps of fixed-precision iterative refinement: r₀ | −2.9e-11 / +3.1e-11 / −3.6e-11 | −6.4e-4 / −0.15 / −38 (**diverges**) |
| Row+column max-equilibration before SuperLU: backward error | 4.3e-16 | 2.9e-16 |
| Row+column max-equilibration before SuperLU: r₀ | −3.1e-10 | −2.4e-7 |
| DEVSIM UMFPACK Newton step (G_NS): r₀ | +4.2e-11 | −1.7e-8 |

**Reading:**
- The L1 path, as pre-registered, factors the unscaled exported matrix. For the T1 right-hand side, which is nonzero only in the drain-contact and circuit rows, this solve is **not backward stable**: equilibrated backward error 1e-3 to 1e-2. Fixed-precision iterative refinement cannot repair it at L = 10 µm, level 3; it diverges.
- Equilibrating before factorization restores backward stability at round-off. The remaining difference to UMFPACK (2.2e-7) is then a forward error at the conditioning of the problem.
- The class I point is therefore a limitation of the **L1 linear-algebra path**, not evidence of a deck or identity defect. The UMFPACK-based quantities (G_NS, L2) agree with the identity to 1.7e-8.
- This explanation is post hoc. A future T1 with an equilibrated L1 solve would be a **new, versioned E definition**, applied to all 72 points and pre-registered before execution.

## 4. The five comparisons (M0, 1COL, L1, L2, L3)

### 4.1 M0 vs 1COL: convergence of G_GCA,h to the continuum M0 conductance

R04 §7.4 with actual h_rep ratios (`results/phase3/C3/analysis_v2/G_mesh_convergence.csv`).
- G_GCA,h·L does not depend on L (same depth lines). The three lengths therefore give identical *relative* results; there are 3 distinct outcomes, not 9.
- The same classification applied to the 2-D L2 conductances is also asymptotic with the same orders (p₁₂₃ = 1.981–1.987).

| V_G | Classification | p₁₂₃ / p₂₃₄ | G_h,finest/G_M0 − 1 | U_disc (rel., estimate) | 2·U_disc test (T2 form) | Richardson − M0 (diagnostic) |
|---|---|---|---|---|---|---|
| 0.8 V | asymptotic | 1.987 / 2.013 | 1.49e-5 | 1.88e-5 | pass | −1.7e-7 |
| 1.25 V | asymptotic | 1.986 / 2.012 | 9.83e-6 | 1.24e-5 | pass | −1.2e-7 |
| 2.0 V | asymptotic | 1.981 / 2.000 | 1.39e-5 | 1.76e-5 | pass | −2.0e-7 |

**Analysis correction (v1 → v2):**
- The first analysis (`results/phase3/C3/analysis/`, commit `896d9d9`, kept unchanged) computed the observed order as ln(d₂₃/d₁₂)/ln r₁₂.
- That formula is exact only for a constant refinement ratio. With the actual, slightly non-constant ratios (1.4124, 1.4137, 1.4147) it overestimated p by about 0.01. v1 gave p₁₂₃ = 1.991–1.997, U_disc 1.24–1.87e-5, Richardson −4.7e-8 to −1.05e-7.
- v2 (`analysis_v2/`, commit `9a88e17`) solves the general three-level relation, including the q(p) term (Celik et al. 2008).
- The error was found by a new unit test on an exact h² sequence (`phase3/tests/test_c3_analysis.py`).
- Classifications and pass/fail outcomes are unchanged. The T1 verdict tables of v1 and v2 are bitwise identical.

**Continuum identity vs its numerical verification:**
- The continuum statement G = (qμ_n/L)∫n₀ dy for a laterally uniform V_D = 0 state is a mathematical consequence of the model equations.
- T1 checks its discrete analogue on each mesh (§2).
- This table checks that the discrete quantities converge, at second order, to the independently computed continuum value M0.
- After Richardson extrapolation the difference is −1.2e-7 to −2.0e-7. That is 60–110× below U_disc, but larger than M0's own estimated U_GCA,rel (1.7e-8, Package B B-3). The residual is therefore not explained by M0's tolerance. A plausible source is the higher-order terms neglected by the extrapolation (p₁₂₃ and p₂₃₄ differ by up to 0.026), but this is not demonstrated. It is a diagnostic and is not resolved further here.

### 4.2 1COL vs L1, 1COL vs L2, L1 vs L2

These are §2–§3 above:
- 1COL vs L1 is r₀,L1;
- 1COL vs L2 is r₀,L2;
- L1 vs L2 is the cross term inside E_lin,L2. It ranges from ~1e-10 (L = 1 µm, level 0) to 2.7e-6 (the class I point).

Diagnostic only: r₀ evaluated with DEVSIM's own Newton step (G_NS) is ≤ 9.1e-8 on all 72 points.

### 4.3 L3 cross-check (finite-δ extrapolation, level 0, 18 points)

| Result | Value |
|---|---|
| Status | 18/18 `ok` |
| \|r₀,L3\| | ≤ 1.1e-12 |
| E_L3 | 1.5e-13 to 2.1e-9 |
| Consistency with L2, \|r₀,L2 − r₀,L3\| ≤ 3(E_L2 + E_L3) | 18/18 |
| Largest consistency ratio diff/limit | 0.22 (L = 3 µm, R-I(b), V_G = 1.25 V) |

The finite-δ ratio follows r(δ) ∝ δ² where it is resolved (V_G = 1.25 and 2.0 V; log-slope 2.000–2.025). At V_G = 0.8 V, r(δ) is at round-off (|r| ≲ 1e-13), so the 0.8 V points in `T1_overview.png` (right panel) show noise, not a δ² law.

## 5. Boundary diagnostics (both models, all 72 points)

| Check (Research_08 §5) | Result |
|---|---|
| Potential jump at coincident oxide/Si pairs (triple points) | **0.0 V** on every point |
| Discrete Gauss closure \|Σ contact charges − Σ space charge\| / \|Σ space charge\| | ≤ **4.3e-14**. DEVSIM's contact-charge sign convention is not documented; the matching sign was the same on all 72 points, and the opposite sign gives ≈ 2, as expected of a pure sign |
| Single-entry (Dirichlet) oxide rows besides the gate | R-I(a): **2** (the two intended triple-point ties); R-I(b): **0** |
| Lateral uniformity of ψ at V_D = 0 | R-I(a): ≤ 1.1e-16 V; R-I(b): ≤ 1.7e-11 V (bottom corners, §6) |
| Lateral uniformity of n (relative) | ≤ 3.8e-15 |

### 5.1 Decoupling diagnostic (R-I(a) only; s_couple = 1 vs 0)

Pre-registered: |ΔG_L1|/G ≤ E_L1. **34/36 pass. 2 fail:**

| Point | \|ΔG\|/G | E_L1 | Ratio |
|---|---|---|---|
| L = 1 µm, level 3, V_G = 2.0 V | 1.27e-8 | 7.49e-9 | 1.70 |
| L = 3 µm, level 1, V_G = 1.25 V | 2.58e-9 | 2.34e-9 | 1.10 |

**Interpretation:**
- Both changes are at the level of the L1 linear-solve noise (§3.1). This suggests that E_L1 underestimates that noise somewhat.
- They are reported as diagnostic failures. The T1 verdict rule does not include this diagnostic, so the verdicts are unaffected.
- The largest decoupling change overall (2.3e-6, at the class I point) passes only because E_L1 is large there.

### 5.2 C-1(b) criteria B3 and B5 evaluated on the T1 runs

The C-1(b) review (`C1b_Final_Acceptance_Review_2026-10-30.md`) deferred B5 to this report.

**B5** — the linear-solve part of E_lin (|G_L1 − G_NS|/G plus the section-vs-contact difference) must be ≤ 1e-8. It is **met on 41/72 points and not met on 31/72:**

| L (µm) | Level 0 | Level 1 | Level 2 | Level 3 |
|---|---|---|---|---|
| 1 | 0/6 exceed | 0/6 | 0/6 | 2/6 |
| 3 | 0/6 | 1/6 | 2/6 | 3/6 |
| 10 | 6/6 | 5/6 | 6/6 | 6/6 |

**B3** — equilibrated backward error ≤ 1e-13. It was established for the Newton right-hand side. For the T1 right-hand side, the post-hoc diagnostic (§3.1) shows 9.1e-4 and 9.4e-3 on the two points examined. **B3 does not transfer to the L1 T1 solve.**

**Consequence:**
- The L1 path, as pre-registered, meets B5 on every point only for L = 1 µm at levels 0–2 and L = 3 µm at level 0 (24 points).
- Elsewhere it is an independent cross-check with a larger, reported linear-solve error, not a ≤ 1e-8 linear-response tool. This covers all of L = 10 µm, level 3, and part of L = 3 µm levels 1–2.
- The verdicts are unaffected, because this error is part of E by pre-registration.
- Addendum A to the C-1(b) review records this.

## 6. R-I(b) boundary verification (WS4)

Evidence: the exported Newton system at L = 1 µm, level 0, V_G = 2.0 V (`results/phase3/C3/rib_inspect/`), plus the T1 diagnostics on all 36 R-I(b) points.

| Requirement | Evidence | Result |
|---|---|---|
| **Ownership.** Which contact owns which equation | Contact equation lists: body = PotentialEquation (ψ = ψ_b0); gate = PotentialEquation; source and drain = ElectronContinuityEquation only (φ_n Dirichlet; no ψ equation); `oxc_s`/`oxc_d` = **no equations** (inert in R-I(b)) | As specified, with one exception: the bottom corners (below) |
| **Interface equations at the triple points** | Si node (0,0) and (L,0): the row has 2 nonzeros (ψ_Si − ψ_ox continuity). The Si bulk Poisson row is moved into the coincident oxide row, which has 7 nonzeros coupling oxide and Si neighbours, i.e. the summed displacement-flux balance | Interface continuity row **active** at both triple points (the Research_08 §5 requirement) |
| **Continuity** | Potential jump at every coincident pair | 0.0 V (36/36) |
| **Displacement flux** | No ψ equation on the lateral Si or oxide boundaries → zero normal D (natural condition). Discrete Gauss closure | ≤ 4.3e-14 (36/36) |
| **Conservation** | V_D = 0: section-sum vs contact current inside E_lin. Finite bias (preview, V_D = 50 mV): δ_KCL ≤ 4.4e-15, δ_cont ≤ 7.4e-15 | At round-off where the linear solve is accurate |
| **Refinement** | T1 on levels 0–3 at L = 1, 3, 10 µm | **36/36 PASS** for L1 and L2 |
| **No artificial constraints** | Single-entry oxide rows | **0** |
| **Lateral uniformity** | Max ψ spread vs level, L = 1 µm: 5.0, 3.4, 2.4, 1.7 (×1e-12 V); L = 10 µm: 17, 16, 14, 13 (×1e-12 V) | Located at the bottom corners; decreasing with refinement |

### Open item R-I(b)-1: bottom-corner ownership

**What happens:**
- With the default mesh, the two bottom corners (0, x_max) and (L, x_max) belong to the source/drain contacts, not to the body contact. In R-I(b) those contacts carry no ψ equation. The corner nodes therefore solve a Poisson row with natural conditions on both edges, instead of ψ = ψ_b0.
- Their ψ deviates from ψ_b0 by 5.0e-12 V at L = 1 µm, level 0.
- Adding the corners to the body contact (`--body-corners 1`) fixed the source corner only. The drain corner remained a bulk row. The exported matrix shows that DEVSIM assigns a node shared by two contacts to one of them only.

**Impact (estimate, not a bound):**
- body_corners = 1 vs default, 9 points (L = 1 µm levels 0–1, L = 10 µm level 0): Re Y changes by ≤ 7.5e-12 relative; T1 verdicts identical; ψ spread unchanged.
- This is a deviation from the stated R-I(b) bottom condition on two nodes. It is far below every E in this report.

**Fix (C-4):** end the lateral contacts one node above the bottom, so that the corners belong to the body. Then verify again from the exported matrix.

## 7. R-I(a) vs R-I(b): two specified reservoir models

Per Research_08 §7, S_res is a sensitivity between **these two** models, not a bound over admissible reservoirs.

**At V_D = 0 the two models agree.** Both reproduce the same laterally uniform equilibrium. |Re Y^(a) − Re Y^(b)|/G ≤ 5.4e-8 over all 36 pairs (largest at L = 10 µm, level 3, V_G = 2.0 V, where E_L2 sums to 2.9e-6). T1 cannot distinguish the reservoir models, by construction.

**At finite bias the models differ** (diagnostic preview, `fb_preview.py`, V_D = 50 mV). S_res(a,b) = I_a/I_b − 1:

| L (µm) | V_G = 0.8 V | V_G = 1.25 V | V_G = 2.0 V | Levels available |
|---|---|---|---|---|
| 1 | 1.24e-9 | 1.08e-3 | 6.08e-5 (level 2) | 0, 1, 2 (S_res changes ≤ 1.4 % from level 0 to 2) |
| 3 | 1.57e-10 | 1.32e-4 | 7.10e-6 | 0 |
| 10 | 1.47e-11 | 1.22e-5 | 6.56e-7 (level 1) | 0, 1 (changes ≤ 3 %) |

**Reading (preview):**
- S_res decreases about as L⁻² and peaks in moderate inversion.
- At L = 1 µm, S_res (1.1e-3) is about twice the magnitude of each model's own same-mesh discrepancy vs the GCA (−5.3e-4 for a, +5.5e-4 for b). The sign of "the GCA error" there is set by the reservoir model.
- This is the D1 framing risk recorded in the novelty matrix and in the finite-bias design (FB-1). It is not resolved by C-3.

## 8. Conclusions for WS3/WS4

1. **R-I(a) T1:**
   - 35/36 PASS for both L1 and L2. 1 class I point (L = 10 µm, level 3, V_G = 2.0 V), caused by the L1 linear-algebra path (§3.1).
   - Not a FAIL. Not reclassified.
   - Within the stated resolving power E ≤ 10⁻⁶, the discrete identity holds on every other configuration of the pre-registered matrix.
2. **R-I(b):**
   - 36/36 PASS. The boundary equations are verified from the exported matrix.
   - One documented two-node deviation at the bottom corners, with an estimated impact of ≤ 7.5e-12.
3. **Mesh convergence of G to M0:** asymptotic second order (p = 1.98–2.01) for all V_G. Differences are within 2·U_disc (ratio difference/U_disc = 0.79). Richardson difference 1–2e-7. These are analysis v2 values, after the order-formula correction in §4.1.
4. **Not met, or limits:**
   - B5 on 31/72 points.
   - B3 for the T1 right-hand side (2 points examined).
   - 2/36 decoupling diagnostics.
   - Without equilibration, the L1 path meets the ≤ 1e-8 linear-solve criterion on every point only for L = 1 µm levels 0–2 and L = 3 µm level 0.
5. **C-4 items:**
   - versioned L1 definition with equilibrated factorization, pre-registered and re-run on all points;
   - R-I(b) bottom-corner ownership fix;
   - R-I(c) implementation and T1;
   - the decoupling-diagnostic tolerance reconsidered as a versioned definition.

## Files

- **Code:**
  - `phase3/c3/t1.py` (driver, pre-registered definitions)
  - `phase3/c3/t1_analysis.py` (aggregation, R04 §7.4)
  - `phase3/c3/rib_inspect.py`
  - `phase3/c3/l1_refine_diag.py` (post-hoc)
  - `phase3/c3/fb_preview.py`
- **Raw results:**
  - `results/phase3/C3/t1/` (12 JSON + logs + RUNS.txt)
  - `results/phase3/C3/t1_rib_bodycorners/`
  - `results/phase3/C3/rib_inspect/`
  - `results/phase3/C3/l1_diag/run1_normwise|run2|run3/`
  - `results/phase3/C3/fb/`
- **Analysis:**
  - `results/phase3/C3/analysis_v2/` (cited): `t1_table.csv`, `t1_summary.json`, `G_mesh_convergence.csv`, `T1_overview.png/.svg`, `RUN.txt`
  - `results/phase3/C3/analysis/` (v1, kept; superseded for §4.1 only)
- **Tests:**
  - `phase3/tests/test_c3_analysis.py`
  - `results/phase3/C3/regression_final/`
