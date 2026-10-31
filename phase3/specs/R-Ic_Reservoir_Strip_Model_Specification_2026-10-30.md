# R-I(c) — Reservoir-Strip Reference Model: Specification (C-4 Task D)

**Date:** 2026-10-30. **Source requirement:** R04 v2.1 §4.3 ("Optional R-I(c) family: reservoir strips of length ℓ_res ∈ {0.1, 0.3, 1}·W_dep at each end, with uniform φ_n and natural ψ conditions. It broadens the sampled set of specifications. It still does not produce a bound").
**Implementation:** `phase3/c4/deck4.py` (model `c`); registered in `phase3/specs/C4_T1v2_Preregistration_2026-10-30.md` §2.

## 1. Definition

**Coordinates.** DEVSIM x runs along the channel and y is depth (the R04 (y, x) pair with axes swapped).
- Channel: 0 < x < L.
- Source strip: −ℓ_res ≤ x ≤ 0.
- Drain strip: L ≤ x ≤ L + ℓ_res.

**ℓ_res.**
- ℓ_res = f·W_dep, with W_dep = √(2ε_Si·2ψ_F/(qN_A)) = 0.30393 µm (`meshmap.wdep_max_cm`).
- Pre-specified: f ∈ {0.1, 0.3, 1}.
- f = 3 is used only in the ℓ_res diagnostic, labelled exploratory.

**Material and physics (identical to R-I(a)/(b)):**
- Si with N_A = 10¹⁶ cm⁻³ in the strips too.
- Holes frozen; constant μ_n = 1350 cm²/Vs; Boltzmann statistics; no G–R; ψ intrinsic-referenced.

**Oxide and gate.** The 20 nm oxide and the gate (ψ = V_G + ψ_b0) cover the strips. The oxide lateral faces at x = −ℓ_res and x = L + ℓ_res are homogeneous Neumann.

**Electron condition in the strips.** At every Si node of a strip above the bottom row, φ_n = 0 (source) or φ_n = V_D (drain):
- Contact `source`/`drain` carries the electron equation n = n_i exp((ψ − bias)/V_T). The drain strip is on circuit node `drain_bias`.
- These contacts carry **no** ψ equation, so ψ in the strips is the bulk Poisson solution ("natural ψ").

**Strip bottom row.** Contact `bs`/`bd` carries ψ = ψ_b0 (R04 §4.2 bottom condition) **and** φ_n = bias.

**Channel bottom.** ψ = ψ_b0 (contact `body`), with zero normal electron flux.

**Outer lateral Si faces** (x = −ℓ_res, L + ℓ_res): ψ homogeneous Neumann. φ_n is set by the strip condition.

**Triple points.** These are at the outer ends. The one-node oxide contacts carry no equation, and the Si interface-continuity row is active (Research_08 §5, R-I(b) rule).

**Channel length.** L is the distance between the inner strip edges, as for R-I(a)/(b). G_GCA,h uses the same L and the same 1COL depth lines.

**Terminal current.** −V1.I, equal to the electron current into `drain` + `bd`. KCL is checked against `source` + `bs`.

## 2. Mesh

- Channel lateral map: identical to R-I(a)/(b) (uniform, 50·√2^k cells).
- Each strip: uniform, N_k = round(10·√2^k) cells.
- Depth and oxide maps: D-1 maps, identical for all models.

Strip spacing at level 0: 3.0 nm (f = 0.1), 9.1 nm (f = 0.3), 30.4 nm (f = 1), 91 nm (f = 3). At the channel/strip junction the lateral spacing jumps, e.g. 20 nm → 3 nm at L = 1 µm, f = 0.1. This is accepted. Convergence is checked by refining both together.

## 3. Interpretation decided before use, and its limits

**Gated strips.** R04 does not say whether the gate covers the strips. It is taken to cover them, so that:
- the V_D = 0 state is laterally uniform; and
- the zero-bias identity holds with channel length L, which T1 then verifies.

Ungated strips would add a built-in lateral transition at V_D = 0. That makes the reference structurally different (closer to R-II), and it is not tested.

**Limits.** R-I(c) is one more **specified** reservoir model:
- Its results, together with R-I(a) and R-I(b), describe sensitivity **within the tested family**.
- They do not bound the effect of other admissible reservoirs (doped n⁺ regions, finite reservoir conductivity, non-equilibrium injection).
- None of a, b, c is "the" physical reference.

**Limit cases:**
- f → 0 reproduces the corrected R-I(b) node for node (strip of one column plus the dual bottom corner).
- f → ∞ is a semi-infinite, perfectly conducting, gated reservoir.

## 4. Verification required (C-4)

- Contact equation lists.
- Bottom-row equations from the exported matrix.
- Continuity (0 V), Gauss closure, single-entry oxide rows (expected 0), lateral uniformity at V_D = 0.
- T1 v2 on L ∈ {1, 3, 10} µm, levels 0–3 (f = 0.3), and L ∈ {1, 10} µm (f = 0.1, 1).
- T0/T2 compatibility: the laterally uniform 2-D state equals the 1COL solution. The T0/T2 comparison of 1COL with M0 is established in C-2.
- Finite-bias conservation (δ_KCL, δ_cont) in the ℓ_res diagnostic.
