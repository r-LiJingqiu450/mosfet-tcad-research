# C-4 Results — Reservoir-Strip-Length Diagnostic and TC-1 Mechanism Test

**Date:** 2026-10-31 (UTC).
**Registered in:**
- `phase3/specs/C4_T1v2_Preregistration_2026-10-30.md` §4 (ℓ_res diagnostic);
- `phase3/reports/C4_Small_Bias_Perturbation_Analysis_2026-10-30.md` §4 (TC-1).

**Neither is FB-1.** No H-C decision is taken, and no σ is interpreted as a reservoir-independent error.

**Data:**
- `results/phase3/C4/fb4/` — 24 runs, all OK.
- `results/phase3/C4/tc1/` — 6 runs, all OK.
- `results/phase3/C4/analysis_post/` — CSV, summary, figure.

Runs used the bounded checkpointing runner (`phase3/c4/jobrunner.py`): 0 failures, 0 timeouts, 0 skips.

## 1. ℓ_res diagnostic

**Setting:** L = 1 µm, V_G = 1.25 V, V_D = 50 mV, levels 0–3.

**Metric:** ε_h = I_GCA,h/I_X,h − 1 (same-mesh, GL8). U_disc and extrapolation are per R04 §7.4, with the corrected observed-order formula.

**Conservation:** δ_KCL ≤ 1.7e-14 and δ_cont ≤ 2.3e-14 at all 24 states.

| Reference X | ε_h, levels 0 / 1 / 2 / 3 | Classification (p₁₂₃ / p₂₃₄) | ε extrapolated ± U_disc | Continuum metric ε_M0, extrapolated ± U |
|---|---|---|---|---|
| R-I(a) | −5.299 / −5.219 / −5.177 / −5.154 (×1e-4) | asymptotic (1.72 / 1.88) | **−5.126e-4 ± 3.5e-6** | −5.126e-4 ± 1.6e-5 |
| R-I(b) corrected | +5.518 / +5.612 / +5.657 / +5.680 (×1e-4) | asymptotic (1.98 / 2.08) | **+5.703e-4 ± 2.9e-6** | +5.705e-4 ± 1.6e-5 |
| R-I(c) f = 0.1 | +2.299 / +2.360 / +2.388 / +2.402 (×1e-4) | asymptotic (2.01 / 2.15) | **+2.416e-4 ± 1.8e-6** | +2.417e-4 ± 1.5e-5 |
| R-I(c) f = 0.3 | +5.578 / +6.202 / +6.502 / +6.650 (×1e-5) | asymptotic (2.02 / 2.11) | **+6.80e-5 ± 1.8e-6** | +6.81e-5 ± 1.5e-5 |
| R-I(c) f = 1 | −2.119 / −1.029 / −0.452 / −0.165 (×1e-5) | asymptotic (2.01 / 1.83) | **+1.2e-6 ± 3.6e-6** (not resolved: \|ε\| < 2U) | +1.3e-6 ± 1.6e-5 |
| R-I(c) f = 3 (**exploratory**) | −7.24 / −4.43 / −2.55 / −1.44 (×1e-5) | **monotone non-asymptotic** (1.53 / 1.14) | no extrapolation; finest −1.44e-5, U = 1.7e-4 (inconclusive) | asymptotic, +8e-7 ± 3.2e-5 |

**Pairwise sensitivity** S(X,Y) = I_X/I_Y − 1 between specified models, level 3:

| Pair | S |
|---|---|
| a vs b | +1.084e-3 |
| a vs c(f = 1) | +5.14e-4 |
| b vs c(f = 1) | −5.69e-4 |
| c(0.3) vs c(1) | −6.8e-5 |
| c(1) vs c(3) | −1.3e-5 (still changing with level) |

**Reading** (within the tested family and at this single bias point only):
1. **Resolved and mesh-independent.** For a, b, c(0.1) and c(0.3) the discrepancies are 60–200× their U_disc. The same-mesh and continuum metrics agree after extrapolation to ≤ 0.2 %.
2. **Strongly reservoir-dependent.**
   - Within the gated-strip family, ε falls monotonically from +5.70e-4 (f → 0, which is R-I(b)) through +2.42e-4 and +6.8e-5 to an unresolved value at f = 1 (|ε| < 7e-6).
   - R-I(a) sits at −5.13e-4.
   - The spread across admissible ideal reservoirs, ~1.1e-3, is **more than 100× the smallest member**.
3. **The reservoir-symmetric average has no special status.** σ = (ε^a + ε^b)/2 = +2.9e-5 ± 3.2e-6 (half the sum of the two U_disc) is resolved. The longest resolved gated strip gives ε^c(1) = +1.2e-6 ± 3.6e-6, which is not. So σ is not the value approached by the "long reservoir" member of the family, and **it cannot be read as a reservoir-independent reduction error**.
4. **Long-strip behaviour is not established.**
   - f = 3 is non-asymptotic: the strip mesh is 10·√2^k cells over 0.91 µm, i.e. 91 nm at level 0, which is too coarse.
   - Whether ε^c(ℓ) converges to ~0, or to a small nonzero value, as ℓ/W_dep → ∞ is **open**.
   - A future test needs a fixed strip cell size and longer strips.

## 2. TC-1 (mechanism test of Proposition 2)

**Setting:** V_G = 1.25 V; V_D ∈ {0.5, 1, 2, 4, 8} mV; levels 0–2; L ∈ {1, 3} µm. Weighted cubic fit, σ scaled by χ²/dof.

| L | Ends | c₁ (V⁻¹), levels 0 / 1 / 2 | \|c₁\|/σ (min over levels) | c₂ (V⁻²), level 2 |
|---|---|---|---|---|
| 1 µm | ab | 0.033595 / 0.033567 / 0.033539 | 6.2e5 | +0.0876 |
| 1 µm | ba | −0.033595 / −0.033567 / −0.033539 | 7.7e5 | −0.0767 |
| 1 µm | aa | −1.3e-7 / −1.4e-7 / −1.3e-7 | **3.1–3.7** | −0.1633 |
| 1 µm | bb | +1.9e-7 / +1.7e-7 / +1.8e-7 | **3.1–4.0** | +0.1764 |
| 3 µm | ab | 0.0037369 / 0.0037356 / 0.0037341 | 5.2e5 | +0.0079 |
| 3 µm | ba | −0.0037369 / −0.0037356 / −0.0037341 | 4.9e5 | −0.0103 |
| 3 µm | aa | −2.5e-8 / −2.2e-8 / −2.2e-8 | **3.3–3.6** | −0.0217 |
| 3 µm | bb | +2.2e-8 / +2.1e-8 / +2.4e-8 | **3.3–3.7** | +0.0193 |

**Registered verdict: `INCONCLUSIVE`.**
- The rule requires |c₁| ≤ 3σ for aa and bb at every level. It is exceeded slightly, at 3.1–4.0σ, at every level.
- It is not "falsified", because no symmetric deck reaches 10σ.
- The asymmetric part of the rule is met: >10σ everywhere, and the level-1→2 change of c₁ is 0.04–0.08 %.

**Post-hoc observations** (labelled as such; they do not change the verdict):
- The symmetric-deck linear term is 5–7×10⁻⁶ of the asymmetric one at every level and both L. It is consistent with the truncation bias of a cubic fit (χ²/dof ~ 10⁷) rather than a physical O(V_D) term.
- c₁(ab) + c₁(ba) = 0 to ≤ 2e-6 relative, as the mirror relation of the end-layer problem implies.
- c₁(ab) scales as L^−1.999 (levels 0–2). This matches the end-layer scaling of Proposition 3.
- c₂ of the symmetric decks scales as L^−1.70 → −1.84 with refinement.
- c₂ is approximately additive in end contributions: c₂(ab) + c₂(ba) − c₂(aa) − c₂(bb) is ≤ 1.4 % of |c₂(aa)| at L = 1 µm.

**What a decisive TC-1 would need** (for a future, newly registered test): a "zero linear term" criterion robust to fit-model truncation. For example:
- a quartic fit, or a smaller-δ ladder; and
- a bound expressed relative to the asymmetric c₁, e.g. |c₁,sym| ≤ 10⁻³·|c₁,asym| at every level.

That rule is **not** applied here, because it was formulated after seeing the data.
