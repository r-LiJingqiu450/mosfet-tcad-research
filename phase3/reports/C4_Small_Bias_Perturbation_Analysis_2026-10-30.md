# Small-Drain-Bias Perturbation Analysis of the GCA Discrepancy in Ideal-Reservoir References

**Date:** 2026-10-30 (C-4, Task F). **Status:** the project's own derivation.

**Not checked against Ward 1988/1990.** Their full texts were not obtainable in this environment. The result is elementary and may well be known. **No originality is claimed for it.** Its purpose is to explain the numerical observations and to give a sharp, falsifiable test (TC-1).

## 1. Setting

**Notation:** x is depth and y runs along the channel. The channel is 0 < y < L.

**Equations** (F0 physics: holes frozen, constant μ_n, no generation/recombination):

  −∇·(ε∇ψ) = q(p − n − N_A) in Si,  −∇·(ε_ox∇ψ) = 0 in oxide,  ∇·(n∇φ_n) = 0 in Si,
  n = n_i e^{(ψ−φ_n)/V_T},  p = n_i e^{−ψ/V_T}.

**Boundary conditions:**
- ψ = ψ_g at the gate and ψ = ψ_b0 at the bottom;
- zero normal electron flux at the interface and the bottom;
- φ_n = 0 on the source reservoir set and φ_n = V_D on the drain reservoir set.

**Mirror-symmetric reservoir class:**
- Geometry, doping and permittivities are symmetric under y ↦ L − y.
- The drain reservoir set is the mirror image of the source set.
- The electrostatic end conditions at the two ends are of the same type:
  - homogeneous Neumann (R-I(b), outer ends of R-I(c)); or
  - Dirichlet data equal to the 1-D slice solution at the local φ_n (R-I(a)).
- R-I(a), R-I(b) and R-I(c) (all ℓ_res) belong to this class.

**Terminal current.** By conservation and the zero-flux top/bottom conditions, the current I is the same through every cross-section of the channel. Hence

  I = (qμ_n/L) ∫∫_channel n ∂φ_n/∂y dx dy.   (1)

**Expansions in V_D:**
- φ_n = Σ V_D^k u_k;
- ψ = ψ₀ + Σ V_D^k ψ_k;
- n = n₀(1 + V_D η₁ + V_D² η₂ + …), with η₁ = (ψ₁ − u₁)/V_T and η₂ = (ψ₂ − u₂)/V_T + η₁²/2.
- All u_k with k ≥ 2 vanish on the reservoir sets, because φ_n there is exactly 0 or V_D.

## 2. Results

**Proposition 1 (first order; R04 §5.2).** ∇·(n₀∇u₁) = 0 with n₀ = n₀(x) gives u₁ = y/L in the channel. Therefore I₁ = (qμ_n/L)N_s(V_G, 0) for every reservoir in the class. This is the T1 identity.

**Proposition 2 (second order in the current).** From (1), the V_D² term is

  I₂ = (qμ_n/L) ∫∫ [n₀ ∂u₂/∂y + n₀η₁/L] dx dy = (qμ_n/L²) ∫∫ n₀ η₁ dx dy.

The first part integrates to n₀[u₂]₀^L = 0.

ψ₁ solves the linearized Poisson problem

  −∇·(ε∇ψ₁) + (q/V_T)(n₀ + p₀)ψ₁ = (q/V_T) n₀ u₁   (oxide: −∇·(ε_ox∇ψ₁) = 0; ψ₁ = 0 at the gate and bottom).

Let g(x) = ∂ψ_slice/∂φ_n be the 1-D linearized slice response. Then ψ₁^G = g(x)·u₁(y) solves the channel PDE exactly, because ∂²u₁/∂y² = 0.

Write ψ₁ = ψ₁^G + δ₁. Subtracting the laterally uniform part g/2 leaves a problem whose data are antisymmetric under y ↦ L − y; uniqueness then makes δ₁ antisymmetric. Hence ∫∫ n₀ δ₁ = 0, and

  I₂ = (qμ_n/L²) ∫∫ n₀ (g − 1) (y/L)/V_T = (qμ_n/L)·½ dN_s/dφ_n = I₂^GCA.

**Consequence:** for every mirror-symmetric ideal reservoir, I_2D − I_GCA = O(V_D³), i.e.

  ε(V_D) := I_GCA/I_2D − 1 = c^X(V_G, L)·V_D² + O(V_D³).

**Proposition 3 (structure of c^X; scaling argument, not a proof).** The V_D³ current contains two kinds of terms:
- terms quadratic or bilinear in δ₁: ∫∫ n₀(g−1)u₁δ₁ and ∫∫ n₀δ₁²/2;
- the second-order potential correction δ₂.

δ₁ is an end-layer function. Its lateral extent ℓ_x is set by the electrostatics (W_dep, t_ox) and is independent of L for L ≫ ℓ_x; its amplitude is ~ g·ℓ_x/L. δ₂ also contains an interior part forced by ∂²ψ₂^G/∂y² = g₂/L², of relative size (W/L)².

Both the reservoir-dependent end contribution and the interior contribution therefore scale as L⁻² at leading order: **c^X = [A(V_G) + B^X(V_G)]·W²/L² + o(L⁻²)**. B^X depends on the reservoir through δ₁, which differs between R-I(a), (b) and (c).

**Implications** (if Proposition 3 holds):
1. No reservoir-independent "GCA reduction error" is defined at leading order. The interior part A is an outer-expansion construct and cannot be separated from B^X by L-scaling.
2. The average of two reservoir models (σ in the FB-1 design) has no special status. **This supports the C-4 instruction's caution.**
3. FB-1's hypothesis H-C ("a reservoir-independent σ exists") is not well posed as stated.

## 3. Consistency with the numerical evidence so far

| Observation | Proposition | Status |
|---|---|---|
| T1: G_2D = G_GCA on every mesh, R-I(a)/(b) (C-3), (c) (C-4) | 1 | consistent (verified to E ≤ 1e-6) |
| r(δ) ∝ δ² with log-slope 2.000–2.025, no linear term, both models, L = 1 and 10 µm (C-3 L3 ladders) | 2 | consistent |
| S_res(a,b) and α ∝ L^−1.95 (1.25 V, L = 1 → 10 µm) | 3 | consistent |
| ε_h at 50 mV, L = 1 µm, 1.25 V: a −5.3e-4, b +5.5e-4, c(ℓ = 0.3 W_dep) +5.6e-5 (level 0) | 3 (B^X reservoir-dependent) | consistent; level-0 only for c |

**A mechanism test is needed.** These agreements are post-hoc: Proposition 2 was derived after the data were seen. TC-1 tests something the data did not show.

## 4. TC-1 — pre-registered falsification test of the mechanism

**Prediction.** Proposition 2 depends on mirror symmetry. A **mirror-asymmetric** reservoir pair makes δ₁ lose its antisymmetry. That gives ∫∫ n₀δ₁ ≠ 0 and hence a **linear term**: ε = c₁V_D + O(V_D²) with c₁ ≠ 0.

**Decks.** All four use the same corrected mesh (C-4 `b` geometry, dual bottom corners). Only the lateral ψ treatment differs:

| Deck | Source end | Drain end |
|---|---|---|
| aa | Dirichlet ψ = 1COL(φ_n) | Dirichlet ψ = 1COL(φ_n) |
| bb | Neumann | Neumann |
| ab | Dirichlet (1COL(0)) | Neumann |
| ba | Neumann | Dirichlet (1COL(V_D)) |

The triple-point oxide contact follows its end: Dirichlet tie or inert.

**States:**
- L ∈ {1, 3} µm; V_G = 1.25 V; levels 0–2.
- V_D ∈ {0.5, 1, 2, 4, 8} mV.
- Newton τ = 1e-13 V, repeated at 1e-15 V.
- GCA side: 4-point Gauss–Legendre over 1COL slices (3-point value kept as the quadrature estimate).

**Estimator.** Weighted least squares r(δ) = c₁δ + c₂δ² + c₃δ³ on the five points. Weights 1/e_i², where e_i = |r(τ) − r(τ/100)| + |I_GCA,4 − I_GCA,3|/I_GCA + 1e-14. The parameter uncertainty is the LS covariance scaled by max(1, χ²/dof).

**Decision rules (fixed now):**
- **Supported** if, at every level and both L:
  - |c₁| ≤ 3σ(c₁) for aa and bb; and
  - |c₁| > 10σ(c₁) for both ab and ba;
  - with c₁ converged across levels (relative change from level 1 to 2 ≤ 10 %).
- **Falsified** if either:
  - a resolved linear term appears for aa or bb (|c₁| > 10σ at two or more levels); or
  - ab and ba show no resolved linear term (|c₁| ≤ 3σ at every level).
- **Inconclusive** otherwise.

Under every outcome, no claim beyond this mechanism test is made.

### Registration note (added at commit, before any registered TC-1 run)

The decision rules above were written before the one development run below and are committed unchanged. That run was a level-0 smoke test of `phase3/c4/tc1.py` at L = 1 µm (`results/phase3/C4/dev_smoke/tc1/`). It is kept, but **it is not part of the registered TC-1 matrix**. Its result:

| Deck | c₁ (V⁻¹) | \|c₁\|/σ | c₂ (V⁻²) |
|---|---|---|---|
| ab | +3.36e-2 | 6.7e5 | +0.083 |
| ba | −3.36e-2 | 8.4e5 | −0.082 |
| aa | −1.3e-7 | 3.4 | −0.168 |
| bb | +1.9e-7 | 4.0 | +0.171 |

Under the rule as written, the symmetric decks at level 0 alone (3.4σ and 4.0σ > 3σ) already rule out a "supported" outcome. The rule is **not** changed. The registered runs are evaluated with it, whatever that implies.
