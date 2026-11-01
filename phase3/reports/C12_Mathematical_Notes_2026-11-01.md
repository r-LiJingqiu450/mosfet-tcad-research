# C-12 Mathematical Notes — the Liang–Taur Model, the End Correction and Two-Sided Current Bounds

**Date:** 2026-11-01. **Companion to:** `C12_MultiGate_Literature_Assessment_2026-11-01.md`.

**Notation.**
- [LT04] and [L06] are the Liang–Taur paper and Liang's dissertation (see the assessment, §1).
- Symmetric double gate, Δφ₁ = Δφ₂ = Δφ, gate insulators of thickness t_i and permittivity ε_i, film thickness t_si.
- k = π/λ₁, V_T = kT/q.
- s = πt_si/(2λ₁), so that kx ranges over [−s, s] inside the film.
- A₀ = E_g/2q + Δφ − V_g is the source barrier above the long-channel potential.

**Numerical checks** are in `phase3/c12/lt_endcorr.py` (K1–K7, K9) and `phase3/c12/lt_vt_rolloff.py` (K8). Tests are in `phase3/tests/test_c12_lt_identities.py`.

---

## N1. The LT "conjugate function" is the ε-weighted eigenfunction

Separating u = X(x)·e^{±ky} in ∂ₓ(ε∂ₓu) + ∂_y(ε∂_yu) = 0 gives the Sturm–Liouville problem

(εX′)′ + k²·ε·X = 0,  X = 0 at both gate electrodes,  X and εX′ continuous at the interfaces.

Its weight function is ε(x), so eigenfunctions for different k are orthogonal under ∫ε·X_m·X_n dx. This is L06 p. 29.

For the n = 1 mode in the top insulator (ξ = x + t_i + t_si/2 ∈ [0, t_i]):
- [LT04] eq. (9) gives u_L1 = [cos(πt_si/2λ₁)/sin(πt_i/λ₁)]·sin(πξ/λ₁);
- the printed g₁ (foot of LT04 p. 1388) is g₁ = [sin(πt_si/2λ₁)/cos(πt_i/λ₁)]·sin(πξ/λ₁).

Hence g₁/u_L1 = tan(πt_si/2λ₁)·tan(πt_i/λ₁) = ε_i/ε_si, by eq. (8) with n = 1:

ε_si·tan(πt_i/λ₁) = ε_i·cot(πt_si/2λ₁).

In silicon g₁ = u_L1 = cos(πx/λ₁). Therefore **g₁ = ε(x)·u_L1/ε_si**, and LT04 eq. (11) is the ε-weighted projection of L06 eq. (2.23). **[K1: residual ≤ 9×10⁻¹⁵.]**

## N2. The closed form b₁ (LT04 eq. 12) follows from eq. (11)

The boundary data is ψ(x, 0) − υ = A₀ in silicon, and rises linearly from 0 at the gate to A₀ at the interface in each insulator. With λ = λ₁:

**Numerator.**
- Silicon: A₀·∫cos(πx/λ)dx = A₀·(2λ/π)·sin(πt_si/2λ).
- Both insulators: 2·A₀·[sin(πt_si/2λ)/cos(πt_i/λ)]·(1/t_i)·∫₀^{t_i} ξ·sin(πξ/λ) dξ = 2A₀·sin(πt_si/2λ)·[(λ/π)²·tan(πt_i/λ)/t_i − λ/π].
- Sum: 2A₀·(λ/π)²·sin(πt_si/2λ)·tan(πt_i/λ)/t_i.

**Denominator.**
- Silicon: ∫cos² = t_si/2 + (λ/2π)·sin(πt_si/λ).
- Insulators: [sin(πt_si/λ)/sin(2πt_i/λ)]·t_i − (λ/2π)·sin(πt_si/λ).
- Sum: t_si/2 + t_i·sin(πt_si/λ)/sin(2πt_i/λ).

The ratio is LT04 eq. (12) (= L06 eq. 2.25) at V_ds = 0. **[K2: ≤ 10⁻¹⁴ relative, for 5 geometries.]**

## N3. The long-channel end correction implied by eqs. (13) + (16)

**Definition** (C-10): G = G_1D·L/(L − 2ΔL), with G_1D = qμW·N_∞/L and N_∞ = t_si·n_i·e^{(V_g−Δφ)/V_T}.

At V_ds → 0, LT04 eq. (16) gives 1/G = (qμW)⁻¹·∫₀ᴸ dy/N(y). With the lowest mode (13):

N(y)/N_∞ = M(β·h(y)),  M(z) = (1/2s)·∫_{−s}^{s} e^{z·cos θ} dθ,  β = b₁/V_T,
h(y) = [sinh k(L−y) + sinh ky]/sinh kL.

Hence 2ΔL(L) = ∫₀ᴸ g(β·h(y)) dy, with g(z) = 1 − 1/M(z). For L → ∞, substituting z = β·e^{−ky}:

**ΔL = (λ₁/π)·F(β, s),  F(β, s) = ∫₀^β g(z) dz/z.**

**Properties [D]:**
- (a) **Thin-film limit.** As s → 0, M(z) → e^z and F → ∫₀^β (1 − e^{−z}) dz/z = Ein(β). This is the C-10 form, and the C-11 DG predictor. **[K3: 5×10⁻⁸ at s = 10⁻³.]**
- (b) **Ordering.** cos θ ≤ 1 implies M(z) ≤ e^z, so g(z) ≤ 1 − e^{−z} and F(β, s) ≤ Ein(β), with equality only at s = 0. The pure-Ein predictor therefore always overestimates LT's ΔL: by +2.1 % to +7.2 % for t_si = 5–25 nm [N].
- (c) **Asymptotics.** F = ln β + C(s) + O(∫_β^∞ dz/(z·M)), with

  C(s) = ∫₀¹ g dz/z − ∫₁^∞ dz/(z·M(z)).

  C(0) = Ein(1) − E₁(1) = γ. C(s) decreases with s: C(0.5) = 0.536, C(0.8) = 0.473, C(1.2) = 0.351 [N].

  By Laplace's method, M(z) ≈ e^z·(2s)⁻¹·(2π/z)^{1/2} for s√z ≫ 1, so the remainder is O(β^{1/2}·e^{−β}).
- (d) **Sensitivity.** d ln ΔL / d ln b₁ = g(β)/F ≈ 1/F ≈ 0.27–0.33 for barriers of 0.3–0.5 V [K7]. **Errors in the boundary amplitude enter ΔL only logarithmically.**
- (e) **Finite L.** ΔL(L)/ΔL(∞) = 0.978–0.991 at L = 3λ₁, and 1 to within 10⁻⁴ for L ≥ 5λ₁ [K4].

## N4. Two-sided bounds on the subthreshold current at any V_ds (frozen potential)

**Slotboom form.** With Φ = e^{−φ_n/V_T}, the electron current is J_n = σ_S·∇Φ, where σ_S = qμ·n_i·V_T·e^{ψ/V_T}. Continuity without recombination gives ∇·(σ_S∇Φ) = 0.

**Boundary conditions:**
- Φ = 1 on the source face;
- Φ = e^{−V_ds/V_T} on the drain face;
- zero normal flux at the insulator interfaces.

**Linear problem.** In subthreshold the mobile charge does not affect ψ, so ψ is fixed and the problem is linear. Then

I = G_S·(1 − e^{−V_ds/V_T}),

with G_S the conductance of the σ_S medium between the two faces.

**Upper bound (Dirichlet principle).** Use the trial function Φ = Φ(y), i.e. φ_n constant across the film. Then G_S ≤ W·[∫dy (∫σ_S dx)⁻¹]⁻¹, so

**I ≤ I_slice = μW·kT·(1 − e^{−qV_ds/kT}) / ∫dy (∫n_i e^{qψ/kT} dx)⁻¹,  which is L06 eq. (2.30) = LT04 eq. (16).**

**Lower bound (Thomson principle).** Use the trial current J = J(x)·ŷ, i.e. independent rows. Then G_S ≥ W·∫dx (∫dy/σ_S)⁻¹, so I ≥ I_rows.

**Long channel.** ΔL_rows = (λ₁/π)·⟨Ein(β·cos θ)⟩_θ ≤ ΔL ≤ ΔL_slice = (λ₁/π)·F(β, s). The two bounds coincide in the thin-film limit.

**Numbers [K6, K8]:**
- end-correction bracket: 0.14–1.85 % (t_si = 5–25 nm);
- total current: I_rows/I_slice ≥ 0.9948 for t_si = 10 nm and L = 20–1000 nm at V_ds = 0.05 V and 1 V.

**Scope.** These are rigorous bounds on the transport approximation for a given potential. They say nothing about the error in the potential itself.

## N5. Why higher-order modes do not matter for ΔL

1. **Where ΔL is decided.** The integrand g changes from ≈ 1 to ≈ 0 near y*, where b₁·e^{−ky*} ≈ V_T, i.e. e^{−ky*} = β⁻¹. For y < y*, g = 1 − O(e^{−β·e^{−ky}}) whatever the higher modes do.
2. **Size of the next mode there.** Mode n contributes b_n·e^{−k_n·y*} = b_n·β^{−λ₁/λ_n}. For the symmetric family the next mode is n = 3, with λ₁/λ₃ ≈ 3, so its contribution relative to V_T is ≈ (b₃/b₁)·β⁻².
3. **Coefficient decay.** The b_n fall off as O(n⁻²), because the boundary data is continuous with a kink at the interfaces.

**[K5]** Full series (300 symmetric modes) against one mode:
- ≤ 0.02 % for SiO₂ insulators;
- 0.18–0.36 % for the thick high-k case t_i/t_si = 0.54, which is outside LT's interpolation condition anyway.

## N6. Long-channel V_t roll-off: an end term plus a second-order remainder

**Setting.** Let z_s = β_s·e^{−ky} and z_d = β_d·e^{−k(L−y)}, with β_s = b₁/V_T and β_d = c₁/V_T. The exact centre amplitude is z_s + z_d + O((β_s + β_d)·e^{−kL}).

**Decomposition:**

∫₀ᴸ g(z_s + z_d) dy = ∫₀ᴸ g(z_s) dy + ∫₀ᴸ g(z_d) dy + ∫₀ᴸ [g(z_s + z_d) − g(z_s) − g(z_d)] dy.

**Bounding the terms.**
- *Last term.* Since g(0) = 0, the second-difference bound gives |g(a + b) − g(a) − g(b)| ≤ ab·sup|g″|. Here z_s·z_d = β_s·β_d·e^{−kL} uniformly in y, so the last integral is ≤ L·sup|g″|·β_s·β_d·e^{−kL}. (g″ is bounded on [0, ∞).)
- *First two terms.* ∫₀ᴸ g(z_s) dy = ΔL_s − ∫_L^∞ g(z_s) dy, and the tail is ≤ ⟨cos θ⟩·β_s·e^{−kL}/k. The drain term is analogous.

**Result:**

∫₀ᴸ dy·N_∞/N(y) = L − ΔL_s − ΔL_d + R,  R = O((L + 1/k)·β_s·β_d·e^{−πL/λ₁}).

Note that β_s·β_d·e^{−kL} = (Δψ_min/2V_T)², with Δψ_min = 2√(b₁c₁)·e^{−kL/2} the central barrier lowering.

**V_t roll-off.** The constant-current roll-off of L06 eq. (3.10) is

ΔV_t = −V_T·ln[L/(L − ΔL_s − ΔL_d)] + O(V_T·(Δψ_min/2V_T)²).

- The first-order barrier lowering is part of the overlapping tails. Adding L06 eq. (3.7) on top therefore double-counts it.
- **[K8]** The remainder is +0.034 mV at 70 nm and +1.64 mV at 50 nm (V_ds = 50 mV; λ₁ = 17.7 nm), scaling as (Δψ_min/2V_T)².
- The end term alone reproduces the full integral to < 0.001 mV at 100 nm.

**The Taylor formula (3.15) vanishes faster than any power of 1/L.** D₁ ∝ e^{−πL/4λ₁}, so D₁L → 0. Since erf(a) = (2/√π)·a + O(a³), the logarithm in (3.15) → ln[1 + O(D₁²L²)] → 0. The minimum-potential term (3.7) is also exponentially small.

**Crossover.** The regimes cross where Δψ_min = V_T, at L_× = (2λ₁/π)·ln(2√(b₁c₁)/V_T). Because ΔL_s + ΔL_d = (2λ₁/π)·[ln √(β_s·β_d) + C(s)], this is L_× = ΔL_s + ΔL_d + (2λ₁/π)·(ln 2 − C(s)): the end regions meet.

## N7. The gate-all-around analogue (sketch, not computed)

For a cylinder of radius R, the lowest mode is J₀(r/λ_c), where λ_c solves [Oh00] eq. (6) (≈ (D_Si + 2κt_ox)/4.810 for thin oxide).

The same steps give ΔL = λ_c·∫₀^β [1 − 1/M_c(z)] dz/z, with

M_c(z) = (2/R²)·∫₀^R exp(z·J₀(r/λ_c)/J₀(0))·r dr,

and β the projected amplitude. The thin-wire limit is again Ein(β).

The 2-D analytical nanowire solution with this structure is Yu, Yuan, Song & Taur (IEEE T-ED 56(10), 2009) [E3]. Not checked numerically here.

## N8. Order of magnitude of gate-induced depletion at the source/drain edge (U1)

With N_D = 10²⁰ cm⁻³ and a gate-to-edge difference ΔV ≈ 0.4 V:

W ≈ (2ε_si·ΔV/(q·N_D))^{1/2} = (2 × 1.04×10⁻¹² × 0.4 / (1.6×10⁻¹⁹ × 10²⁰))^{1/2} cm ≈ 2.3 nm.

This is comparable with t_si/2 for t_si ≤ 5 nm, consistent with L06 §4.2 (the edge error grows as t_si/t_i decreases).

**[K9]** A rough estimate from L06 Fig. 4.2(c), values read by eye: b₁ is overestimated by ≈ 10 %, which changes the long-channel ΔL by ≈ −3 %.

**Note:** degeneracy and Thomas–Fermi screening, which would reduce W, are ignored here.
