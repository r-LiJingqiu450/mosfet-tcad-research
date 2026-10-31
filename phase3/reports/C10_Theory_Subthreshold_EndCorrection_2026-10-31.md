# C-10 Theory Note — The Zero-Bias Subthreshold End Correction: Current-Integral Identity, Similarity, and Far-Field Asymptotics

**Date:** 2026-10-31. **Branch:** `phase3/C10-subthreshold-theory-2026-10-31`.

**Labels:**
- **[P]** proved (under the stated assumptions);
- **[C]** conditional identity (exact given an explicitly stated approximation);
- **[V]** verified numerically within the F0 model;
- **[O]** observation;
- **[I]** interpretation;
- **[E1]** full text read, **[E2]** abstract only, **[E3]** title/metadata only, **[T]** standard textbook result not re-read in this phase.

**Model (F0):**
- 2-D, zero drain bias, abrupt n⁺ wells;
- Boltzmann statistics; holes in equilibrium with the body; constant mobility μ;
- ψ measured from the intrinsic level, with ψ_b0 = −V_T ln(N_A/n_i) in the neutral p-bulk;
- ψ_bi = V_T ln(N_D/n_i) on the well, so the junction step is ψ̂_bi·V_T = V_T ln(N_D N_A/n_i²).

## 1. The end correction is the end correction of the subthreshold current integral

**Definitions.**
- The zero-bias channel conductance per unit width is G(L).
- The long-channel sheet conductance is σ∞ = qμ∫n∞(y)dy, where n∞ is the far-field (one-column) electron profile.
- The end correction is defined operationally by G = σ∞/(L − 2ΔL) for L ≫ the end-layer size. C-6/C-9 showed that ΔL is L-independent there [V].

**(1a) [P] Equilibrium-weighted network (C-6 P4).**
- At zero bias the linear response of the drift–diffusion system is that of a conductor with local conductivity s(x, y) = qμn₀(x, y), where n₀ is the equilibrium density.
- The Scharfetter–Gummel flux has ∂J/∂ψ = 0 at equilibrium. The discrete conductance is therefore that of a resistor network with equilibrium edge weights.

**(1b) [P] Rayleigh bounds (C-6 P5).**
- Shorting all vertical lines (slices) gives an upper bound G ≤ G_slice. Here 1/G_slice = ∫dx/σ(x), with σ(x) = ∫s(x, y)dy.
- Cutting all vertical links (rows) gives a lower bound G ≥ G_rows. Here G_rows = ∫dy[∫dx/s(x, y)]⁻¹.
- Hence ΔL_rows ≤ ΔL ≤ ΔL_slice, with

  ΔL_slice = ∫₀^∞ [1 − σ∞/σ(x)] dx    (per end, long-channel limit).

**(1c) [P] Prefactor-free form in subthreshold.**
- With n₀ = (n_i²/N_A)·e^{(ψ−ψ_b0)/V_T}, and electrons too few to affect Poisson's equation:

  σ(x)/σ∞ = ∫e^{ψ(x,y)/V_T}dy / ∫e^{ψ∞(y)/V_T}dy

  and ΔL depends only on the electrostatic potential.
- The row form, ∫₀^L e^{−ψ(x,y)/V_T}dx at fixed y, is the denominator of the classical subthreshold current integral. This is Fjeldly & Shur 1993, eq. (2) [E1]. Its long-channel end correction is ΔL_rows.

**(1d) [V] The bracket is tight in subthreshold.** In the subthreshold states (α∞ ≥ 15; the C-8 states and the 22 C-10 states), ΔL_slice − ΔL_rows ≤ 1.2 % of ΔL, except 5.4 % at ψ̂ = 12 for device 2. ΔL lies between them, computed from the source-side current (E10-1, C-10 report §5).

**Consequence.** At V_DS → 0 the operational conductance end correction and the end correction of the classical subthreshold current integral are the **same object** to within the bracket width (typically ≤ 1.2 %). Fjeldly–Shur's "L − x_s − x_d" (their eq. 8) is an approximation to exactly this quantity. Their x_s is the 1-D abrupt-junction depletion width (their eq. 9) evaluated at the interface potential step ψ_bi − ψ_s∞.

## 2. Similarity theorem

**[P] Statement.** Make the following assumptions:
- F0 physics with Dirichlet wells (deck W);
- electrons negligible in Poisson's equation (subthreshold, α∞ ≫ 1).

Then lengths in units of L_D = (ε_Si V_T/(qN_A))^{1/2} and potentials in units of V_T give a dimensionless boundary-value problem with the parameters:
- ψ̂_g = V_G/V_T (the decks' V_FB = 0 convention), equivalently ψ̂ = (ψ_s∞ − ψ_b0)/V_T;
- ψ̂_bi;
- t̂ = (ε_Si/ε_ox)t_ox/L_D;
- ĥ = h_w/L_D;
- the domain extents in units of L_D.

By (1c) the conductance ratio carries no other parameter. Hence

ΔL/L_D = F(ψ̂, ψ̂_bi, t̂, ĥ, …).

**Proof sketch.**
- Poisson's equation becomes ∇̂²ψ̂ = e^{−(ψ̂−ψ̂_b0)} − 1 + (electron term) in the silicon, and ∇̂²ψ̂ = 0 in the oxide. The oxide thickness enters only through the continuity of εψ_y, i.e. through t̂.
- The Dirichlet data are ψ̂_bi on the well face (0 ≤ ŷ ≤ ĥ), ψ̂_g on the gate and ψ̂_b0 at the bottom.
- The electron term is (n_i/N_A)²e^{ψ̂−ψ̂_b0}. It is the only place n_i/N_A enters, and it is dropped in the electron-free limit.
- The conductance ratio follows from (1c).

**Correction after P10-1 (post hoc).** The electron-free assumption fails next to the junction face.
- In deck W the face carries ψ = ψ_bi = V_T ln(N_D/n_i), so the electron density there is n_face = N_D.
- Relative to N_A this is n_face/N_A = N_D/N_A = e^{ψ̂_bi}(n_i/N_A)², which is 1e4 for device 1 and 1e2 for device 2.
- Those electrons screen part of the junction step. **n_face/N_A is therefore an additional dimensionless parameter**, and the corrected statement is ΔL/L_D = F(ψ̂, ψ̂_bi, t̂, ĥ, N_D/N_A, …).
- Device 2' of P10-1 matched t̂, ĥ and ψ̂_bi but had n_face/N_A = 1. It differed from device 1 by +7.5 % and +6.3 %, and P10-1 was **REJECTED**.
- An exact image that also matches n_face/N_A (n_i raised to 1e12 cm⁻³; `phase3/c10/probe10ni.py`; a mathematical test, not a device) reproduces device 1 within **+0.8 % and +0.6 %**. This is [V] for the corrected statement.

**Consequence for C-8/C-9.**
- Device 1 (1e16 cm⁻³, 20 nm) and device 2 (1e18, 2 nm) have **identical t̂ = 1.4675**: the oxide was scaled with L_D.
- They differ only in ψ̂_bi (36.84 vs 41.45) and in ĥ (12.3 vs 123).
- C-8 also compared them at equal α∞, not at equal ψ̂. Because α∞ contains (n_i/N_A)², equal α∞ means ψ̂ = 21.0 vs 30.5.
- **Approximate equality of ΔL/W_d for these two devices is therefore expected from similarity. It is not evidence of a universal constant.**
- **Test:** P10-1 checks the theorem with an exact similarity image of device 1. P10-2 measures the separate effects of ψ̂_bi and ĥ.

## 3. Far-field asymptotics

**[C] Assumptions.** Far from the junction (x ≫ near-junction region) the perturbation u = ψ − ψ∞ is small and obeys the linearised equation. Its slowest-decaying mode dominates:

u(x, y) ≈ B·Φ₀(y)·e^{−x/λ},  Φ₀(0) = 1.

λ is set by the depletion approximation:
- the two-layer eigenproblem tan(kW_d) = −(ε_Si/ε_ox)tan(k t_ox), λ = l₀ = 1/k;
- this is the "generalised scale length" of Frank, Taur & Wong (IEEE EDL 1998), λ_FTW = π l₀ [T];
- with Boltzmann holes the effective depletion edge is soft, and λ exceeds l₀ by a few per cent [O].

**The surface field factor.**
- The electron layer is thin (≈ V_T/E_s ≪ W_d), so σ(x)/σ∞ = e^{u_s/V_T}·E∞/E_s(x).
- Continuity of the normal displacement at the interface gives ε_Si u_y(0⁺) = C_ox u_s. Hence E_s/E∞ = 1 − r·u_s/V_T, with

  r = C_ox V_T/Q_dep,∞.

**Integral.** Substituting z = (B/V_T)e^{−x/λ}:

ΔL = λ·∫₀^{B/V_T} [1 − e^{−z}(1 − r z)] dz/z + R_near = λ[Ein(B/V_T) + r(1 − e^{−B/V_T})] + R_near,

where Ein(z) = γ + ln z + E₁(z) and γ = 0.5772.

**About R_near.**
- R_near is the error of replacing the true near-junction profile by the extrapolated tail.
- Where u ≫ V_T both integrands are ≈ 1. So R_near comes only from higher modes in the transition region u ≈ V_T, and is small.

**For B ≫ V_T:**

ΔL ≈ λ[ln(B/V_T) + γ + r].    (3.1)

**[I] Interpretation of (3.1).** ΔL is the distance at which the far-field S/D potential tail, extrapolated to the junction, has fallen to V_T·e^{−(γ+r)} ≈ 0.51 V_T. Four consequences:
- The scale is the 2-D scale length λ, not a depletion width. t_ox enters through λ.
- The bias dependence is weak: λ shrinks as W_d shrinks, while ln(B/V_T) grows as ψ_s∞ falls. ΔL(ψ̂) therefore has a broad stationary region. The C-8 "saturation" over α∞ = 15–150 (only 60 mV of ψ_s) is this stationary region.
  - **[O, post hoc]** This explains E8-1. With ΔL nearly constant and l₀ shrinking 3–4 %, κ = ΔL/l₀ rises by ≈ 0.14. The ln(B/V_T) term adds ≈ 0.1–0.2. The prediction is Δκ ≈ +0.3; observed +0.33/+0.40. The C-7 "H-log" law, κ ≈ ln α∞ + 1.6 (Δκ = +2.3 per decade), had the wrong logarithm: ln α∞ is linear in ψ_s, while (3.1) is logarithmic in ψ_bi − ψ_s∞.
- The 1-D depletion width x_s (Fjeldly–Shur eq. 9) measures where the 1-D junction potential reaches the channel value. It ignores the gate-controlled tail beyond that point. (3.1) counts the tail out to ≈ 0.5 V_T, so ΔL > x_s is expected whenever λ·ln(B/V_T) is comparable to x_s.
- **B is a connection coefficient.** It is fixed by the full nonlinear near-junction problem: the junction step, the lateral depletion of the junction face below the channel depletion edge, and the oxide fringe.
  - Define β = B/(ψ_bi − ψ_s∞). No a-priori formula for β is derived in this phase.
  - Because ΔL depends on B only logarithmically, a 20 % error in β changes ΔL by only ≈ 0.18λ (≈ 5 %).

## 4. Status

| Item | Status |
|---|---|
| (1a)–(1c) | [P] |
| (1d) | [V] (C-8 states and P10) |
| Similarity theorem | [P] once n_face/N_A is included. The registered, electron-free form was **rejected** by P10-1 (+6–7.5 %); the corrected form is [V] to 0.6–0.8 % (post-hoc exact image) |
| (3.1) | [C]. P10-5, with the r term, was **rejected** at its 5 % tolerance: errors +1.2 to +5.7 %, median +2.7 %. [O, post hoc] Without the r term, ΔL = λ_fit·Ein(B/V_T) holds to −1.7 to +2.8 % (median +0.1 %) on all 24 states with α∞ ≥ 100. Equivalently, κ_fit − ln(B/V_T) = 0.44 to 0.65 (median 0.574 ≈ γ) |
| E8-1 explanation | [O, post hoc] |
| β | Empirical in this phase |
