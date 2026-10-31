# Zero-Bias End Correction of the Linear-Regime MOSFET Conductance — Theory Note

**Date:** 2026-10-31. **Status:** the project's own derivations (level P). Each statement is marked as one of:
- **[proved]** a complete short argument is given;
- **[sketch]** an argument with a stated gap;
- **[hypothesis]**.

**Setting:**
- F0 physics: Boltzmann statistics, frozen holes, constant μ_n, no generation/recombination.
- Channel 0 < x < L. n⁺ wells as in Ward's truncation ([WOC90] (2.4a), (2.5)): ψ = ψ_bi and φ_n = 0 / V_D on the well edges.
- 2-D oxide; gate ψ = V_G + ψ_b0. Bottom contact ψ = ψ_b0, with zero electron flux.

## P4 [proved]: exact electrostatic decoupling of the zero-bias conductance

**Continuum.** J_n = −qμ_n n∇φ_n. At equilibrium ∇φ_n = 0, so the first variation is δJ_n = −qμ_n n₀∇δφ_n. Hence δφ_n solves ∇·(n₀∇δφ_n) = 0 with the contact data. The O(V_D) current is the conductance of the medium σ₀ = qμ_n n₀(x, y). It does not depend on δψ: electrostatic feedback first enters at O(V_D²).

**Scharfetter–Gummel discretization.** With Δ = (ψ₀ − ψ₁)/V_T and B(t) = t/(eᵗ − 1), the DEVSIM edge flux can be written as

J_e = qμ_n·EIL·V_T·n_i·B(Δ)·e^{ψ₀/V_T}·(e^{−φ₁/V_T} − e^{−φ₀/V_T}).

The bracket vanishes at equilibrium, so ∂J_e/∂ψ = 0 there, and δJ_e = qμ_n·EIL·n₀B(Δ)·(δφ₀ − δφ₁).

The discrete zero-bias conductance is therefore **exactly** that of a resistor network:
- each Voronoi edge has conductance c_e = qμ_n·EdgeCouple·EIL·n₀B(Δ);
- n₀B(Δ) = n₁B(−Δ) holds at equilibrium.

**Numerical confirmation.** An independent network assembly agrees with DEVSIM's full coupled-Jacobian linear response to 2e-9 relative. That is within the latter's adjoint error estimate of 1e-8 (`results/phase3/C6/w1`, `cross_check`).

**Consequences:**
1. The zero-bias GCA error is purely a property of the **equilibrium** electron distribution n₀(x, y), plus the 2-D geometry of the current flow in it.
2. The project's T1 identity (Prop. 1) is the special case n₀ = n₀(x), laterally uniform.
3. Ward's §6.4 replaces n₀ by an asymptotic construction and the flow by slices (6.33).

## P5 [proved]: two-sided variational bounds

Let G be the network conductance between the two well edges. Define:
- **G_slice**: every vertical edge conductance replaced by +∞ (each column becomes an equipotential);
- **G_rows**: every vertical edge conductance set to 0.

By Rayleigh's monotonicity law,

  **G_rows ≤ G ≤ G_slice**.

In the continuum, the Dirichlet principle with trial potentials φ = φ(x) gives the upper bound, and Thomson's principle with purely lateral trial currents gives the lower bound.

G_slice is Ward's quasi-1-D current hypothesis (6.33) evaluated with the **exact** 2-D charge:
- 1/G_slice = ∫ dx / (qμ_n N(x)), with N(x) = ∫ n₀ dy.

**Consequence:** with exact charge, Ward's slice treatment always **overestimates** the end enhancement of the conductance. The current-spreading correction ΔL_slice − ΔL ≥ 0 has a definite sign.

## P6 [sketch]: affine length law

Hypothesis of the sketch: the equilibrium n₀(x, y) approaches the laterally uniform 1-D profile exponentially fast away from each well, with a decay length ℓ_E. The supporting argument is that the linearized Poisson operator −∇·ε∇ + q(n₀ + p₀)/V_T is coercive in the strip with the gate and bottom Dirichlet conditions.

Under that hypothesis, the conduction problem separates for L ≫ ℓ_E into two semi-infinite end problems, each with a finite **end resistance**, plus exponentially small interaction terms. Therefore

  **1/G(L) = (L − 2ΔL)/σ_∞ + O(e^{−L/ℓ}),  with σ_∞ = qμ_n N_∞,**

and ΔL(V_G) is independent of L.

**Gaps:**
1. Coercivity gives exponential decay of the electrostatic perturbation in an energy norm. Transferring it to n₀ = n_i e^{ψ/V_T} pointwise needs an L∞ estimate, which is not proved here.
2. The transverse conduction spectrum for the weight n_∞(y) can have very small eigenvalues: the weakly coupled bulk-minority layer is an example. The corresponding slow modes carry conductance shares of order n_bulk/n_inv ~ 10⁻¹⁴ and are dynamically irrelevant, but they make the "exponentially small" remainder non-uniform.

**Contrast with Ward's ordering.** [WOC90] (3.1) leaves an algebraic O(ε²) interior correction in general. At V_D → 0 that correction vanishes identically, because the leading outer potential is y-independent. So the remainder after ΔL is beyond all orders in ε. This is testable: ΔL(L) must be constant for L ≫ ℓ.

## P7 [sketch]: equivalent threshold-shift form

Suppose the end enhancement is a local threshold lowering δV(x), so that qN(x) ≈ C_ox(V_G − V_T∞ + δV(x)). Then the slice expression gives

  ΔL ≈ ∫ δV/(V_G − V_T∞ + δV) dx  →  A/(V_G − V_T∞),  with A = ∫ δV dx at large overdrive.

The end correction is then indistinguishable from an **L-dependent linear-regime threshold shift** ΔV_T,lin = −2A/L, which is **algebraic in L**. The barrier-controlled (minimum-surface-potential) threshold roll-off of quasi-2-D theory is **exponential** in L (Liu et al. 1993, E3).

**[hypothesis]** For abrupt n⁺ junctions, the long-channel tail of the threshold extracted from linear-regime I–V is therefore dominated by the algebraic end term, not by the exponential roll-off.

Testable in the 2-D data through:
- ΔL·N_∞ ≈ const at large V_G;
- a comparison of ΔV_T,lin with the shift of the channel-centre surface potential.

**Test result (C-6 W1, L = 4 µm): the simple form of P7 is rejected.**
- ΔL·N_∞ rises monotonically through 2.3e12, 2.47e13, 3.66e13, 4.79e13 and 5.79e13 nm·cm⁻² for V_G = 1.0, 1.25, 1.5, 2.0 and 3.0 V. It does not saturate.
- The two-term form ΔL = ΔL_g + A/N_∞, fitted at 1.25 and 2.0 V, mispredicts 1.5 V by −13 % and 3.0 V by +39 %.
- The end correction is therefore neither a constant length nor a constant threshold shift.
- P6 (affine in L) still holds: at each V_G the correction is an L-independent length. So the algebraic 1/L character of the end term in linear-regime extraction stands. Its V_G dependence is not a single-parameter law.

## Observation (from the 2-D data, not derived): structure of the end layer

The 2-D surface-potential excess ψ_s(x) − ψ_s,∞ has two parts:
1. **A short layer** next to the n⁺ edge, with a 1/e length of 2–6 nm in inversion (21 nm at V_G = 1.0 V).
2. **A long tail.** Its local decay length falls from about 150–170 nm (fitted over 150–450 nm) to about 112–124 nm (fitted over 300–800 nm). Its extrapolated amplitude falls from about 130 mV at V_G = 1.0 V to about 1 mV at 3.0 V.

The tail length is comparable to the quasi-2-D natural length √(ε_Si t_ox W_dep/ε_ox) ≈ 135 nm (Young 1989; Liu et al. 1993, E3), so the tail is consistent with charge sharing.

Ward's assumed profile ([W88] p. 121) has only the oxide-fringe length 2t_ox/π = 12.7 nm and **no depletion-scale tail.** The inversion charge depends exponentially on ψ_s near threshold. This missing tail is the dominant error of [W88] §6.4 near threshold (C6 report, §1).
