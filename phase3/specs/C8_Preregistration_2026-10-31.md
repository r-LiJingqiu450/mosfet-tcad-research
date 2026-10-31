# C-8 Pre-registration: Theory Test of κ(α∞) and Extraction-Consequence Analysis

**Date:** 2026-10-31. Committed **before** any C-8 simulation.

## E8-1 — Competing hypotheses for the weak-inversion branch of κ = ΔL/l₀

**Development trend** (C-7, α∞ = 1.2–7.6): κ ≈ ln α∞ + 1.6.

**Physical argument for saturation (stated before testing).**
- In weak inversion the electrons are electrostatically passive except where N(x) is large near the junction.
- That screened region shrinks to zero as N∞ → 0.
- So, deep in subthreshold, N(x)/N∞ ≈ exp(δψ⁰_s(x)/V_T)·E⁰_s∞/E⁰_s(x), where δψ⁰ and E⁰ come from the **electron-free** end problem and are independent of N∞.
- Hence κ should approach an electron-free limit that depends on V_G only through the slowly varying junction step (ψ_bi − ψ_s∞)/V_T and W_d. **κ cannot keep growing like ln α∞.**

**Competing hypotheses:**
- **H-log:** κ = ln α∞ + 1.6 continues. It predicts κ(α∞ = 15, 50, 150) = 4.31, 5.51, 6.61.
- **H-sat:** κ(150) − κ(15) is small. The only growth is through the junction step: ln[(ψ_bi−ψ_s)/V_T] grows by 0.34 (1e16) and 0.55 (1e18) between α∞ = 15 and 150. Times a tail-length factor ~1.3, that gives ≈ 0.45 and 0.71.

**States** (never simulated), deck W, abrupt junction, level 0:
- 1e16 cm⁻³ / 20 nm, L = 2 µm, V_G = 0.9294 / 0.8532 / 0.7835 V;
- 1e18 cm⁻³ / 2 nm, L = 1 µm, V_G = 1.2218 / 1.1475 / 1.0800 V.

Both sets correspond to α∞ = 15, 50, 150 (`results/phase3/C8/subth/preregistration_data.json`).

**Rule.** Δκ = κ(α∞ = 150) − κ(α∞ = 15), evaluated per structure.

| Outcome | Condition |
|---|---|
| **H-sat SUPPORTED** | Δκ ≤ 0.9 for both structures |
| **H-log SUPPORTED** | Δκ ≥ 1.8 for both structures |
| **INCONCLUSIVE** | Anything else |

## E8-2 — Consequences for extraction methods (exact consequences of the C-6 affine law; predictions of sign and form)

**Data.** A fine V_G sweep of the scaled device (N_A = 1e18, t_ox = 2 nm, deck W, L = 1 µm, V_G = 1.20…2.00 V in steps of 0.04 V). By P6 the zero-bias conductance at any length is

G(L, V_G) = σ∞(V_G)/(L − 2ΔL(V_G)),

with no series resistance. The mobility is constant, so every extracted artefact is purely electrostatic.

**Predictions.**

**(a) Terada–Muta.** The intersection of the 1/G-vs-L lines at two gate voltages V₁ < V₂ gives:
- a spurious series resistance R_sd,TM = 2[ΔL(V₁) − ΔL(V₂)]/[σ∞(V₂) − σ∞(V₁)] **> 0**;
- an extracted ΔL_TM **> ΔL(V₁) > ΔL(V₂)**.

These are exact algebraic consequences, stated before computing them.

**(b) Linear-extrapolation (maximum-g_m) threshold.** It shows an **algebraic** roll-off ΔV_T(L) ∝ −1/L at lengths where the exponential quasi-2-D roll-off is negligible (L/l₀ ≥ 40).

**(c) Split-C–V mobility.**
- The measured channel charge (gate-modulated, ∫C_gc dV_G) contains an end excess per end ΔL_C ≥ 0.
- Inequality: ΔL_C ≥ ΔL whenever the gate-modulated charge excess dominates. Then μ_split is **underestimated** by about 2(ΔL_C − ΔL)/L.
- The sign is to be checked numerically. ΔL_C is defined from the gate-modulated charge, not the total column charge, because the junction spill-over charge is not gate-modulated.

No pass/fail threshold applies to E8-2. It is a computation of consequences. Only its stated signs are falsifiable.
