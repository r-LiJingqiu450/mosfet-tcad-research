# C-7 Realism Test — Registration

**Date:** 2026-10-31. Written and committed **before any realism state was run.** The only prior runs are the two regression smoke states in `results/phase3/C7/dev`, which reproduce C-6 to 1e-4 nm.

## Question

Does the near-threshold source/drain end correction ΔL (C-6: 370 nm per end at V_G = 1.0 V for an abrupt junction with no overlap) survive realistic junction and gate geometry, and field-dependent mobility?

## Matrix

All states are zero bias, at mesh level 0 unless stated.

**Deck R** (resolved n⁺ wells, N_D = 1e20, h_w = 0.5 µm, l_w = 0.3 µm), L = 2 µm, V_G ∈ {1.0, 1.25, 2.0} V. Six (ov, σ) pairs (nm):

| ov | σ | Description |
|---|---|---|
| 0 | 0 | Baseline |
| 10 | 0 | Overlap |
| 30 | 0 | Overlap |
| 0 | 10 | Lateral erfc grading |
| 0 | 30 | Lateral erfc grading |
| 30 | 30 | Combined |

That gives 18 states.

**Deck W** ψ_bi sensitivity: ±50 mV, V_G ∈ {1.0, 1.25}, L = 2 µm. That gives 4 states.

**Checks**, for (ov, σ) = (30, 30) at V_G = 1.0:
- level 1, L = 2;
- level 0, L = 4.

Every state also yields the Lombardi-mobility ΔL by post-processing (exact under C-6 P4).

## Quantity

ΔL per end relative to the **metallurgical** channel length (N_D = N_A at x = 0, L):
- constant mobility, 1COL reference;
- Lombardi mobility, 1COL-mobility reference.

## Decision rules (fixed now)

**S1 — survives.** At V_G = 1.0 V, all of the following hold:
- for the combined case (30, 30), ΔL ≥ 0.25·W_dep (76 nm) with constant μ **and** with Lombardi μ;
- ΔL(1.0 V)/ΔL(2.0 V) ≥ 2 for that case, i.e. strong gate dependence persists.

**S0 — does not survive.** At V_G = 1.0 V for any realistic case, ov ∈ {10, 30} with σ ∈ {10, 30}: ΔL < 0.1·W_dep (30 nm), or ΔL changes sign.

**Intermediate.** Otherwise.

**Interpretation limits:**
- ΔL here is measured from the metallurgical junction. With grading, part of the donor tail lies on the n side and is conductive. That shift is an ordinary "lateral diffusion" contribution, not the intrinsic electrostatic end effect. Both are reported. The **gate dependence** ΔL(1.0) − ΔL(2.0) is the quantity least affected by where L is defined.
- Fermi–Dirac statistics and band-gap narrowing are represented only as a ±50 mV built-in-potential sensitivity.
