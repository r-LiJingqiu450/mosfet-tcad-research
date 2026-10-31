# C-9 Pre-registration: Independent Audit of the C-8 Claims and a Physical-Realism Test

**Date:** 2026-10-31. Committed before the C-9 runs.

**Disclosed exception:** two development smoke states in `results/phase3/C9/dev`, at 1e18 cm⁻³ / 2 nm, L = 1 µm, V_G = 1.2914 V:
- deck W, no halo: ΔL = +43.75 nm, which reproduces C-7;
- deck R with a 2e18 halo: ΔL = −706 nm.

**Device and matrix:** scaled device (N_A = 1e18, t_ox = 2 nm, N_D = 1e20, F0 physics); **direct 2-D solves** at each length (no affine extrapolation). V_G ∈ {1.08, 1.24, 1.32, 1.40, 1.48, 1.56, 1.64, 1.72, 1.80, 1.90, 2.00} V, with these decks:

| Deck | Lengths (µm) |
|---|---|
| W (Ward truncation) | 0.5, 1, 2 |
| R (resolved abrupt wells) | 0.5, 1, 2 |
| R + halo 2e18 (x_h = 20 nm, σ = 15 nm, depth 15 nm, σ_y = 15 nm) | 0.5, 1, 2 |
| R + halo 5e17 | 0.5, 1 |

## Hypotheses and rules

**H9-1, affine length law at short L.**
- Decks W and R without halo, V_G ≥ 1.24 V.
- **SUPPORTED** if |ΔL(0.5 µm) − ΔL(2 µm)| ≤ 1 nm at every V_G. Otherwise **REJECTED** in the range where it fails.

**H9-2, C-8 split-C–V proxy vs terminal charge.**
- At V_G = 1.32 V, L = 1 µm, decks W and R, constant μ.
- Compare the mobility error computed from the terminal charge Q_SD (S/D contact displacement charge + electrons, integrated from V_lo = 1.08 V) with the C-8 electron-excess proxy.
- **EQUIVALENT** if they differ by ≤ 20 % of the proxy error. Otherwise the proxy is **NOT EQUIVALENT**.

**H9-3, halo realism.**
- Halo 5e17, L = 1 µm, V_G = 1.32 V.
- **PREDICTED:** ΔL < 0, i.e. the halo-induced reverse end effect dominates the intrinsic +40 nm.
- The 2e18 case is already disclosed above (−706 nm at 1.2914 V) and is reported, not tested.

**H9-4, maximum-g_m threshold with field-dependent (Lombardi) mobility, deck W.**
- The g_m peak is now interior. ΔV_T,lin(L) = V_T(L) − V_T,GCA-mobility(∞).
- **PREDICTED:** ΔV_T(1 µm) ∈ [−15, −3] mV and ΔV_T(0.5)/ΔV_T(1) ∈ [1.7, 2.3], i.e. 1/L.
- **REJECTED** otherwise.
- Sensitivity to the upper sweep limit (1.80 vs 2.00 V) is reported.
