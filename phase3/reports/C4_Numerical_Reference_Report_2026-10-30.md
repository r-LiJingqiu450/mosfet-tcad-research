# C-4 Numerical Reference Report — T1 v2, Corrected R-I(b), R-I(c)

**Date:** 2026-10-30/31. **Branch:** `phase3/C4-2026-10-30`, from `e4e4528` (C-3 head, verified).
**Pre-registration:** `phase3/specs/C4_T1v2_Preregistration_2026-10-30.md`, commit `409e521`. It was committed before any registered run. The registered code (`phase3/c4/t1v2.py`, `deck4.py`) is unchanged from that commit through the last run (checked with `git diff`).
**Data:**
- `results/phase3/C4/t1v2_{B,C,D,Dprime}/` — 64 runs, all exit 0.
- `results/phase3/C4/analysis_t1v2/` — table, summary, figure.
- Development runs made before registration are kept in `results/phase3/C4/dev_smoke/`.

All E values are sums of **estimates**, not bounds. Thresholds are unchanged from R04 §8.1.

---

## 1. Task A — preservation

| Item | Result |
|---|---|
| C-3 head | `e4e4528e05ee…` |
| C-3 archive | `SHA256SUMS` re-verified (13 files); bundle `git bundle verify` OK |
| Ancestry | c909bf1 → e4e4528 |
| `main` | `c64c463`, untouched |
| Golden / production S1–S3 | Untouched (verified again in the final regression, §6) |
| C-4 branch | Non-main. No push: no write access, not retried. Preserved by bundle, ZIP and manifests (archive files) |

**Historical v1 runs.** All 72 were reproduced bit-for-bit inside the v2 runs:
- the unscaled-SuperLU G equals the v1 G_L1 at 72/72 points;
- the AC Re Y equals the v1 Re Y at 72/72 points.

The v1 NON-RESOLVING record (R-I(a), L = 10 µm, level 3, V_G = 2.0 V) is kept unchanged.

## 2. Task B — T1 v2: equilibrated matrix-export solver and L2-native error estimate

### 2.1 Verdicts (primary rule, registered E budgets)

| Deck | Points | L1v2 PASS | L2v2 PASS | Conservative (E_wc) L2v2 | max \|r₀\| L2v2 | max E L2v2 | max E L1v2 |
|---|---|---|---|---|---|---|---|
| R-I(a) (C-3 deck) | 36 | 36 | 36 | 33 PASS, 3 class I | 9.1e-8 | 2.1e-7 | 6.0e-7 |
| R-I(b) C-3 deck (historical) | 36 | 36 | 36 | 33 PASS, 3 class I | 9.1e-8 | 2.2e-7 | 5.2e-7 |
| R-I(b) corrected | 36 | 36 | 36 | 33 PASS, 3 class I | 8.1e-8 | 2.2e-7 | 5.5e-7 |
| R-I(c) f = 0.3 | 36 | 36 | 36 | 33 PASS, 3 class I | 5.4e-8 | 1.7e-7 | 5.6e-7 |
| R-I(c) f = 0.1 (L = 1, 10 µm) | 24 | 24 | 24 | 21 PASS, 3 class I | 7.0e-8 | 1.8e-7 | 5.6e-7 |
| R-I(c) f = 1 (L = 1, 10 µm) | 24 | 24 | 24 | 21 PASS, 3 class I | 7.0e-8 | 1.7e-7 | 5.4e-7 |
| **Total** | **192** | **192** | **192** | **174 PASS, 18 class I** | 9.1e-8 | 2.2e-7 | 6.0e-7 |

- No FAIL on any path or variant.
- **Under the worst-case rounding term**, every L = 10 µm, level-3 point (3 per deck variant) is NON-RESOLVING: E_wc = 1.1–1.9e-6, all other points PASS. This conservative verdict is reported as registered. It is not primary.

### 2.2 The former NON-RESOLVING point (R-I(a), L = 10 µm, level 3, V_G = 2.0 V)

| Solution | r₀ | Equilibrated backward error η | Error estimate |
|---|---|---|---|
| x_U, unscaled SuperLU (= v1 L1) | −2.68e-6 | 9.4e-3 | (v1: E = 9.3e-6, class I) |
| **x_E, equilibrated SuperLU (L1 v2)** | −2.42e-7 | 1.0e-15 | E_L1v2 = 6.0e-7 (sec 3.1e-7, dwr 2.2e-7, eq 3.9e-8, round 2.4e-8) → **PASS** |
| x_M, UMFPACK on the exported system | −1.75e-8 | 4e-15 | |
| x_NS, DEVSIM Newton step | −1.75e-8 | 2e-14 | |
| **x_AC, DEVSIM AC (L2 v2)** | −1.75e-8 | 4e-15 | E_L2v2 = 8.4e-8 (eq 2.5e-8, sec 2.5e-8, dwr 9.3e-9, round 2.4e-8) → **PASS** |

**Reading:**
- v2 resolves the point by both independent paths.
- The cause identified in C-3 is confirmed and removed: the unscaled factorization was not backward stable (η = 9.4e-3).
- The remaining L1v2–UMFPACK difference (2.2e-7) is a forward error at this conditioning (κ_G = 1.2e10). It lies within E_L1v2.

### 2.3 Error-estimate validation and conditioning

| Check | Result |
|---|---|
| Backward error η, equilibrated system (B3v2 ≤ 1e-13) | x_E ≤ 1.6e-15, x_M and x_AC ≤ 4.2e-15, x_NS ≤ 2e-14 at all 192 points. **B3v2: 192/192.** x_U up to 8.8e-2 |
| Does the estimate predict a real error? Equilibrated SuperLU − UMFPACK (two backward-stable solvers) vs E_dwr + E_round of both | Within the estimate at **192/192** (max ratio 0.93, median 0.72). The predicted sign is correct at 192/192 |
| DWR prediction for the *unstable* x_U | Unreliable point by point: ratio of predicted to actual −57 to +39, median 0.8. Reported. The residual of a non-backward-stable solution is not evaluated accurately enough in double precision, consistent with the C-3 finding that fixed-precision refinement diverges |
| Output condition κ_G = \|λ\|ᵀ(\|A\|\|x\| + \|rhs\|)/\|G\| | 1.5e7 (L = 1, level 0) to 1.2e10 (L = 10, level 3). Grows ≈ L² and ≈ 2× per level. Dominated by the surface electron-continuity rows |
| L1v2 vs L2v2 consistency, \|Δr₀\| ≤ 3(E₁ + E₂) | 192/192 (max ratio 0.14) |
| AC vs UMFPACK solution | ≤ 3.2e-13 relative (same factorization family; ω² effect ≤ 3.8e-13 = E_freq) |
| Decoupling (R-I(a), equilibrated path) | **36/36** (v1: 34/36). Max \|ΔG\|/G = 1.1e-7 |
| **B5v2** (linear-solve part of E_L1v2 ≤ 1e-8) | **Not met at 96/192 points.** Met at every L = 1 µm point (max 3.1e-9) and at L = 3 µm level 0. Not met at L = 3 µm levels 1–3 (24 of 36) and at every L = 10 µm point (max 5.6e-7) |

**Why B5v2 fails, and what it means:**
- The failures are not solver instability: B3v2 holds everywhere.
- They are the forward-error consequence of κ_G ≈ 1e9–1e10. No backward-stable double-precision solve can deliver 1e-8 there.
- The 1e-8 figure is therefore unattainable for L ≥ 3 µm on the finer meshes. The T1 resolving power (E ≤ 1e-6) is still met.
- Two routes would improve resolving power there: extended-precision residuals, or a better-conditioned formulation (scaled quasi-Fermi variables). Neither is needed for E ≤ 1e-6.

### 2.4 L2-native estimate

E_L2v2 uses only the AC solution and the UMFPACK adjoint:
- E_freq;
- E_eq;
- section current vs contact current of the AC solution;
- the dual-weighted residual of the AC system;
- the RSS rounding term.

It no longer contains the v1 cross term |Re Y − G_L1|. At the point where that term made v1 non-resolving, E_L2v2 = 8.4e-8 vs v1 E_L2 = 2.7e-6. **L1–L2 disagreement is now only a diagnostic**, and it is satisfied everywhere (§2.3).

## 3. Task C — R-I(b) bottom corners

### 3.1 What the equations were

From the exported matrix at L = 1 µm, level 0 (`bottom_rows` in every T1 v2 record):

| Deck | Corner ψ row (nnz) | Corner φ_n row (nnz) | Meaning |
|---|---|---|---|
| R-I(b) C-3 | 4 (bulk Poisson row) | 2 (drain: 3) | Corner belongs to the lateral contact, which has no ψ equation. ψ is a natural-BC unknown, not ψ_b0 |
| R-I(b) corrected | **1 (ψ = ψ_b0)** | 2 / 3 (φ_n = bias) | One-node contacts `bs`/`bd` carry **both** conditions R04 §4.2 assigns to that node: bottom ψ Dirichlet, and lateral φ_n Dirichlet "over the full depth". No node is shared by two contacts, which DEVSIM cannot represent |

The neighbouring bottom node has a ψ row of nnz 1 (body) and an electron row of nnz 8 (zero-flux bulk row). This is as intended in both decks.

### 3.2 Verification

| Check | C-3 deck | Corrected |
|---|---|---|
| Lateral ψ spread at V_D = 0 (max) | 1.7e-11 V (at the bottom corners) | **1.1e-16 V** |
| Bottom ψ − ψ_b0 (max) | 1.7e-11 V | **0** |
| Gauss closure | 4.3e-14 | 2.1e-15 |
| Continuity jump | 0 | 0 |
| Single-entry oxide rows | 0 | 0 (no unintended constraint) |
| Contact equations | source/drain: φ_n only | source/drain: φ_n only; bs/bd: ψ + φ_n; oxc: none |
| T1 v2 (L1v2/L2v2), 36 points | 36/36 PASS | 36/36 PASS |
| \|G_corrected − G_C-3\|/G vs combined linear-solve estimate | — | ≤ **0.59** of the estimate at every point: **not resolved** (max raw difference 3.1e-8, at L = 10 µm) |
| Refinement | spread 5.0 → 1.7e-12 V (L = 1 µm, levels 0→3) | 1e-16 V at all levels |
| Finite bias, 50 mV, L = 1 µm, level 0, V_G = 1.25 V | I = 0.014588270572574346 A/cm | 0.014588270572574344 A/cm (relative change 1.4e-16) |

**Does the correction change the intended physical model?**
- No. It makes the deck implement the stated R-I(b) specification at two nodes where the C-3 deck did not.
- The effect on every reported quantity is below resolution.
- The C-3 estimate (≤ 7.5e-12) stands. The C-4 T1 and diagnostic results use the corrected deck.

## 4. Task D — R-I(c) T1 and deck verification

The specification is in `phase3/specs/R-Ic_Reservoir_Strip_Model_Specification_2026-10-30.md`.

| Check | Result (f = 0.1, 0.3, 1) |
|---|---|
| Contact equations | source/drain strips: φ_n only (ψ natural); bs/bd: ψ_b0 + φ_n; body: ψ_b0; gate: ψ; oxc: none |
| Strip bottom rows (exported matrix) | ψ row nnz 1, φ_n row nnz 2–3 at every strip bottom node examined; channel bottom nodes ψ nnz 1 / electrons nnz 8 |
| Unknowns, level 3 | 246,116 (f = 0.3, any L) vs 176,508 without strips |
| Continuity / Gauss / single-entry rows | 0 V / ≤ 2.9e-14 / 0 |
| Lateral uniformity at V_D = 0 (incl. strips) | ≤ 1.1e-16 V |
| T0/T2 compatibility | The 2-D V_D = 0 state equals the 1COL solution: uniform to 1e-16 V; T1 identity at E ≤ 2.2e-7. With C-2's T0/T2 (1COL vs M0, asymptotic, within 2·U_disc), R-I(c) inherits T0/T2 |
| T1 v2 | 84/84 PASS on both paths; conservative verdict class I only at L = 10 µm, level 3 |
| Units / current definition | Terminal current = −V1.I = drain + bd electron current (equal to 2e-15 relative); KCL ≤ 5e-15 (ℓ_res diagnostic) |

The ℓ_res dependence of the finite-bias discrepancy is reported separately (`C4_lres_and_TC1_Results`). It is a scientific question, not a deck check.

## 5. Status of the reference decks after C-4

| Deck | T1 (E ≤ 1e-6) | Boundary verification | Open items |
|---|---|---|---|
| R-I(a) | 36/36 both paths (v2) | as C-3 | none (the v1 class I is explained and resolved in v2) |
| R-I(b) corrected | 36/36 | complete (§3) | none |
| R-I(c), 3 lengths | 84/84 | complete (§4) | ungated-strip variant not defined or tested |

**Limits that remain:**
- The conservative (worst-case rounding) verdict is non-resolving at L = 10 µm, level 3, for every deck.
- B5v2 (≤ 1e-8) is unattainable at L ≥ 3 µm on fine meshes because of conditioning.
- T1 remains a zero-bias identity. Finite-bias correctness rests on conservation (δ_KCL, δ_cont ≤ 3e-14) and mesh refinement.

## 6. Regression

See the archive verification and `results/phase3/C4/regression_final/`.
