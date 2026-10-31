# C-4 Pre-registration — T1 version 2, corrected R-I(b), R-I(c), and the ℓ_res diagnostic

**Date:** 2026-10-30 (committed before any registered run). **Branch:** `phase3/C4-2026-10-30`.
**Code:** `phase3/c4/t1v2.py` (its module docstring is the binding algorithm text), `phase3/c4/deck4.py`, `phase3/c4/run_t1v2.sh`.

**What does not change:**
- T1 v1 (`phase3/c3/t1.py`, results `results/phase3/C3/t1/`) stays the historical record. Its NON-RESOLVING point (R-I(a), L = 10 µm, level 3, V_G = 2.0 V) stays classified as such.
- The v2 results are a **separate, versioned experiment** (`t1v2-2026-10-30`).
- **Thresholds are unchanged:**
  - PASS if |r₀| ≤ 3E and E ≤ 10⁻⁶;
  - NON-RESOLVING (class I) if E > 10⁻⁶;
  - FAIL otherwise.

## 1. Algorithm (Task B)

Per V_D = 0 state, one linear system A x = rhs, with rhs = (b(V1 = δ) − b(0))/δ. That rhs has a single nonzero, in the V1 circuit row (checked). Five solutions are computed:

| Label | Solver | Role |
|---|---|---|
| x_U | SciPy SuperLU, unscaled (exactly the v1 L1 expression) | comparison; reproduces v1 G_L1 |
| x_E | **Equilibrated SuperLU**: A_s = D_r A D_c (row max, then column max); y = A_s⁻¹(D_r rhs); x = D_c y | **L1 v2** |
| x_M | UMFPACK 5.1 on the exported (A, rhs) via `devsim.umfpack.umfshim` | reference/diagnostic |
| x_NS | DEVSIM's own undamped Newton step at V1 = δ, ÷ δ | diagnostic (= v1 G_NS) |
| x_AC | DEVSIM small-signal AC at f₁ = 1 Hz and f₂ = 0.1 Hz; system (A + jωA_d) x = rhs | **L2 v2** (primary) |

**Adjoint:** Aᵀλ = c with G = cᵀx. λ_M is computed by UMFPACK (transpose) and λ_E by the equilibrated factorization.

**Error components** (estimates, not bounds; definitions in the docstring):

| Component | Definition | Meaning |
|---|---|---|
| E_dwr | \|λᵀ r\|/G | First-order estimate of the solution error in G |
| E_round | u‖λ∘s‖₂/G, with s = \|A\|\|x\| + \|rhs\| | Probabilistic rounding uncertainty of that estimate (Higham & Mary 2019) |
| E_cond | κ_G u | Worst-case counterpart of E_round |
| E_sec | | Section-vs-contact current of the same solution |
| E_eq | | Change of r₀ under τ_eq ÷ 100 |
| E_freq | | (AC only) |

**Budgets:**
- E_L1v2 = E_eq + E_sec + E_dwr + E_round, using x_E and λ_E only.
- E_L2v2 = E_freq + E_eq + E_sec + E_dwr + E_round, using x_AC and λ_M only.
- **Neither budget uses the other path.** The v1 cross term |Re Y − G_L1| is dropped from L2.
- The conservative sensitivity verdict uses E_wc (E_cond in place of E_round). It is reported for every point but is not primary.

**Why E_round and not E_cond is primary** (development evidence, `results/phase3/C4/dev_smoke/`; produced before registration and kept):
- At R-I(a), L = 10 µm, level 2, V_G = 2.0 V:
  - the actual equilibrated-SuperLU − UMFPACK difference is −4.5e-8;
  - the DWR prediction is −3.4e-8;
  - E_round = 1.7e-8;
  - E_cond = 6.6e-7.
- At L = 1 µm, level 0: actual −4.45e-11, predicted −4.31e-11.
- κ_G is dominated by cancellation in the surface electron-continuity rows. The worst-case alignment of all row errors is not observed.

**Diagnostics** (not part of the verdict):
- the L1v2–L2v2 consistency check;
- B3v2 (η(x_E) ≤ 1e-13);
- B5v2 (linear-solve part of E_L1v2 ≤ 1e-8);
- DWR validation on x_U (predicted vs actual U − M);
- decoupling (R-I(a), equilibrated path);
- the bottom-row equations;
- continuity, Gauss closure, single-entry oxide rows, lateral uniformity.

## 2. Decks (Tasks C, D)

| Model | Deck | Notes |
|---|---|---|
| `a` | R-I(a), C-3 deck, same builder call | — |
| `b_c3` | R-I(b), C-3 deck (historical: bottom corners carry no ψ condition) | — |
| `b` | **R-I(b) corrected.** Lateral contacts end one row above the bottom. Each bottom corner is a one-node contact carrying ψ = ψ_b0 **and** φ_n = bias | Both conditions R04 §4.2 states for that node. No node is shared by two contacts |
| `c` | **R-I(c):** gated Si strips of length ℓ_res = f·W_dep, f ∈ {0.1, 0.3, 1}, W_dep = 0.30393 µm, at both ends | Same doping and physics. Every strip node above the bottom row has φ_n = bias, with ψ natural. The strip bottom row has ψ = ψ_b0 and φ_n = bias. Outer lateral boundaries are Neumann. Strip mesh uniform, N = round(10·√2^k) cells. Channel mesh identical to a/b. f → 0 reproduces `b` |

**R-I(c) interpretation recorded now:** R04 §4.3 does not say whether the gate covers the strips. Here it does, so the V_D = 0 state is laterally uniform and the zero-bias identity applies with the same channel length L. Ungated strips would be a structurally different (R-II-like) reference and are not tested.

## 3. Registered run matrix

| Set | Models | L (µm) | Levels | V_G (V) | Points |
|---|---|---|---|---|---|
| T1v2-B | a, b_c3 (with x_U) | 1, 3, 10 | 0–3 | 0.8, 1.25, 2.0 | 72 |
| T1v2-C | b | 1, 3, 10 | 0–3 | 0.8, 1.25, 2.0 | 36 |
| T1v2-D | c, f = 0.3 | 1, 3, 10 | 0–3 | 0.8, 1.25, 2.0 | 36 |
| T1v2-D′ | c, f ∈ {0.1, 1} | 1, 10 (3 µm if time allows; stated either way) | 0–3 | 0.8, 1.25, 2.0 | 48 (+24) |

## 4. ℓ_res diagnostic (Task D; not FB-1, no H-C decision)

**Scope:**
- L = 1 µm only, V_G = 1.25 V only, V_D = 50 mV only.
- Levels 0–3.
- Models a, b, and c with f ∈ {0.1, 0.3, 1}. f = 3 is added and labelled **exploratory, not pre-specified**.

**Quantities:**
- Same-mesh discrepancy ε_h^X = I_GCA,h/I_X,h − 1 (GL8 over 1COL slices, as in FB-1).
- Its R04 §7.4 classification (corrected order formula) and extrapolated value with U_disc.
- S_res(X, Y) for every pair.

**Pre-stated reading rules:**
- If ε^c varies with f beyond 2·U_disc, the GCA discrepancy is reservoir-specification dependent **within this family**. That says nothing outside the family.
- No average of models is called a reservoir-independent reduction error.
- Three or four tested reservoirs never establish independence.
