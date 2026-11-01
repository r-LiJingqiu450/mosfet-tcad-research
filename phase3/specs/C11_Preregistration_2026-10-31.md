# C-11 Pre-registration — can the subthreshold end correction be predicted from 1-D quantities?

**Date:** 2026-10-31. **Branch:** `phase3/C11-connection-coefficient-2026-10-31`, from C-10 `da913b5`.

**Status:** written and committed before the predictors are evaluated on the test set.

**Test set:** the 24 deep-subthreshold states (α∞ ≥ 100) in `results/phase3/C10/analysis/c10_result.json`, with ΔL from the source-side current.
- These are existing 2-D solutions. The predictors below have **no free parameter**, and their outputs on these states have not been computed.
- **Disclosure:** while debugging `phase3/c11/bpredict.py`, H-lin was run on one state outside the test set (device 2, α∞ = 50). There it gave λ_lin = 12.26 nm (2-D asymptotic local decay length ≈ 12.4 nm) and B_lin/V_T = 16.8 against the C-10 window fit B_fit/V_T = 29.6. ΔL_pred was 41.7 nm against 50.0 nm (−17 %).
- That state also showed that the local decay length of the 2-D surface tail falls from ≈ 30 nm (x ≈ 6 nm) through 15 nm (u ≈ V_T) to ≈ 12.4 nm (u < 1e-4 V). The conductance-relevant region u ≈ V_T is therefore pre-asymptotic.

## Predictors

| Model | Definition | Code |
|---|---|---|
| **H-lin** (linearised eigen-projection) | λ and B from the 1-D transverse eigenproblem and the ε-weighted projection of the junction-face data; ΔL_pred = λ·Ein(B/V_T) | `phase3/c11/bpredict.py` (docstring gives the full definition) |
| **H-q2D** (classical nonlinear quasi-2-D surface-potential equation, as model M7 step 3 of C-7) | ε_Si W_eff ψ_s'' = Q_s(ψ_s) − C_ox(ψ_g − ψ_s), with W_eff = l₀²(C_ox + C_d)/ε_Si (linear decay rate l₀ by construction), ψ_s(0) = ψ_bi (+ shift), ψ_s(∞) = ψ_s∞. Then ΔL_pred = ∫₀^∞ [1 − e^{−u/V_T}(1 − r u/V_T)] dx, with u = ψ_s − ψ_s∞ and r = C_ox V_T/Q_dep,∞ | `phase3/c11/c11_eval.py` |

## Decision rules

Each predictor is applied to every state, with a two-sided tolerance of ±5 %.

| Outcome | Condition |
|---|---|
| **SUPPORTED** | ≥ 22 of 24 states within ±5 % and none outside ±10 % |
| **REJECTED** | median |error| > 5 % |
| **INCONCLUSIVE** | otherwise |

**Expectation**, stated before evaluation: H-lin is REJECTED (systematically low B). H-q2D is unknown; it was refuted near threshold in C-7 but never tested in deep subthreshold.

**Descriptive (no decision rule):** the asymptotic amplitude B_∞ of each 2-D tail, from a fit in its far tail (1e-6 V < u < 1e-4 V). The ratios B_∞/B_lin and B_fit/B_∞ are reported.

**Only if a predictor is SUPPORTED:** its predictions for new, unseen devices will be pre-registered in an addendum before those states are run.
