# Finite-Drain-Bias Discrepancy Experiment — Design for Scientific Falsifiability (FB-1)

**Date:** 2026-10-30. **Status:** design only. **Not authorized for execution** (it is part of Package E, which requires independent review).
**Diagnostic basis:** a deliberately minimal set of states:
- `results/phase3/C3/fb/`: L ∈ {1, 3, 10} µm at level 0; L = 1 µm at levels 1–2; L = 10 µm at level 1; V_G ∈ {0.8, 1.25, 2.0} V; V_D = 50 mV;
- the L3 δ-ladders of the T1 runs (V_D = 0.25–8 mV, level 0).

All diagnostic numbers are previews without a full uncertainty analysis.

## 1. What the diagnostics already show (and why the design must change)

| Observation (preview) | Value | Consequence for the design |
|---|---|---|
| r(δ) = I_2D,h/I_GCA,h − 1 grows as δ² with no linear term | log-slope 2.000–2.025 for δ = 0.25–8 mV wherever resolved (V_G = 1.25, 2.0 V; both models; L = 1, 10 µm) | The natural small-V_D quantity is the **coefficient** c = lim ε/V_D², not ε at one V_D |
| The two reservoir models give discrepancies of **opposite sign and nearly equal size** at L = 1 µm | V_D = 50 mV, V_G = 1.25 V: ε_h^(a) = −5.30e-4, ε_h^(b) = +5.52e-4 | At L = 1 µm "the GCA reduction error" is not defined independently of the reference's reservoir model |
| Reservoir-antisymmetric part α = (ε^(a) − ε^(b))/2 decays ≈ L⁻² | 1.25 V: −5.41e-4 → −6.08e-6 (L = 1 → 10 µm); 2.0 V: −3.00e-5 → −3.19e-7 | End/reservoir effect with L⁻² scaling (ideal reservoirs, no access resistance) |
| Reservoir-symmetric part σ = (ε^(a) + ε^(b))/2 does **not** decay with L at fixed lateral cell count | 2.0 V: −3.3e-6 (L = 1) vs −4.2e-6 (L = 10); 1.25 V: +1.1e-5 vs −1.6e-5 | Either a physical reservoir-independent reduction term, or a lateral discretization error (dx/L is fixed at 1/50). **Decided by refinement, §1.1** |
| Same-mesh comparison removes most discretization error | V_G = 0.8 V: same-mesh ε_h ≈ 6e-10 vs continuum-M0 ε ≈ 1.2e-4 (L = 1 µm, level 0) | The same-mesh metric is the primary measurement; the continuum comparison is secondary |
| The discrepancy peaks near threshold | ε_h at V_D = 50 mV, L = 1 µm: 6e-10 (0.8 V), 5e-4 (1.25 V), 3e-5 (2.0 V) | Moderate inversion is the most informative bias regime; subthreshold is a negative control |

### 1.1 Is σ physical or numerical? (lateral/depth refinement, L = 1 and 10 µm)

FILLED_FROM_LEVEL_RUNS

## 2. Definitions (pre-registered for FB-1)

- **References.** X ∈ {a, b, c}. R-I(a): lateral ψ = 1COL data. R-I(b): lateral ψ Neumann. **R-I(c), new and required:** reservoir strips of length ℓ_res ∈ {0.1, 0.3, 1}·W_dep with uniform φ_n and natural ψ (R04 §4.3). All use the triple-point treatment of Research_08 §5.
- **Current discrepancy (primary, same mesh).** ε_h^X(V_D) := I_GCA,h(V_D)/I_X,h(V_D) − 1.
  - I_X,h is the 2-D drain current on level h.
  - I_GCA,h = (qμ_n/L)∫₀^{V_D} N_s,h^FV(φ)dφ from 1COL slices on the identical depth lines, by 8-point Gauss–Legendre; 6-point for the quadrature estimate.
- **Continuum discrepancy (secondary).** ε^X(V_D) := I_M0/I_X,1 − 1 with I_X,1 the finest level. It is reported with U_num of I_X (R04 §7.4) and is interpretable only where resolved.
- **Small-bias coefficient.** c_h^X := ε_h^X/V_D², estimated by a weighted fit over the δ-ladder V_D ∈ {1, 2, 4, 8} mV. The quadratic law is a test (H-A), not an assumption.
- **Reservoir decomposition.**
  - σ_h := (ε_h^a + ε_h^b)/2 (reservoir-symmetric);
  - α_h := (ε_h^a − ε_h^b)/2 (reservoir-antisymmetric);
  - S_res(a,b) := I_a/I_b − 1 = −2α + O(α²).
  - With R-I(c): σ_h^(X,Y) for every pair.
- **Numerical uncertainty (estimates, R04 §7.1).**
  - U_disc: R04 §7.4 applied to the ε_h (not I) sequence on 4 levels with depth and lateral refined together. Actual h_rep ratios; classification asymptotic / non-asymptotic / oscillatory.
  - U_iter: τ ÷ 100 repeat on ≥ 10 % of states.
  - U_quad: GL 8-point vs 6-point.
  - U_lin: δ_cont and the cross-solver difference.
  - U_num := U_disc + U_iter + U_quad + U_lin (reporting convention, not a bound).
- **Finite-current metric.** ε at V_D = 50 mV (linear region, V_D ≤ 0.5·V_Dsat checked per point).
- **Low-current metric.** δV_G(I′) = V_G,GCA(I′) − V_G,X(I′) in mV at fixed I′ on a ≤ 0.01 V V_G grid, with U_VG = U_rel,I·mV_T (R04 §6). Used for V_G ≤ V_T,cc.
- **Resolved.** abs(ε) > 2·U_num, and on the same deck and level: T1 passed, attained Newton updates ≤ 1e-12 V, δ_cont ≤ 1e-8.

## 3. Hypotheses that the experiment can falsify

| ID | Hypothesis | Falsified if |
|---|---|---|
| H-A | For both reservoir models ε_h ∝ V_D² at small V_D (no linear term) | The log-slope of abs(ε_h) vs V_D on {1, 2, 4, 8} mV differs from 2 by more than its propagated uncertainty, at any resolved point |
| H-B | α (the reservoir-antisymmetric part) scales as L⁻² | Two-point log-slopes between L ∈ {1, 2, 4, 10} µm differ from −2 beyond uncertainty |
| H-C | A reservoir-independent reduction discrepancy exists: σ is resolved, survives refinement, is equal for all three reservoir models within U_num, and decays with L | σ unresolved after refinement (**numerical**); or σ^(a,b) ≠ σ^(a,c) beyond U_num (**reservoir-dependent**); or σ does not decay with L |
| H-D | The discrepancy peaks in moderate inversion (V_G − V_T,cc ≈ 0.1–0.3 V) | The maximum of abs(c) over V_G ∈ {1.0, 1.1, 1.25, 1.5, 2.0} V is not in that range |

## 4. Minimal discriminating experiment (for authorization)

- Lengths: L ∈ {1, 2, 4, 10} µm.
- Gate: V_G ∈ {1.0, 1.25, 2.0} V, plus the control 0.8 V.
- Drain: V_D ∈ {2, 8, 50} mV.
- Reservoir models: a, b, c (ℓ_res = 0.3·W_dep only).
- Mesh: levels 0–3.

That is 4 × 4 × 3 × 3 × 4 = **576 bias states** (about 2–3 solver calls each), on DEVSIM plus the 1COL slices. Estimated wall time **6–9 h** on 2 cores, from measured level-3 timings (~3 min per V_G per mode at L = 1 µm).

**Cheaper first stage (stop/go):** L ∈ {1, 10} µm, V_G = 1.25 V, V_D ∈ {8, 50} mV, models a, b, c, levels 0–3. That is **48 states, about 1 h**. It decides H-C. If σ is numerical or reservoir-dependent, the full FB-1 is **not** run and D1 is reframed (§6).

## 5. Inclusion, exclusion, failure, inconclusive

- **Include:** points satisfying the resolved-rule prerequisites (T1 pass, Newton, δ_cont); abs(I_X) ≥ 10·U_abs; V_D ≤ 0.5·V_Dsat.
- **Exclude** from claims, but report: subthreshold points with abs(ε) below 10× the linear-solve noise floor (negative control only); any state with failed Newton or δ_cont > 1e-8 (class N; reported, never dropped).
- **Failure (class P):** T1 failure on any deck used; σ growing with L beyond 2·U_num; an R-I(c) T1 failure.
- **Inconclusive (class I):** non-asymptotic U_disc classification of ε_h; abs(ε) ≤ 2·U_num; fewer than 3 resolved lengths for an L-scaling statement.

## 6. What existing literature predicts, and what would be new

- **Predicted by prior art:**
  - leading-order GCA with O(ε²) aspect-ratio corrections (Ward 1988 ch. 1, E1) — consistent with α ∝ L⁻²;
  - end/interface-layer effects on the averaged channel conductivity (Ward, Odeh & Cohen 1990, **E2 only**);
  - linear-regime I–V including end effects (Ward 1988 §6.4, **E3 only**).
  - Therefore H-A and H-B may already be derivable from Ward. **Their numerical confirmation alone would be routine.**
- **What would be a meaningful new contribution (only after the Ward full texts are checked):**
  - (i) If H-C holds: a verified, uncertainty-qualified **reservoir-independent** reduction coefficient c_σ(V_G, L), with its bias dependence, compared quantitatively with Ward's asymptotics, including where the asymptotics fail.
  - (ii) If H-C fails (σ numerical or reservoir-dependent): a demonstration that GCA "reduction error" measured against idealized 2-D references is **dominated by the reference's reservoir model** at practically relevant L. Any reduced-model validation against such references is then conditional on an arbitrary boundary choice. This is a negative, methodological result: relevant to the verification of compact and reduced models, and modest in novelty.
- **Not new:** confirming that GCA error decreases with L; exponential short-channel behavior; 1/L access terms (Research_03).

## 7. Decision gate

Run the 48-state stop/go stage only after:
- (a) the independent review of A–D (R05 rule 6);
- (b) R-I(c) is implemented and passes T1;
- (c) the Ward 1990 full text (and ideally thesis §5–6) has been read, so that H-A/H-B can be checked against existing derivations before computing.
