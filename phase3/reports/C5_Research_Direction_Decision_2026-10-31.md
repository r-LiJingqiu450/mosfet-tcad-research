# Research-Direction Decision after the Ward Full-Text Audit

**Date:** 2026-10-31. **Branch:** `phase3/C5-literature-2026-10-31`.

**Basis:**
- `C5_Ward_FullText_Audit_2026-10-31.md` (the audit; equation and page references there);
- the C-3/C-4 results;
- `results/phase3/C5/ward64/ward64_predictions.json` (a-priori evaluation of Ward's published §6.4 formulas).

**No simulation was run for this decision.**

## Summary of decisions

| Item | Decision |
|---|---|
| D1 (GCA error decomposition against idealized 2-D references) | **`NOVELTY_REJECTED`, and rejected as a research direction** |
| K3a, K3b, K3c, Prop. 1–3 as contributions | **Rejected.** Elementary, or implied by Ward's ε-ordering and Brews' remark |
| K7 (symmetry rule), H-R, H-S, FB-1 | **Not to be pursued.** Physically subleading model class |
| K1/K4 (verification methodology) | Methodology only. Fallback for a short computational note |
| **New direction C-W:** verification and decomposition of the O(ε) source/drain end correction to the linear-regime conductance (Ward 1988 §6.4) against verified 2-D drift-diffusion with resolved n⁺ wells (R-II) | **Recommended for a small, gated investment. `NOVELTY_UNPROVEN`**, pending literature gate LG-1 (§2) |
| Smallest decisive experiment | **E-W1** (§4): 48 zero-bias linear-response states. Needs separate authorization; not run |

---

## 1. Which of our proposed findings are already established?

| Finding | Where it stands | Decision |
|---|---|---|
| ε ∝ L⁻² scaling of the discrepancy coefficients (K3b, Prop. 3) | Implied by the ε² ordering of Ward's outer expansion: [WOC90] (3.1a,b); [W88] (5.5)–(5.6). End layers are O(ε) wide | **Established in substance** |
| GCA comparisons depend on the end-region treatment (K3c, qualitative) | Brews' interpretation: v̄_ds is the bias "a few depletion widths away from the drain" ([W88] p. 89; [WOC90] p. 1108). Ward 1990 explicitly leaves the inner-region influence "under investigation" (p. 1108). [W88] §6.4 makes the end layers control the current at O(ε) | **Established qualitatively** |
| No φ_n boundary layers; leading φ_n depends only on y (premise of Prop. 1) | [W88] p. 86, p. 89; [WOC90] p. 1108 | **Established** |
| GCA = leading outer term (Pao–Sah) | [WOC90] (3.6); [W88] (5.9)–(5.13) | **Established**; our I_GCA is (3.6) with exact 1-D charge |
| End effects as an intrinsic, gate-bias-dependent effective-length reduction | Charge sharing near junctions is the classical short-channel threshold mechanism. Gate-bias-dependent extracted channel length and its interpretation by 2-D-model I–V are in [T00] (E2) | **Established phenomenon.** The asymptotic form ([W88] (6.34)) is unpublished in journal form for the MOSFET (E1 for [WOC90] p. 1108; E3 for [RW91], [WRO94]) |
| Exact first-order (T1) identity for ideal reservoirs (Prop. 1) | Not stated by Ward, who never uses such reservoirs; elementary | Not established but **not a contribution** |
| ε ∝ V_D², no linear term (K3a, Prop. 2) | Absent from both works; no small-V_D expansion exists there | Not established; **elementary**; rejected |
| Mirror-symmetry selection rule (K7, TC-1) | Absent | Not established; elementary; physically marginal (below) |

**The decisive finding of the audit** is not that Ward published our results. He did not. It is that Ward's formulation exposes **why our results are physically subleading**:
- Every R-I deck has a **laterally uniform equilibrium**. That is exactly why Prop. 1 holds and the C-3/C-4 discrepancies start at O(ε²V_D²).
- With physical n⁺ source/drain wells ([WOC90] (2.5); [W88] (6.1), Fig. 6.1), the equilibrium has O(ε) end layers. The GCA error is then an **O(ε) conductance correction present at V_D → 0** ([W88] (6.33)–(6.34)).
- Ward's own formulas, evaluated for our F0 parameters (audit §3.2), give 1 − G_GCA/G_2D ≈ **0.20 at L = 1 µm and V_G = 1.25 V**. C-4 measured |ε| ≈ **5 × 10⁻⁴** at V_D = 50 mV in the idealized class.
- That gap of 10²–10³ survives every interpretive and asymptotic uncertainty stated in the audit.
- The idealized-reservoir discrepancy is therefore a reference artefact of an order that is invisible in real devices. **D1 cannot be repaired by better numerics; its object is the wrong quantity.**

## 2. Does a scientifically meaningful, potentially original contribution remain?

**Yes, one candidate.** It is narrower than the original ambition and must pass a literature gate.

**C-W.** For a long-channel enhancement MOSFET with abrupt n⁺ source/drain:
1. quantify, with verified 2-D drift-diffusion (R-II), the O(ε) source/drain end correction to the linear-regime conductance, expressed as an L-independent ΔL_end(V_G);
2. test Ward's unpublished asymptotic theory ([W88] (6.33)–(6.34)) against it;
3. decompose the end correction into three parts:
   - (i) depletion-charge reduction near the junction (charge sharing; Ward's α₀ₗ(ỹ));
   - (ii) oxide-fringe surface-potential relaxation (Ward's v_s(ỹ), phenomenological);
   - (iii) 2-D current spreading into the well, which Ward's quasi-1-D (6.33) omits by assumption;
4. state the resulting validity length L\*(η, V_G) that replaces the unverified ε ≤ 0.01 rule ([W88] p. 16; [WOC90] p. 1104).

**Why it is scientifically meaningful:**
- It tests a specific, falsifiable published theory that its author explicitly left unevaluated:
  - "reduced to quadrature" ([W88] p. 141);
  - inner regions "left to future investigators" (p. 138);
  - "under investigation" ([WOC90] p. 1108).
- The pre-registered predictions already exist: audit §3.2 and the predictions JSON.
- The decomposition (i)–(iii) separates a charge effect from a transport effect. That is the reduced-model-vs-higher-fidelity error decomposition the project set out to do, now on a physically defined quantity.

**What could make it non-original (literature gate LG-1, before any simulation):**
- [W92] *SIAM J. Appl. Math.* 52:112–139 (E2) does the analogous small-drain-bias construction for **buried-channel** FETs. Must be read in full to see whether it treats the MOSFET case or compares with 2-D numerics.
- [RW91] NASECODE VII (1991) and [WRO94] IMA Vol. 59 (1994): E3. Must be checked for an evaluated MOSFET end correction.
- [T00] (E2) and the classical channel-length-extraction and charge-sharing literature (to be identified and read; not verified in this audit) must be checked for a **verified** closed-form or asymptotic ΔL(V_G) for abrupt junctions.

**Significance, stated candidly:**
- Even if C-W passes LG-1 and E-W1, the contribution is a verified asymptotic theory of a known phenomenon for a classical (long-channel, abrupt-junction) device. Its main users are reduced-model and compact-model developers and the asymptotic-device-theory community.
- **IEEE T-ED acceptance is unlikely** unless the result changes how channel length or mobility is extracted or validated. That is not yet demonstrated.
- Realistic venues: *Solid-State Electronics*, *IEEE J-EDS*, *J. Computational Electronics*, or, for the asymptotic-verification angle, *SIAM J. Appl. Math.* / *Eur. J. Appl. Math.*
- The target-venue decision should be taken after E-W1, not before.

## 3. Which research direction deserves further investment?

| Option | Verdict | Reason |
|---|---|---|
| D1 and its follow-ups (FB-1, H-R, H-S, K7 tests) | **Stop** | Physically subleading quantity, reference-dependent at its own order, and implied by Ward's ordering where not elementary |
| **R-II via C-W** | **Invest, in two gated steps (LG-1, then E-W1)** | The only candidate with a physically defined quantity, an existing falsifiable theory, and pre-computed predictions. It reuses the project's strongest asset |
| Methods note on K1 (code verification by the zero-bias identity with validated adjoint error estimates) | Fallback only | Solid but assembled from standard components; low significance |
| Broad R-II parameter sweeps or modern short-channel devices | Not now | ε ~ 1 there; no asymptotic theory to test; becomes routine TCAD |

**The project's strongest asset is the T1 v2 machinery**: a verified, adjoint-error-estimated zero-bias linear-response computation of G. It is exactly the tool C-W needs, because **the decisive quantity is a V_D → 0 conductance**. No finite-bias run is required.

## 4. The smallest decisive experiment: E-W1

**Status:** design only. Needs separate authorization, because it introduces the R-II deck class, previously flagged as needing Package E authorization. A formal pre-registration document is to be frozen before any run; the predictions below are already fixed in `ward64_predictions.json`.

**Decks** (F0 physics unchanged):

| Deck | Definition | Role |
|---|---|---|
| **R-II(h)** | Abrupt box n⁺ wells, N_D = 10²⁰ cm⁻³, depth h. Gate edge at the metallurgical junction (no overlap: Ward's geometry). Wells extend 0.5 µm beyond the gate edge, with top ohmic contacts carrying φ_n = 0 / V_D. 2-D oxide; no gate over the wells. Body contact at x_max = 1.2 µm | Physical reference |
| **R-I(w)** | Channel only. End lines carry ψ = ψ_bi(n⁺) and φ_n Dirichlet: Ward's domain truncation, [WOC90] (2.4a), (2.5) | Isolates the equipotential-well assumption |

**Quantities:**
- G_2D: zero-bias linear-response conductance (T1 v2 x_AC path; E_dwr and E_round reported).
- G_GCA = μN_s/L from the 1COL slice at V_D = 0.
- **ΔL_2D(L) = (L/2)(1 − G_GCA/G_2D)**, with mesh extrapolation and U_disc (R04 §7.4).
- G_M2 = μ/∫₀^L dy/N_2D(y), using the 2-D equilibrium sheet density N_2D(y). This is Ward's (6.33) with exact 2-D charge.
- Ward's M1 = ΔL_W (pre-computed).

**Matrix (48 zero-bias states):**
- Core: R-II(h = 0.5 µm) × L ∈ {1, 2, 4, 8} µm × V_G ∈ {1.25, 2.0} V × levels 0–3 = **32**.
- Sensitivity at L = 2 µm: R-II(h = 1.0 µm) and R-I(w) × 2 V_G × levels 0–3 = **16**.
  - Well-depth insensitivity must be tested because Ward's justification (p. 144) depends on the missing p. 142.

**Cost:** one equilibrium Newton solve plus one linear solve per state. The estimate is **4–8 CPU-hours**, to be confirmed by a single level-0 smoke state before the matrix. Run through `phase3/c4/jobrunner.py` (checkpointed, bounded).

**Main numerical risk:** the junction–interface corner. Ward finds a singular field there (x^{−1/3}, (6.36)–(6.37); the 90° analogue in §6.3). The n⁺ Debye length is ~0.4 nm. Reduced observed order is possible; the R04 §7.4 classification must be passed, not assumed.

**Pre-registered predictions** (Ward [W88] (6.34) at v̄_ds → 0; base convention; ~±30 % asymptotic uncertainty):

| V_G | ΔL_W per end | 1 − G_GCA/G_2D at L = 1 / 2 / 4 / 8 µm |
|---|---|---|
| 1.25 V | 98 nm (variants 88–107) | 0.196 / 0.098 / 0.049 / 0.025 |
| 2.0 V | 47 nm (variants 37–50) | 0.093 / 0.047 / 0.023 / 0.012 |

Also predicted:
- the sign: G_2D > G_GCA;
- the trend: ΔL(1.25 V) > ΔL(2.0 V);
- the structure: ΔL_2D(L) → constant, with an O(1/L) remainder.

**Proposed decision rules (to be frozen at registration):**

**D-a, O(ε) structure.** Fit ΔL_2D(L) = ΔL_∞ + b/L over L ∈ {2, 4, 8} µm (L = 1 µm as a check). Required:
- ΔL_∞ > 5U;
- |ΔL_2D(8) − ΔL_2D(4)| ≤ 0.15·ΔL_∞ + 2U.

**D-b, intrinsic.** |ΔL_2D(h = 1.0) − ΔL_2D(h = 0.5)| ≤ 0.10·ΔL_2D at L = 2 µm, both V_G.

**D-c, Ward magnitude.** ρ = ΔL_∞/ΔL_W:
- **CONFIRMED** if ρ ∈ [0.7, 1.3] at both V_G and the trend matches;
- **DEFICIENT** if ρ ∈ [0.3, 0.7) ∪ (1.3, 3];
- **REFUTED** otherwise, or if the sign or trend is wrong.

**D-d, mechanism.**
- Δ_spread = ΔL from (G_2D vs G_M2) on R-I(w) is resolved if > 3U.
- Δ_well = ΔL_2D(R-II) − ΔL_2D(R-I(w)) is resolved if > 3U.

**Outcome → action:**

| Outcome | Interpretation | Action |
|---|---|---|
| D-a or D-b fails | No intrinsic O(ε) end correction at accessible L | **Abandon C-W** (criteria A1/A2) |
| D-c CONFIRMED, Δ_spread and Δ_well unresolved | Ward's quadrature is right and complete. The only result is "an unpublished 1988 formula works" | Short note at most; **stop T-ED investment** (A3) |
| D-c DEFICIENT or REFUTED with Δ_spread, Δ_well or the profile ablation resolved and attributable | Identified failure of a named assumption: (6.33) quasi-1-D current, the p. 121/139 profile, the equipotential well, or leading-order charge sharing | **Continue:** this is the publishable core. Next, a separately authorized extension in N_A, t_ox and λ₊ to establish the decomposition and L\*(η, V_G) |
| D-c CONFIRMED with resolved Δ_spread/Δ_well | Ward's charge model is right; the transport and well corrections are separately quantified | Continue with the decomposition paper, as in the previous row |

## 5. Evidence that would justify abandoning C-W

| # | Evidence | When it can arise |
|---|---|---|
| **A0** | [W92], [RW91], [WRO94], [T00] or the charge-sharing/ΔL literature already contains an asymptotic or closed-form ΔL_end(V_G) for an enhancement MOSFET with abrupt n⁺ junctions, **compared against 2-D drift-diffusion** | LG-1, before any simulation |
| **A1** | No L-independent ΔL_2D up to L = 8 µm (D-a fails): the O(ε) structure is not observable at accessible lengths | E-W1 |
| **A2** | ΔL_2D changes by > 10 % with well depth (D-b), or proves sensitive to contact placement: the quantity is not intrinsic | E-W1 |
| **A3** | Ward CONFIRMED with no resolved mechanism split: nothing beyond confirmation | E-W1 |
| **A4** | The R-II reference cannot be verified within ~2× the cost estimate: | E-W1 smoke test or levels 0–3 |
| | - non-asymptotic mesh behaviour at level 3; | |
| | - conservation worse than 1e-10; | |
| | - linear-response error estimate not ≪ ΔL signal | |
| **A5** | After E-W1, a scaling check of ΔL_end to thin-oxide, high-doping parameters shows the effect is < 1 % for every test-structure length used in practice. Then it has no practical consequence for T-ED readers. This is not abandonment of the science, but **retargeting the venue** | After E-W1 |

## 6. What is not done and remains forbidden

- No simulation was launched for this decision.
- No golden baseline, production S1/S2/S3 physics or main branch was touched.
- FB-1, H-R/H-S tests, Package E sweeps and the manuscript are **not** started.
- E-W1 requires:
  - authorization;
  - a frozen pre-registration;
  - independent review of this audit (R05 rule 6), especially the §6.4 transcription and the two interpretive choices in `phase3/c5/ward64_quadrature.py`.
