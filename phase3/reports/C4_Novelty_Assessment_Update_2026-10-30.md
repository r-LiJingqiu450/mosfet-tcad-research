# C-4 Claim-versus-Prior-Art Assessment (update of the C-3 matrix)

**Date:** 2026-10-30/31. **Overall novelty status: `NOVELTY_UNPROVEN`.**

**Several individual claims are downgraded or withdrawn below.** The D1 framing ("the GCA reduction error measured against an idealized 2-D reference") is **withdrawn as posed**. The reason is not that prior art establishes it, which is unknown because the Ward full texts were inaccessible. The reason is that C-4 shows the quantity is **not reference-independent** (Proposition 3; ℓ_res diagnostic; TC-1).

**Evidence levels:**
- E1: full text read.
- E2: abstract.
- E3: metadata.
- P: the project's own derivation, which is not literature evidence.

## 1. Inputs since C-3

| Input | What it shows | Effect on claims |
|---|---|---|
| Ward 1990 / thesis ch. 5–6 | **Not read** (C4_Ward_FullText_Audit: access blocked) | Nothing in their territory can move to POTENTIAL |
| Propositions 1–2 (P) | For mirror-symmetric ideal reservoirs the GCA current is exact through O(V_D²). ε ∝ V_D² follows from an elementary symmetry argument | **The δ² law is not a contribution**, whether or not Ward states it |
| Proposition 3 (P, scaling argument) | The leading coefficient has reservoir-dependent and interior parts of the same order, L⁻² | No reservoir-independent reduction error at leading order |
| TC-1 (registered verdict INCONCLUSIVE) | Mirror-asymmetric decks have a linear term c₁ = ±0.0336 V⁻¹ (L = 1 µm) and ±0.00374 V⁻¹ (L = 3 µm): ∝ L^−2.00, c₁(ab) = −c₁(ba) to 1e-6. Symmetric decks have \|c₁\| at 5–7×10⁻⁶ of that, but at 3.0–4.0σ of a cubic fit with large χ². The registered "≤ 3σ" criterion therefore fails | The mechanism of Prop. 2 is strongly indicated post hoc but **not** confirmed by the registered rule |
| ℓ_res diagnostic (L = 1 µm, 1.25 V, 50 mV) | ε moves monotonically from R-I(b) (+5.70e-4) through R-I(c) f = 0.1 (+2.42e-4) and f = 0.3 (+6.8e-5) to an unresolved value at f = 1 (\|ε\| < 7e-6), while R-I(a) is −5.13e-4 (C4_lres_and_TC1_Results) | The discrepancy is a property of the reservoir specification, within one family of admissible ideal reservoirs |

## 2. Claim matrix (C-4 status)

| Claim | Prior art (level) | Project evidence | C-4 status | Reason |
|---|---|---|---|---|
| **K1** Discrete zero-bias identity as a TCAD code-verification test, with an equilibrated adjoint-based (DWR) error estimate validated against two backward-stable solvers | Exact-solution code verification (Roache; Oberkampf & Roy, E2). Adjoint/dual-weighted-residual output error estimation is standard numerical analysis (Becker & Rannacher; Giles & Süli; E3 from memory of the field, not re-checked here) | 192 points, 6 decks; estimate covers the cross-solver difference at 192/192 | **Downgraded: methodology, not a T-ED contribution** | Assembly of standard techniques; publishable at most as a verification appendix or methods note |
| **K2** Triple-point and bottom-corner contact treatment in DEVSIM | DEVSIM documentation and issues not searched | Corrected decks verified | **Rejected as a scientific contribution** | Implementation detail |
| **K3a** ε ∝ V_D² with no linear term at small V_D | Unknown (Ward, not read) | Prop. 2 (P); C-3 L3 ladders | **Rejected as novel** | Elementary consequence of mirror symmetry |
| **K3b** Reservoir-antisymmetric part ∝ L⁻² | Unknown | C-3 α ∝ L^−1.95 | **Rejected as novel by itself** | Expected scaling of an end layer (Prop. 3) |
| **K3c** Dependence of the small-V_D coefficient on the end/reservoir condition, quantified with verified uncertainty | Ward §6.4 "I–V in the linear regime including end effects" (E3) is the direct competitor | ℓ_res diagnostic; TC-1 | **Unproven** | Cannot be judged without Ward §6.4. Even if new, it is a statement about idealized references |
| **K4** Same-mesh discrepancy metric | Standard practice (E2) | Same-mesh and continuum metrics converge toward each other with refinement | **Downgraded: methodological device** | — |
| **K5** Validity boundary L\*(η) | Brews et al. 1980 (E2/E3); Kim & Yang 1970 (E3) | None | **Unproven, and weakened** | At small V_D the "GCA error" depends on the reference reservoir, so L\* would too unless a physical reference (R-II) is used |
| **K6** Reduction vs structural error ladder | Research_03 | None | **Unproven; premise weakened** | ε_red is reference-dependent, so the ladder needs a physically structured reference |
| **K7 (new)** Symmetry selection rule: for ideal φ_n-Dirichlet reservoirs, the O(V_D²) current discrepancy vanishes iff the reservoir pair is mirror-symmetric; asymmetric pairs give c₁ ∝ L⁻² with c₁(ab) = −c₁(ba) | Unknown (not searched beyond Ward titles) | Prop. 2 (P); TC-1 (post hoc strongly consistent; registered verdict INCONCLUSIVE) | **Unproven** | Likely elementary or known in asymptotic device theory. Its relevance to real devices (asymmetric S/D, LDD) is untested |
| **D1 as originally framed** — "quantitative decomposition of the GCA error against higher-fidelity 2-D references" | — | Prop. 3; ℓ_res; TC-1 | **Withdrawn as posed** | Its central quantity is not defined independently of an arbitrary reservoir specification at the leading small-bias order |

## 3. What would count as a genuine contribution

A contribution must satisfy all of the following:
1. **Literature:** the full texts of Ward 1988 ch. 5–6 and Ward–Odeh–Cohen 1990 are read, together with the classical GCA-correction literature: Pao–Sah, charge-sheet models, Brews et al., Kim & Yang. The claim must not appear there in substance. Notation or implementation differences do not count.
2. **Physical relevance:** the quantity is defined against a reference whose end conditions are physically meaningful (R-II, n⁺ source/drain), or is shown to be insensitive to the reservoir specification over an admissible family. C-4 shows the second does **not** hold for the idealized family at small V_D.
3. **Independent evidence:**
   - verified numerics (T1 v2-grade references, ≥ 4 refinement levels, asymptotic classification);
   - a theory comparison (asymptotic or perturbative) that predicts the numbers rather than fits them;
   - at least one independent implementation or published data set.
4. **Significance for T-ED readers:** the result must change how compact/reduced models are built or validated. One candidate, not yet established: "reduced-model validation against idealized 2-D references can be dominated by the reference's reservoir model; here is a symmetry rule and a quantitative criterion". Whether that is new and useful depends on (1).

**Current position:** none of K1–K7 meets (1)–(4). **`NOVELTY_UNPROVEN`.**
- Not `NOVELTY_REJECTED`: no publication was shown to establish the remaining candidates (K3c, K7), because the decisive texts could not be read.
- Not `NOVELTY_POTENTIAL`: there is no literature comparison.
