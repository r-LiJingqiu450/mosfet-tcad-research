# Ward 1988 / Ward–Odeh–Cohen 1990 — Full-Text Scientific Audit

**Date:** 2026-10-31. **Branch:** `phase3/C5-literature-2026-10-31`, from C-4 head `d3cc880`.

**Supersedes:** `C4_Ward_FullText_Audit_2026-10-30.md`. That file recorded that access was blocked; it is kept unchanged as a record.

## 0. Sources, evidence level, page map

| Code | Source | Copy audited | Read |
|---|---|---|---|
| **[W88]** | M. J. Ward, *Asymptotic Methods in Semiconductor Device Modeling*, PhD thesis, Caltech, 1988 (DOI 10.7907/08TK-NM84) | User attachment `Ward_jw_1988.pdf`, 156 PDF pages, SHA-256 `a3a00cdd…2dd35`. PDF page ≈ printed page + 7 | Printed pp. 15–16, 83–90, 95, 105–106, 120–123, 131–141, 143–144 as page images (E1). OCR used only to locate passages. **Printed p. 142 is missing from the scan** |
| **[WOC90]** | M. J. Ward, F. M. Odeh, D. S. Cohen, "Asymptotic methods for metal oxide semiconductor field effect transistor modeling", *SIAM J. Appl. Math.* 50(4):1099–1125, 1990 (DOI 10.1137/0150066) | User attachment `WARsiamjam90.pdf`, 27 pp., SHA-256 `6696444b…025d` | pp. 1099–1108 and 1118–1125 as page images, equations checked on images (E1). pp. 1109–1117 (§4 equilibrium asymptotics, §5 start) via OCR only (E1-OCR) |
| [W92] | M. J. Ward, "Singular perturbations and a free boundary problem in the modeling of field-effect transistors", *SIAM J. Appl. Math.* 52(1):112–139, 1992 (DOI 10.1137/0152007) | Crossref abstract only | **E2.** Buried-channel FET; source/drain regions treated with a Howison–King free-boundary solution; "for small drain biases … a current-voltage relation for moderately short devices" |
| [RW91] | L. G. Reyna, M. J. Ward, "The asymptotic solution of some free boundary problems in device modeling and their applications", NASECODE VII, 1991, pp. 232–239 | Ward's publication list | E3 |
| [WRO94] | M. J. Ward, L. G. Reyna, F. M. Odeh, "Some examples of singular perturbation problems in device modeling", IMA Vol. Math. Appl. 59, 1994, pp. 389–406 | Ward's publication list | E3 |
| [T00] | Y. Taur, "MOSFET channel length: extraction and interpretation", *IEEE T-ED* 47(1):160–170, 2000 (DOI 10.1109/16.817582) | IBM Research abstract page | E2. 2-D model I–V used to interpret the extracted channel length as the "injection points" where current spreads from the surface into the bulk source/drain |

**Evidence levels:**
- E1: full text read.
- E2: abstract.
- E3: metadata only.

Every statement below about [W88] or [WOC90] cites a printed page and an equation number. Statements about [W92], [RW91], [WRO94] and [T00] are limited to their stated level.

Ward's scaled variables are used where the audit compares with Ward:
- w = ψ/(V_T lnλ), with ψ measured from the intrinsic level;
- λ = N_A/n_i;
- ε = (L_d/L)(lnλ/λ)^{1/2}, the aspect ratio, with L_d the intrinsic Debye length;
- y ∈ [0, 1] is the channel coordinate; x is depth.

**For the F0 decks:**
- n_i = 10¹⁰ cm⁻³ (`phase3/c1/deck.py`), so λ = 10⁶ and lnλ = 13.82.
- L_d = 40.9 µm, so **ε = 0.152/L[µm]**:
  - L = 1 µm: ε = 0.152;
  - L = 3 µm: ε = 0.051;
  - L = 10 µm: ε = 0.015.
- [W88] p. 16 and [WOC90] p. 1104 adopt ε ≤ 0.01 as the quasi-1-D validity condition. Every C-3/C-4 length is outside it.

## 1. What the two works contain: equation-by-equation record

### 1.1 Formulation ([W88] ch. 1, ch. 5; [WOC90] §2)

| Item | Location | Content | Assumptions | Our F0 decks |
|---|---|---|---|---|
| Scaled DD system | [W88] (1.5a–c), (1.6a,b), p. 15; [WOC90] (2.3a–c) | Poisson with Boltzmann carriers. Electron/hole continuity in quasi-Fermi form. ε appears only as ∂_y → ε∂_y | Boltzmann statistics; SRH term present in (1.5) | Same |
| Steady state without holes or recombination | [W88] (5.1a–c), p. 84 | Hole current and recombination ignored. Hole density in Poisson from φ_p = 0 | "As in all analytical modeling … under moderate source-drain biases" | **Identical to F0** (holes frozen at equilibrium, no G–R, constant μ_n) |
| φ_n boundary conditions | [W88] (5.2a,b), p. 84; [WOC90] (2.4a,b), p. 1104 | φ_n = 0 on the source line BG (y = 0) and v̄_ds/lnλ on the drain line CH (y = 1). ∂φ_n/∂x = 0 at x = 0 and at x = x\*, the scaled n-well depth | n-wells are perfect φ_n reservoirs **at the metallurgical boundary** | Our R-I decks are also φ_n-Dirichlet reservoirs. The location and electrostatics differ (row below) |
| **ψ boundary conditions on BG, CH** | [WOC90] (2.5a,b), p. 1105 | **w = w_bi = ln(N_D/n_i)/lnλ on y = 0, and w = w_bi + v̄_ds/lnλ on y = 1.** The n⁺ well is an equipotential at the built-in level | λ₊/λ ≫ 1 | **Different in kind.** R-I(a): ψ = 1-D slice data. R-I(b): Neumann. R-I(c): gated strips, ψ natural. In all three the **equilibrium is laterally uniform**. In Ward's it is not: w jumps from the channel surface value to w_bi over an O(ε) end layer |
| Oxide | [WOC90] (2.6), p. 1105; [W88] p. 16, printed pp. 105–106 | 1-D oxide drop as a Robin condition, valid **only** for y/ε ≫ 1 and (1−y)/ε ≫ 1 (d = t_ox/L = O(ε)). Second-order oxide expansion in d: valid if ψ_s″d² ≪ 1 ([W88] p. 106) | Mid-channel only | We solve 2-D Laplace in the oxide, which is higher fidelity. Ward's O(d²) oxide terms are part of our "interior" correction (Prop. 3, term A) |
| Current | [W88] (5.3)–(5.4), p. 85; [WOC90] (2.7)–(2.8), p. 1105 | I_ds = channel-averaged tangential current; constant in y by the no-flux conditions | — | Same definition as our eq. (1) (Perturbation Analysis §1) |
| Validity of quasi-1-D analysis | [W88] p. 16; [WOC90] p. 1104 | ε ≲ 0.01, i.e. L ≳ 10 µm for λ = 10⁶. "Although this criterion is rather arbitrary, comparison with numerical simulations shows that it is quite satisfactory in practice" ([W88] p. 16) | No simulation data are shown in the pages read | **Not verified by Ward.** Testable (§3 of the decision report) |

### 1.2 Outer (mid-channel) expansion ([WOC90] §3; [W88] §5 intro)

| Item | Location | Content | Relevance |
|---|---|---|---|
| Regular expansion | [WOC90] p. 1106 and (3.1a,b); [W88] p. 85, (5.5a,b) | w = w⁰ + ε²w¹ + …, φ = φ⁰ + ε²φ¹ + … | **The interior correction is O(ε²) ∝ L⁻²** |
| φ_n⁰ = φ_n⁰(y) | [W88] p. 86, after (5.5b); [WOC90] §3 | The no-flux conditions make the leading quasi-Fermi potential independent of x | The basis of GCA/Pao–Sah |
| Second-order equation | [W88] (5.6), p. 86; [WOC90] (3.2) | ∂²φ¹/∂x² + (…)∂φ¹/∂x = B(x, φ⁰), forced by φ⁰″ and φ⁰′. Interpreted as "one step of Gummel's algorithm" | Set up but **never evaluated for the current**: there is no O(ε²) current correction anywhere in either work |
| Solvability → leading current | [WOC90] (3.3)–(3.6); [W88] (5.9)–(5.13) | I_ds ~ −I_c∫μφ⁰′Q_c lnλ dy as ε → 0 (Pao–Sah form); Q_c defined in (3.7) | **This is our I_GCA.** Our 1COL Gauss–Legendre GCA evaluates (3.6) with the exact 1-D charge instead of the λ → ∞ asymptotic Q_c |
| End conditions of the outer problem | [W88] (5.15), p. 89; [WOC90] p. 1108 | φ_n⁰(0) = 0, φ_n⁰(1) = v̄_ds/lnλ. "In contrast to the potential, boundary layers in φ_n are not anticipated near the source and drain" | Ward's statement is the analytical counterpart of our Prop. 1 premise (u₁ = y/L) |
| **Brews' interpretation** | [W88] p. 89 (citing Brews [5]); [WOC90] p. 1108 | "In conventional GCA modeling … the inner regions near the source and drain are neglected. … v̄_ds is not the bias applied on the drain contact, but rather should be interpreted as the bias a few depletion widths away from the drain." | **Establishes, qualitatively, that GCA-vs-device comparisons depend on what is assumed in the end regions (our K3c)** |
| Open problem stated in 1990 | [WOC90] p. 1108 | "The influence of the potential in the inner regions on the current flow is currently under investigation." | The journal paper does **not** contain the end correction |
| Subthreshold | [W88] (5.16)–(5.18), p. 90 | φ_n⁰ = −ln(1 − (1 − e^{−v̄_ds})y)/lnλ; I_ds ∝ (1 − e^{−v̄_ds}), exact in v̄_ds at leading order in ε | A closed-form V_D dependence, but at leading order in ε only; no end effects |
| Small-θ expansion | [W88] (5.24a,b), p. 95 | Expansion in a **mobility-degradation** parameter θ (c = c₀ + θc₁) | **Not** a small-V_D or small-ε expansion. Not related to Props. 1–2 |
| Large-v_ds breakdown | [W88] p. 120 | The outer expansion fails where εφ_n⁰′ = O(1), i.e. for v̄_ds ≳ lnλ/2 at the drain | The finite-bias limit of the theory. Irrelevant at the C-4 biases (V_D ≤ 50 mV, v̄_ds ≤ 2) |
| Linear/saturation I–V | [WOC90] §6, (6.1)–(6.11), Figs. 10–12 | Leading order only, evaluated by Simpson quadrature | No end correction |

### 1.3 Two-dimensional end regions ([W88] ch. 6; no counterpart in [WOC90])

| Item | Location | Content | Status in Ward |
|---|---|---|---|
| 2-D equilibrium near the drain | [W88] (6.1), p. 121; Fig. 6.1, p. 122 | Poisson in region I (p channel) and region II (n⁺ well, λ₊/λ ≫ 1). Gate edge at the metallurgical junction, no overlap | Formulation |
| **Assumed surface potential** | [W88] p. 121 | "the surface potential … in the vicinity of the drain is **taken to be** w(0,ỹ) = w_s + (w_bi − w_s) exp(−πεỹ/2d)" | **Phenomenological.** The oxide Laplace problem is not solved |
| Semi-infinite well | [W88] p. 121 | n-well depth assumed large enough that the 1-D solution as x → ∞ is the boundary condition. "a rough estimate" deferred to §6.5 | Assumption. Its justification is in §6.5, which is affected by the missing p. 142 (§2) |
| Free-boundary reduction | [W88] pp. 121–131, Fig. 6.2, (6.13)–(6.17) | Leading-order match by C¹ patching: free-boundary problem (u = −1, ∂u/∂n = 0 on the depletion edge), solved by the Polubarinova-Kochina method for 90°, 270° and 360° corners. Plausibility argument via a penalty term, "not proven rigorously" (p. 122) | Leading order in λ only; errors O(1/lnλ) |
| Normal derivatives, 90° corner | [W88] (6.18)–(6.20), pp. 132–133 | Elliptic-integral parametrization (s ∈ (1, ∞)) | Used by §6.4 |
| Inversion layers | [W88] (6.21)–(6.27), pp. 133–136 | d_∞ ~ √2(2 + ln(lnλ)/lnλ)^{1/2}; matching | — |
| **Inner corner** | [W88] (6.28), pp. 136–137 | O(1/lnλ) corner region → Liouville-type problem; "such an explicit solution has not been found" | **Unresolved** |
| Mobile charge near the ends | [W88] (6.29)–(6.32), pp. 137–138 | Composite charge Q_c = Q_cd + Q_cs − Q_m; source/drain symmetric at equilibrium | Leading order; excludes the inner corner |
| Weak inversion | [W88] p. 138 | Not constructible (P–K needs constant data). Numerical solution of the inner regions "is left to future investigators" | **Explicitly open** |
| **§6.4, governing system** | [W88] (6.33), p. 138 | w_xx + ε²w_yy = λ⁻¹(e^{(w−φ_n)lnλ} − e^{−wlnλ}) + 1; **lnλ φ_n⁰′ Q(y) = −I_ds/I_c**; Q(y) = −(lnλ/λ)^{1/2}∫₀^{x\*} e^{(w−φ_n⁰)lnλ} dx | **φ_n = φ_n⁰(y) is kept inside the O(ε) end layers.** That is quasi-1-D current flow with no current crowding or spreading into the well. Called "crucial" on p. 140 |
| §6.4, end profiles | [W88] p. 139 | v_s(ỹ) = w_sl + (w_bi − w_sl)exp(**−πỹ/2d**); v_d analogous with v̄_ds/lnλ added | Phenomenological. **Inconsistent with p. 121, which has the factor ε in the exponent.** As printed on p. 139, the decay length would be 2εt_ox/π, which depends on L and makes the end term O(ε²). The p. 121 form gives an L-independent 2t_ox/π. We take p. 121 as intended; both are evaluated (§4) |
| §6.4, source-end charge | [W88] p. 139 | Q_cs(ỹ) ~ −√2 e^{v_s lnλ/2}(√(1+α_l²) − α_l), with α_l = α₀ₗ(lnλ)^{1/2}λ^{(1−v_s)/2}, α₀ₗ(s) = (√2d_∞ₗ/(π√s))K((s−1)/s), ỹ(s) = (2d_∞ₗ/(π√s))K(1/s). Valid for ỹ lnλ ≫ 1 | The 90° corner is between the interface and the well. α₀ₗ(ỹ) carries the **reduction of depletion charge next to the junction**, i.e. the asymptotic form of charge sharing |
| §6.4, drain-end charge | [W88] p. 140 | Q_cd(ỹ : v̄_ds), with d_∞r = √2(2 + v̄_ds/lnλ + ln(lnλ)/lnλ)^{1/2}. Invalid for large v̄_ds: the depletion region exceeds the well depth, and εφ_n⁰′(1) = O(1). "there appears to be no remedy" | Finite-bias limits stated by Ward |
| **§6.4, result** | [W88] (6.34), p. 141 | Q_u(y) = Q_cs(y/ε) + Q_cd((1−y)/ε : v̄_ds) + Q_m(y) − Q_c(w_sl : 0) − Q_c(w_sr : v̄_ds/lnλ). With (6.33), "the problem of determining the I–V curves in the linear regime has been **reduced to quadrature**" | **Never evaluated, never compared with 2-D numerics, and not in [WOC90].** [W92] (E2) published the analogous construction for **buried-channel** FETs |
| §6.5, 270° corner | [W88] (6.35), p. 141; p. 142 **missing**; (6.36)–(6.37), p. 143; p. 144 | The n⁺–p junction corner away from the interface. Same d_∞. Well-depth estimate on p. 144: depth ≳ 0.5 µm for λ = 10⁶ | Partly unverifiable (§2) |

### 1.4 What neither work contains (checked against all pages read)

1. No expansion in small V_D beyond leading order. No statement that the current is exact through O(V_D) or O(V_D²) for any end condition.
2. No mirror-symmetry or antisymmetry argument about the end layers.
3. No evaluated O(ε²) current correction (interior or end). (5.6)/(3.2) is set up but not used for the current.
4. No analysis of ideal reservoirs placed away from the metallurgical junction, with laterally uniform equilibrium (our R-I(a), (b), (c)).
5. No numerical solution of the 2-D inner regions. No quantitative check of (6.34). No quantitative check of the ε ≤ 0.01 criterion. ([W88] p. 16 asserts agreement with simulations without showing it in the pages read.)

## 2. The missing printed p. 142

| Conclusion | Depends on p. 142? | Cross-check with [WOC90] | Status |
|---|---|---|---|
| §6.3 90° corner, (6.18)–(6.32) | No (pp. 132–138 complete) | [WOC90] has no ch. 6 material | Read in full |
| **§6.4, (6.33)–(6.34)** | **No.** pp. 138–141 are complete and §6.4 ends on p. 141 | [WOC90] p. 1108 confirms the end correction is *not* in the paper | Read in full; **unaffected** |
| §6.5 derivation of the 270° outer normal derivatives leading to (6.36)–(6.37) | **Yes.** The derivation and any definitions between (6.35) on p. 141 and (6.36) on p. 143 are on p. 142 | **Not possible:** [WOC90] does not treat the 270° corner | **Affected.** (6.36)–(6.37) are taken as printed; their derivation is unverified |
| Ward's p. 143 asymptote x(s) ~ −(3d_∞/8)s^{−3/2} as s → ∞ | Possibly (notation fixed on p. 142) | — | **Affected.** With the parameter convention for K, which reproduces Ward's other stated limit ∂u/∂y → −d_∞ as s → 1, (6.36) as printed gives **−(d_∞/8)s^{−3/2}**, a factor 3 smaller. The −1/3 power law ∂u/∂y ∝ x^{−1/3} agrees. Whether (6.36) or the asymptote carries a typo cannot be decided without p. 142. Not reconstructed. Script: `phase3/c5/ward64_quadrature.py`, `check_6_36_large_s` |
| p. 144 well-depth criterion (≳ 0.5 µm at λ = 10⁶) | Yes (built on §6.5) | — | **Affected.** The semi-infinite-well assumption behind §6.3–6.4 (p. 121) is justified only through this estimate. **Consequence for the project:** any 2-D test must verify well-depth insensitivity numerically, not assume it |

No equation from p. 142 is reconstructed anywhere in this audit.

## 3. Comparison with the project's findings (C-3, C-4)

### 3.1 The decisive structural difference

The project's Prop. 1 (T1 identity) assumes ∇·(n₀∇u₁) = 0 with **n₀ = n₀(x) laterally uniform**. Every R-I deck satisfies this by construction:
- (a): Dirichlet slice data;
- (b): Neumann;
- (c): gated strips.

That is why G_2D = G_GCA exactly, and why the C-3/C-4 discrepancies first appear at O(V_D²) (Prop. 2) with L⁻² coefficients (Prop. 3).

In Ward's formulation, the n⁺ wells impose w = w_bi at the ends ([WOC90] (2.5); [W88] (6.1), Fig. 6.1). The equilibrium electron density then has **O(ε)-wide end layers**:
- **surface-potential relaxation** on the oxide fringe scale 2t_ox/π ([W88] p. 121);
- **depletion-charge reduction** next to the junction on the depletion scale d_∞ ([W88] α₀ₗ(ỹ), p. 139).

Through (6.33), these change the resistance ∫dy/|Q_u| at **O(ε) and already at v̄_ds → 0**. So in Ward's physically structured problem:
- Prop. 1 is false;
- the leading GCA discrepancy is a **conductance** correction of relative size O(ε) ∝ 1/L, with one sign (G_2D > G_GCA);
- in engineering language it is an intrinsic, gate-bias-dependent effective-length reduction ΔL(V_G).

### 3.2 Magnitude: a-priori evaluation of Ward's (6.34)

Script and outputs: `phase3/c5/ward64_quadrature.py` → `results/phase3/C5/ward64/ward64_predictions.json`.

This is an evaluation of Ward's published formulas, not a validation. Inputs: F0 parameters, N_D = 10²⁰ abrupt wells, interpretive choices as stated in the script.

| V_G (V) | w_s lnλ | ΔL_W per end (nm), base | Variants: modulus-K / literal p. 139 | Charge-sharing only / fringe only | Predicted 1 − G_GCA/G_2D at L = 1 / 3 / 10 µm |
|---|---|---|---|---|---|
| 1.0 | 13.9 (threshold; §6.4 marginal) | 143 | 157 / 136 | 136 / 35 | 0.29 / 0.095 / 0.029 |
| **1.25** | 17.4 | **98** | 107 / 88 | 87 / 26 | **0.20 / 0.065 / 0.020** |
| 1.5 | 18.6 | 71 | 77 / 60 | 59 / 22 | 0.14 / 0.047 / 0.014 |
| **2.0** | 19.8 | **47** | 50 / 37 | 35 / 16 | **0.093 / 0.031 / 0.009** |

**Comparison with the C-4 numbers** at the same point (L = 1 µm, V_G = 1.25 V):
- R-I(a)/(b) discrepancy at V_D = 50 mV: |ε| ≈ 5 × 10⁻⁴;
- Ward's §6.4 prediction for the n⁺-well end correction: ≈ 2 × 10⁻¹, already at V_D → 0.

**Even allowing Ward's asymptotic error (O(1/lnλ), O(ln lnλ/lnλ): ~±30 %) and the interpretive spread (±20 %), the physically leading end effect is two to three orders of magnitude larger than everything C-3/C-4 measured.** At L = 1 µm, ε = 0.15 and the O(ε) term is ~20 %, so higher-order terms are probably not negligible there. The prediction is cleanest at L ≥ 3 µm.

**Mechanism split** (within Ward's formulas): the charge-sharing part (α₀ₗ, depletion-width scale) dominates the oxide-fringe part (v_s, t_ox scale) by 3–4×. Charge sharing next to the junction is the long-channel limit of the classical short-channel threshold mechanism, so its *existence* is well established. What Ward adds is an asymptotically systematic form of it.

### 3.3 Claim-by-claim

| Project finding | What Ward contains | Literature status | Decision |
|---|---|---|---|
| **Prop. 1 / T1 identity**: G_2D = G_GCA for laterally uniform equilibrium | φ_n⁰ = φ_n⁰(y), no φ_n boundary layers ([W88] p. 86, p. 89; [WOC90] p. 1108). Ward never considers uniform-equilibrium reservoirs; his ends make the identity false at O(ε) (§6.4) | Elementary in that model class | **Not a contribution.** The audit adds: **the identity holds only because R-I removes the physical end layers** |
| **K3a** ε ∝ V_D² (Prop. 2) | Absent (§1.4 items 1–2) | Not in Ward; elementary | **Rejected as novel** (unchanged) |
| **K3b** coefficients ∝ L⁻² (Prop. 3) | ε² ordering of interior corrections ([WOC90] (3.1); [W88] (5.5)–(5.6)); end layers of width O(ε) | **Established in substance**: L⁻² is the ε² ordering of Ward's expansion | **Rejected** |
| **K3c** end-condition dependence of the GCA comparison | Brews' remark ([W88] p. 89; [WOC90] p. 1108); "influence of the inner regions … under investigation" ([WOC90] p. 1108); §6.4 makes the end layers control the O(ε) current | **Established qualitatively.** Our quantitative values concern a model class in which the physically leading O(ε) term is identically zero | **Rejected as a contribution** |
| **K7** mirror-symmetry selection rule (TC-1) | Absent | Not found; elementary; meaningful only inside the ideal-reservoir class | **Unproven, not worth investment.** In physical devices the O(ε) conductance term dominates. Linear-response conductance is source/drain-swap symmetric by reciprocity anyway, so asymmetry enters only at O(V_D) |
| **H-R** reservoir-truncation hypothesis | Absent | — | **Moot for a physical contribution.** It concerns the O(ε²V_D²) residue of an unphysical class |
| **K1/K4** verification methodology | — | Standard components (Roache; Oberkampf & Roy; adjoint DWR) | Methodology only (unchanged) |
| **K5** validity length L\*(η) | Ward asserts ε ≲ 0.01 ([W88] p. 16), "quite satisfactory in practice", without data in the pages read | Open as a *quantified, verified* statement | **Rehabilitated only with an n⁺-well reference.** Ward's §6.4 implies L\* ≈ 2ΔL(V_G)/η, which is strongly gate-bias dependent |
| **D1** (GCA error decomposition against idealized 2-D references) | — | — | **Rejected.** The quantity D1 measures is the O(ε²V_D²) residue of a reference class that deletes the O(ε) end physics Ward identifies as leading. It is reference-dependent at its own order (C-4) and physically subleading by 10²–10³ (§3.2) |

## 4. Assumption audit: Ward vs our decks

| Assumption | Ward | F0 decks | Consequence for a test of Ward |
|---|---|---|---|
| Boltzmann; no holes in current; no G–R; constant μ | Yes ([W88] (5.1)) | Yes | Comparable without change |
| Oxide | Robin, mid-channel only; phenomenological profile at the ends ([W88] p. 121, p. 139) | Full 2-D Laplace | The 2-D reference **tests** Ward's profile. It cannot inherit it |
| S/D | Abrupt n⁺ wells, λ₊/λ ≫ 1, equipotential at w_bi; semi-infinite depth; no gate overlap | Ideal reservoirs, laterally uniform equilibrium | **R-I decks cannot test §6.4.** An n⁺-well deck (R-II) is required |
| φ_n in end layers | φ_n = φ_n⁰(y), quasi-1-D ((6.33); "crucial", p. 140) | Full 2-D continuity | The 2-D reference contains current spreading into the well (the effect [T00] uses to interpret extracted L). Ward's quadrature does not |
| λ → ∞ asymptotics | Leading order in 1/lnλ; inner corner (6.28) unresolved | Exact | Expected agreement only within ~±30 % at lnλ = 13.8 |
| ε → 0 | ε ≲ 0.01 | ε = 0.015–0.15 | The L-dependence must be used to separate O(ε) from O(ε²) |

## 5. Audit conclusions

1. **[WOC90]** contains the formulation, the outer expansion and the leading-order (Pao–Sah) I–V. It **explicitly excludes** the end-region influence on the current (p. 1108).
2. **[W88] §6.4** contains the only end-corrected linear-regime theory in the two works, eq. (6.34) with (6.33). It is:
   - an **O(ε) mobile-charge correction** from n⁺-well end layers;
   - built on a **phenomenological** surface-potential profile;
   - **quasi-1-D** in the current;
   - leading order in 1/lnλ, with the inner corner unresolved;
   - "reduced to quadrature" but **never evaluated or verified**.
   Ward states that numerical treatment of the inner regions is "left to future investigators" (p. 138). [W92] (E2) published the analogous construction for buried-channel FETs.
3. **None of the C-3/C-4 results appears in either work.** Every one of them is either elementary (Props. 1–2, K7) or implied by Ward's ordering (K3b, K3c). More importantly, **they live in a reference class that sets Ward's leading O(ε) end correction to zero by construction.** This is the literature-backed reason to reject D1.
4. The missing p. 142 affects only §6.5: the 270° corner derivation, a factor-3 inconsistency in a quoted asymptote, and the well-depth estimate. It does **not** affect §6.4 or any conclusion above, except that the semi-infinite-well premise of §6.4 must be checked numerically.
