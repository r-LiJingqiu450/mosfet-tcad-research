# Novelty Challenge — Claim-versus-Prior-Art Matrix (D1)

**Date:** 2026-10-30. **Status: `NOVELTY_UNPROVEN`.**
**Rule applied:** no novelty is inferred from the absence of search hits. Evidence levels are kept separate:
- **E1** — full text read (automated extraction, not checked equation by equation);
- **E2** — abstract or repository summary;
- **E3** — metadata, title or table of contents only.

## 1. Prior-art inventory (this session's additions marked ★)

| # | Item | Level | What it establishes (at that level) | Threat |
|---|---|---|---|---|
| P1 | Ward, Odeh & Cohen, "Asymptotic methods for MOSFET modeling," *SIAM J. Appl. Math.* 50(4):1099–1125, 1990, DOI 10.1137/0150066 | **E2**★ (repository abstract and keywords) | Formal matched asymptotics for small-aspect-ratio, heavily doped MOSFETs. Resolves interface and internal layers and **examines how those layers affect the averaged channel conductivity**. I–V relations uniform across bias regimes. Keywords: "channel conductivity", "interface layer". Open-access PDF exists at CaltechAUTHORS (record 0j6kp-4jb87) but could not be retrieved through this environment's fetch proxy | **Critical.** Overlaps any claim about end/interface effects on channel conductance |
| P2 | Ward, PhD thesis, Caltech 1988 | **E1 ch. 1–3** (earlier session); **E3 ch. 4–6**★ (table of contents) | Ch. 5.2 "Linear, saturation regimes and pinchoff"; **§6.4 "The I–V curves in the linear regime including end effects"**; §6.3 equilibrium potential near the drain in strong inversion; §6.5 a 270° corner. Ch. 1: GCA is the leading-order outer solution; ε ≤ 0.01 used as a 1-D criterion | **Critical.** §6.4 is the same physical question as the finite-bias design below |
| P3 | Ward, "Singular perturbations and a free boundary problem in the modeling of field effect transistors," *SIAM J. Appl. Math.* 52(1):112–139, 1992, DOI 10.1137/0152007 | E3★ (title, citation) | Free-boundary asymptotics for FETs; content unknown (likely depletion/pinch-off) | Medium |
| P4 | Ward, Reyna & Odeh, "Some examples of singular perturbation problems in device modeling," IMA Vol. Math. Appl. 59, 1994, 389–406 | E3★ | Review-type | Medium |
| P5 | C.K. Kim & E.S. Yang, "On the validity of the gradual-channel approximation for field-effect transistors," *Proc. IEEE* 58(5):841–842, 1970 (KAIST record; DOI not verified) | E3★ | Title only. A published validity criterion for the GCA exists | **High** for any "validity of GCA" framing |
| P6 | Brews, Fichtner, Nicollian & Sze, *IEEE EDL* 1(1):2–4, 1980, DOI 10.1109/EDL.1980.25205 | E2/E3 (unchanged) | Empirical long-channel (subthreshold) boundary from 2-D simulation | High for L\* |
| P7 | Z. Ren (adv. Y. Taur), *Non-GCA modeling of near-threshold I–V characteristics of DG MOSFETs*, PhD thesis, UC San Diego, 2020 (eScholarship 4j36q4hs) | E2★ | Non-GCA corrections (a source/drain-induced charge term ΔQ) for **short-channel** DG MOSFETs, validated against TCAD. Near threshold is the hard regime. The Taur group also lists a non-GCA saturation model for long-channel DG (E3) | Medium. Different geometry and regime, but "quantify non-GCA corrections vs TCAD" is an active line |
| P8 | P. Habaš, "A physics based analytical MOSFET model with accurate field dependent mobility," *Solid-State Electron.* 33, 923–933, 1990 | E1★ (automated) | Assumes the GCA; error analysis only against its own quasi-1-D numerics; **no** quantification of GCA validity vs 2-D | Low |
| P9 | Vandorpe et al., *Solid-State Electron.* 1972 (2-D MOS numerical analysis) | E3 | 2-D vs 1-D comparisons are old | Low–medium |
| P10 | Lu (adv. Taur), PhD thesis, UCSD 2006 | E1 partial (earlier) | Single-L GCA vs 2-D comparisons are routine | Medium |
| P11 | V&V literature: Roache 1994/1998; Oberkampf & Roy 2010; Eça & Hoekstra 2014; Roy 2005 | E2 | Exact/manufactured solutions for code verification; GCI/least-squares uncertainty; model-form vs numerical error separation | High for any "method" novelty |

**Literature checks still open:**
- **L-1:** full text of P1 and P2 ch. 5–6. **Now concretely locatable:** both are open-access Caltech PDFs. Needed: a manual download and attachment.
- **L-2:** Brews et al. 1980 full text.
- **New L-9:** P5 full text.
- **New L-10:** the Taur-group non-GCA papers (*Solid-State Electron.*, 2019–2023) beyond abstracts.

## 2. Claim-versus-prior-art matrix

| Candidate claim | What prior art already establishes | What remains insufficiently quantified (to our knowledge, search-limited) | Proposed new observation or method | Supporting experiment | Refuting observation | Why it would matter | Evidence still missing |
|---|---|---|---|---|---|---|---|
| **K1** The exact **discrete** zero-V_D identity G_2D,h = (qμ/L)·Σn₀w as a code-verification test of 2-D TCAD reference decks | *Continuum:* linear-region conductance = qμN_s/L for a long channel is textbook (P2 §5.2, E3). *Method:* exact solutions for code verification are standard (P11) | A deck-level identity that exposes boundary-condition implementation errors at the 1e-6 level in a production TCAD code | The test plus its resolving power, and the triple-point defect it exposed (−94 % → round-off) | C-3/C-4 T1 (done: 36 + 36 configurations) | The identity already used as a TCAD verification test in the literature, or the defect class already documented for DEVSIM | Reference-quality 2-D decks are a prerequisite for any GCA-error study | A literature check on TCAD code-verification practice (L-4/L-7). **Novelty is likely modest: a methods note or appendix, not a T-ED contribution** |
| **K2** A correct discrete treatment of Si/SiO₂/contact triple points when the TCAD assembly drops interface pairs at contact nodes | FV interface/contact handling is standard. The DEVSIM behavior is in its source | Whether this pitfall is documented for DEVSIM users | Diagnosis plus a continuity tie that adds no artificial constraint; the open R-I(b) circuit-contact shared-node item | Corner study (C-2), T1 continuity and Gauss checks | DEVSIM docs or issues already describe it | Practical, for TCAD verification | DEVSIM issue tracker search. **Not a scientific contribution by itself** |
| **K3** Small-V_D structure of the reduction discrepancy in the idealized references: r(δ) = c·δ² with no linear term; a reservoir-antisymmetric part ≈ L⁻² that dominates at L = 1 µm; peak near threshold | **P1 (E2):** interface-layer effects on the *averaged channel conductivity* are analyzed. **P2 §6.4 (E3):** linear-regime I–V including end effects. Aspect-ratio corrections at O(ε²) are predicted by the ε-expansion (P2 ch. 1, E1) | Whether P1/P2 give the coefficient, its dependence on the reservoir (end) conditions, and its bias dependence; whether these were ever checked numerically with uncertainty | Verification-aware numerical quantification of c(V_G, L) for two specified reservoir models, separating reservoir-antisymmetric and symmetric parts, with mesh-extrapolated uncertainty | The finite-bias design (Finite_Bias_Discrepancy_Experiment_Design, §4) | P1/P2 already state the δ² law with coefficients and end-condition dependence. **Or** the symmetric part is shown to be purely numerical (lateral discretization) and the antisymmetric part fully determined by the reservoir choice. Then K3 reduces to an artifact of the reference specification | Would tell compact-model developers how much of the "GCA error" seen against 2-D references is a property of the *reference's* boundary model | **Full text of P1 and P2 §5–6 (decisive).** Lateral-refinement data (being run) |
| **K4** A **same-mesh** discrepancy metric (discrete GCA on the identical depth lines vs 2-D) cancels the common discretization error and resolves reduction discrepancies far below the continuum-reference uncertainty | Cancellation of common errors in differences is standard numerical practice. Model-form vs numerical separation is standard V&V (P11) | Its use for reduced-model vs 2-D TCAD comparisons, with an explicit account of the residual (lateral) discretization error | ε_red,h with level-wise extrapolation | fb_preview (L = 1 µm, level 0): 6e-10 vs 1.2e-4 at V_G = 0.8 V; lateral-refinement runs pending | The residual lateral discretization error of the same-mesh metric is as large as the continuum-reference error, i.e. no gain | Makes reduction errors measurable where continuum comparisons are inconclusive | Literature check on "same-discretization model comparison" in TCAD and CFD (L-7). **Likely methodological, not standalone** |
| **K5** An uncertainty-aware validity boundary L\*(η) per bias regime | P6 (empirical subthreshold boundary); P5 (a validity criterion, content unknown); P2 ch. 1 (ε ≤ 0.01) | Above-threshold and bias-resolved boundaries with uncertainty intervals | OC2/OC4 intervals | Not yet authorized (needs T3-A data) | P5/P6 already give bias-resolved boundaries with uncertainty | Practical design guidance | Full texts of P5 and P6; the pilot |
| **K6** Separation of reduction vs structural error (R-I → R-II ladder) | Known components: TLM 1/L access, DIBL/λ scaling (Research_03). P7 models S/D-induced charge for short channels | The size of the reduction part within the total GCA error of a structured device, with uncertainty | Attribution algebra, R04 §5.4 | Package E (not authorized) | ε_red is negligible next to ε_str for all relevant L, which makes the ladder uninformative; **or** S_res(a,b) is comparable to ε_red, which makes ε_red ill-defined (the K3 risk) | Would locate where reduced-model error originates | Everything beyond R-I |

## 3. Assessment

1. **Nothing is novel yet.** K1 and K2 are verification infrastructure. K4 is a methodological device. K5 and K6 have no data.
2. **The most interesting new observation (K3) sits squarely in P1/P2's territory.**
   - Ward's 1990 abstract says the paper analyzes how interface/boundary layers affect the averaged channel conductivity.
   - The thesis has a section on linear-regime I–V **including end effects**.
   - Until those full texts are read, K3 cannot be claimed. The prior-session conclusion that D1 can at most claim a *verification-aware numerical quantification* stands, and is now sharper: the quantification must be shown to go beyond what P1/P2 derive (coefficients, end-condition dependence, uncertainty).
3. **A serious framing risk for D1 itself.** At L = 1 µm the finite-bias discrepancy between M0 and R-I(a) and between M0 and R-I(b) has **opposite sign and nearly equal magnitude**. "The reduction error of the GCA" is therefore not well defined independently of the reference's reservoir model there. If this persists at the lengths of interest, D1's headline quantity ε_red is a property of the chosen reference, not of the GCA. The design below makes this falsifiable.
4. **Status: `NOVELTY_UNPROVEN`.** Moving to `NOVELTY_POTENTIAL` requires:
   - (i) the full texts of P1 and P2 ch. 5–6 (and ideally P5) read, and K3/K4 shown to go beyond them;
   - (ii) the lateral-refinement test showing that the reservoir-insensitive part of the discrepancy is physical, not numerical.

## Sources

- [Ward, Odeh & Cohen 1990 — CaltechAUTHORS record](https://authors.library.caltech.edu/31400/) · [SIAM DOI 10.1137/0150066](https://epubs.siam.org/doi/10.1137/0150066)
- [Ward 1988 PhD thesis (Caltech)](https://thesis.library.caltech.edu/440/) · [PDF](https://thesis.library.caltech.edu/440/3/Ward_jw_1988.pdf)
- [Ward 1992, SIAM DOI 10.1137/0152007](https://epubs.siam.org/doi/10.1137/0152007)
- [M.J. Ward publication list](https://personal.math.ubc.ca/~ward/newpub.pdf)
- [Kim & Yang 1970 — KAIST record](https://dspace.kaist.ac.kr/handle/10203/56781?mode=full)
- [Ren 2020 thesis — eScholarship](https://escholarship.org/uc/item/4j36q4hs)
- [Habaš 1990 — TU Wien PDF](https://www.iue.tuwien.ac.at/pdf/ib_1990/JB1990_Habas_2.pdf)
- [A non-GCA model for ground-plane MOSFETs (SSE)](https://www.sciencedirect.com/science/article/abs/pii/S0038110123001673) · [Non-GCA modeling of near-threshold I–V of DG MOSFETs (SSE)](https://www.sciencedirect.com/science/article/abs/pii/S0038110119306781) — abstracts not retrievable here (E3)

---

## Addendum — status at the close of C-3 (2026-10-30)

The matrix above is unchanged. New evidence since it was written:
- **K1:** T1 on 72 points: 71 PASS per path, 1 class I, 0 FAIL (`C3_T1_RIa_RIb_Report_2026-10-30.md`). Still verification infrastructure.
- **K4:** the lateral-refinement runs it lists as "pending" are done (design §1.1).
  - After level-wise extrapolation, the same-mesh and continuum metrics agree to 0.1 % (V_G = 1.25 V) and 1.2 % (2.0 V) at L = 1 µm.
  - At 0.8 V they do not agree (unresolved).
  - The gain from the same-mesh metric is therefore real where the signal is resolved, but it is not demonstrated in subthreshold.
- **K3:** the reservoir-symmetric part σ is mostly numerical at levels 0–2. Its extrapolated limit is small: ≈ 2.8e-5 at L = 1 µm and ≈ 1.5e-7 at L = 10 µm (V_G = 1.25 V). The antisymmetric part α is mesh-converged and dominant. This strengthens the refutation path listed for K3; it does not decide it (R-I(c) and ≥ 4 levels are needed).
- **Literature:** L-1 (Ward 1990 and thesis ch. 5–6), L-2, L-9 and L-10 are still open. No evidence level changed.

**Status remains `NOVELTY_UNPROVEN`.**
