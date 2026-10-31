# Ward 1988 / Ward–Odeh–Cohen 1990 — Full-Text Audit (C-4 Task E)

**Date:** 2026-10-30. **Status: NOT PERFORMED. The full texts were not accessible in this environment.** No equation, section content or assumption of these works is reported from memory or inferred from titles.

## 1. Access attempts (this session)

| Source | Route | Outcome |
|---|---|---|
| Ward, Odeh & Cohen, *SIAM J. Appl. Math.* 50(4):1099–1125, 1990, DOI 10.1137/0150066 | CaltechAUTHORS record 31400 → file `WARsiamjam90.pdf` (2.8 MB) | The record page was read: title, authors, year, journal and DOI confirmed. The file URL redirects to a 60-second signed S3 URL, which the fetch proxy rejects ("url exceeds the maximum fetchable length"). The shell has no route (CONNECT 403). Not retried further |
| Ward, PhD thesis, Caltech 1988, "Asymptotic Methods in Semiconductor Device Modeling", DOI 10.7907/08TK-NM84 | `thesis.library.caltech.edu/440/3/Ward_jw_1988.pdf` (10 MB scan with text layer) | The tool's text extraction is **truncated at 89,891 characters**, ending mid-sentence in §3.1 (p. 51). Chapters 5–6, including §5.2 (linear/saturation regimes) and §6.4 ("The I-V Curves in the Linear Regime Including End Effects"), are outside the extracted text. Paging with an offset beyond 89,891 returns nothing |
| Alternative copies | Web search for the 1990 title | Only the publisher page and the Caltech records; no other copy |

**Readable from the thesis (pages i–51 only), at evidence level E1 for those pages:**
- the table of contents;
- the Abstract and §1.2 statement that the I–V curve "in a particular bias regime, including source and drain edge effects" is obtained from a uniform expression for the mobile charge.

That is a statement of scope, not a result.

## 2. Consequence for the questions of Task E

| Question | Answer from Ward's full text | What the project can say without it |
|---|---|---|
| Small-drain-bias discrepancy scaling | **Unknown** | The project's own derivation (C4_Small_Bias_Perturbation_Analysis): for mirror-symmetric ideal reservoirs the GCA current is exact through O(V_D²), so ε ∝ V_D². The argument is elementary, so a δ² law **cannot be claimed as new** whether or not Ward states it |
| Higher-order aspect-ratio corrections | **Unknown** (Ward thesis ch. 1 introduces ε = aspect ratio; E1 for ch. 1 from an earlier session) | Scaling argument (Prop. 3): the interior and end contributions are both O(W²/L²) |
| Source/drain boundary-layer contributions | **Unknown** (§6.3–6.4 titles only, E3) | Numerical: the reservoir-dependent part is ∝ L⁻² (C-3), and the end-type dependence is strong (TC-1, ℓ_res diagnostic) |
| Dependence of current corrections on end conditions | **Unknown** | Numerical evidence in C-4 (ℓ_res diagnostic, TC-1) |
| Length scaling of the correction coefficients | **Unknown** | Heuristic L⁻² (Prop. 3); α ∝ L^−1.95 observed (C-3) |

**Required action:** attach the two PDFs. The audit then has to record, for every relevant result:
- section;
- equation number;
- scaling;
- assumptions (doping profile, n⁺ regions vs ideal reservoirs, mobility, quasi-Fermi treatment);
- relation to K1–K6 and Propositions 1–3.

Until then the claims are judged on the project's own analysis, which already downgrades several of them (novelty assessment, C-4).
