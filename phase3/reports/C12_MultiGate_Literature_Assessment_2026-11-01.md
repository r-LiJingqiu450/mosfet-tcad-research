# C-12 Literature Assessment — Multi-Gate Short-Channel Theory and the Zero-Bias End Correction

**Date:** 2026-11-01. **Branch:** `phase3/C12-multigate-literature-2026-11-01`, created from C-11 `ee2f5f3`.

**Scope:** an equation-level review of the five documents the user attached, plus lightweight checks that run entirely inside the published analytical model. No device simulation was run. Compute: under 10 CPU-seconds.

**Companion files:**
- Derivations: `phase3/reports/C12_Mathematical_Notes_2026-11-01.md`.
- References: `phase3/reports/C12_references.bib`.
- Scripts: `phase3/c12/lt_endcorr.py` and `phase3/c12/lt_vt_rolloff.py`.
- Tests: `phase3/tests/test_c12_lt_identities.py`.
- Data: `results/phase3/C12/`.

**Labels:**
- [V] verified in the source: full text read, page and equation given.
- [D] mathematical deduction made here. Proved in the Notes, and checked numerically where marked [N].
- [N] numerical check inside the Liang–Taur model. It is **not** a device simulation and **not** physical validation.
- [H] hypothesis. [S] speculation.
- [E2] abstract or tool-summarised text only. [E3] metadata only.

---

## 0. Bottom line

1. **The zero-drain-bias end correction we proposed for double-gate devices is a special case of an established model.** It follows from Liang & Taur (2004), eq. (16) with eq. (13), in the limits V_ds → 0 and L ≫ λ₁.
   - In those limits the end correction per end is ΔL = (λ₁/π)·F(b₁/V_T, s). Here b₁ is the published eq. (12) and s = πt_si/2λ₁.
   - In the thin-film limit s → 0, F tends to Ein. That limit is exactly the predictor λ_DG·Ein(B_proj/V_T) that C-11 §7.8 proposed to pre-register [D, N].
   - The nanowire (gate-all-around) analogue exists: Yu et al., 2009 [E3].
   - The concept itself is older: Fjeldly–Shur 1993 for bulk, and Liang's dissertation (2006) pp. 54–55 for the double gate.
   - **The proposed multi-gate end-correction test should not be run.** It would re-validate a 2004 model in the regime where that model is most accurate.

2. **The dissertation reproduces the 2004 paper; it does not correct its physics.** Chapters 2–3 reprint it in part.
   - It adds a general asymmetric-insulator solution (Appendix 1), compact threshold-voltage formulas (Ch. 3), a limitations analysis (Ch. 4) and hardware correlation (Ch. 5).
   - It also makes two editorial corrections. (i) It restores the charge q that is missing from eqs. (15)–(16) of the paper. (ii) It states the projection as an ε-weighted orthogonality, which is identical to the paper's "conjugate function" [V, D].
   - Its stated dominant error is the idealised potential at the source/drain edge (§4.2, §4.4). It identifies that error but does not resolve it.

3. **One reformulation survived this review (O1).** Inside the Liang–Taur model, the constant-current V_t roll-off **equals the end-correction term** −V_T·ln[L/(L − ΔL_s − ΔL_d)].
   - The remainder is of second order in the central barrier lowering, measured in kT/q [D, N].
   - Wherever that lowering is well below kT/q, the roll-off is therefore algebraic in L (≈ 1/L). It is not set by the minimum potential. Here that means L ≳ 3.5 λ₁; the lowering equals kT/q at 2.3–2.6 λ₁.
   - The compact formulas in the dissertation miss it: at L = 100 nm the minimum-potential formula (3.7) gives −0.2 mV and the Taylor formula (3.15) −3.2 mV, against −8.6 mV from the full integral [N]. Liang's own Fig. 3.14 shows the same gap against her 2-D simulation.
   - **Originality is weak.** The 1/(L − x_s − x_d) structure is Fjeldly–Shur's (bulk), and the double-gate form is a corollary of Liang–Taur. Significance is modest.

4. **Recommendation.**
   - Abandon the multi-gate end-correction test as formulated in C-11.
   - Do not build a multi-gate simulation framework now.
   - Treat the end-correction programme as closed as the main publication line, unless O1 passes a short, decisive literature gate (§8).

---

## 1. Sources and bibliographic verification

| # | Verified identity (from the PDF itself) | Notes |
|---|---|---|
| A | X. Liang, *Analytical Modeling of Short Channel Effects in Double Gate MOSFET*, Ph.D. dissertation, UC San Diego, 2006. 142 PDF pages; committee chair Y. Taur; eScholarship ID qt49f030x3. | Chapters 2–3 are "in part, a reprint" of item B (p. xix). The dissertation cites B as **vol. 51, no. 9, September 2004** (pp. xix–xx, 41, 73). **B is in fact no. 8, August 2004**, per the PDF running heads [V]. |
| B | X. Liang and Y. Taur, "A 2-D analytical solution for SCEs in DG MOSFETs," *IEEE T-ED* 51(8): 1385–1391, Aug. 2004. DOI 10.1109/TED.2004.832707. | Equations (5), (6), (9), (11)–(13), (16) are images. Read from rasterised pages 1386–1388. |
| C | D. J. Frank, Y. Taur, H.-S. P. Wong, "Generalized scale length for two-dimensional effects in MOSFET's," *IEEE EDL* 19(10): 385–387, Oct. 1998. | Eq. (1) and the double-gate form, read from rasterised p. 386. |
| D | S.-H. Oh, D. Monroe, J. M. Hergenrother, "Analytic description of short-channel effects in fully-depleted double-gate and cylindrical, surrounding-gate MOSFETs," *IEEE EDL* 21(9): 445–447, Sep. 2000. | Eqs. (1)–(6) read from rasterised pages. |
| E | Y. Taur, "Analytic solutions of charge and capacitance in symmetric and asymmetric double-gate MOSFETs," *IEEE T-ED* 48(12): 2861–2869, Dec. 2001. | One-dimensional theory only. |

The SHA-256 of each file read is listed in `results/phase3/C12/README.md`. The PDFs are not committed, for copyright reasons.

---

## 2. The established results most relevant to us (deliverable 1)

### R1. The boundary-value problem and its idealisations [V]
B, eqs. (1)–(2) on p. 1386; A, eqs. (2.1)–(2.10) on pp. 20–22.
- **Field equation.** Laplace's equation ∂ₓ(ε∂ₓψ) + ∂_y(ε∂_yψ) = 0 in silicon and insulator. This holds in subthreshold for an undoped or lightly doped body. The threshold shift from N_a is qN_a·t_si/2C_ox: about 1 mV at 10¹⁶ cm⁻³ (B footnote 1; A p. 19).
- **Source and drain.** Abrupt junctions. Dirichlet conditions ψ = E_g/2q and E_g/2q + V_ds hold across the **whole film** at y = 0 and y = L.
- **Gates.** Dirichlet ψ = V_g − Δφ.
- **Insulator gaps.** The potential in the gate–S/D insulator gaps is **interpolated linearly**. B states this is valid "when t_i/t_si < 1, e.g., 0.3" (p. 1386).

### R2. The eigenstructure [V], and its identity across papers [D]
B eqs. (7)–(8), p. 1386–1387: ε_si tan(πt_i/λ_n) = ε_i tan(nπ/2 − πt_si/2λ_n).
- Odd n gives the symmetric modes. Even n is needed for asymmetric gates; footnote 2 of B points out that FTW's equation covers only odd orders.
- λ_n/λ₁ ~ 1/n.

C eq. (1), p. 386: 0 = ε_Si tan(πt/Λ₁) + ε_I tan(πw/Λ₁). This is for bulk, with the body treated as a metal at the depletion depth w.

C p. 387 gives the symmetric double-gate form: 1 = (ε_Si/ε_I)·tan(πt/Λ₁)·tan(πw/2Λ₁), with w the silicon thickness.

D eq. (3), p. 446: κ·tan(t_Si/2λ)·tan(t_ox/λ) = 1, with decay ∝ exp(±x/λ). D attributes it to FTW. For the cylinder, D eq. (6) uses J₀/Y₀, with λ_Cyl ≈ (D_Si + 2κt_ox)/4.810 (p. 447).

**For the symmetric double gate these are one eigenvalue condition in three notations** [D]: FTW p. 387 ≡ D eq. (3) ≡ B eq. (8) with n = 1, and Λ₁ (FTW) = λ₁ (Liang–Taur) = π·λ (Oh). FTW eq. (1) is the bulk counterpart.

**Check [N].** Our solver reproduces every λ₁ the dissertation quotes:

| Case | Our value | Quoted |
|---|---|---|
| Fig. 4.2(c) | 11.77 nm | 11.7 nm |
| Fig. 4.5(a) | 33.57 nm | 33.6 nm |
| Fig. 4.5(b) | 33.94 nm | 33.9 nm |
| Fig. 4.5(c) | 39.04 nm | 39.0 nm |

### R3. The projection coefficients [V, D, N]
B eqs. (10)–(12) and g₁(x) at the foot of p. 1388; A eqs. (2.23)–(2.26), pp. 29–30.

**Conjugate function.** B calls the eigenfunctions "not orthogonal" and builds a "conjugate function set" g_n (p. 1387). The dissertation states instead that they are orthogonal with weight ε(x) (A p. 29).
- Using eq. (8), the printed g₁ equals ε(x)·u_L1(x, 0)/ε_si exactly [D].
- The numerical residual is ≤ 9×10⁻¹⁵ [N].
- The two statements are therefore the **same projection**.

**Closed form.** The closed form b₁ in eq. (12) matches the projection (11) of the boundary data:
- analytically [D];
- to 10⁻¹⁵ by direct projection and to 10⁻¹⁴ by quadrature of (11) with the printed g₁ [N].

**Link to C-11.** C-11's H-lin predictor (`phase3/c11/bpredict.py`) is this same projection, including the linear insulator interpolation, transplanted to the doped bulk [D].

### R4. The potential and current formula [V]
B eq. (13), lowest mode; B eqs. (14)–(16); A eqs. (2.27)–(2.30), pp. 33–34.

I_ds = μW·kT·[1 − exp(−qV_ds/kT)] / ∫₀ᴸ dy (∫ n_i e^{qψ/kT} dx)⁻¹  (A eq. 2.30)

- **Transport assumption.** The quasi-Fermi potential φ_n is constant across the film. A justifies this with DESSIS contours (Fig. 2.8, p. 34).
- **Units.** B's (15)–(16) omit the charge q (Q_i is written without q, and the prefactor as μW·kT/q). A's (2.29)–(2.30) are dimensionally correct [V].

### R5. Validation and limitations [V]
**Validation:**
- B: potential contours within 2 %; I–V curves, V_t roll-off and DIBL agree with 2-D ISE "with no fitting parameters" (pp. 1388–1391).
- B p. 1389: higher-order terms "contribute little to the total integral".

**Limitations (A):**
- §4.1 (pp. 75–76): the n = 2 term changes roll-off by < 5 %, and including it raises the current slightly.
- **§4.2 (pp. 76–82): the source/drain-edge simplification is the dominant error.**
  - Field penetration from the gates pulls the edge potential towards the gate potential.
  - This overestimates b₁ and c₁, hence the current and the roll-off.
  - The error grows as t_si/t_i decreases and as V_ds increases.
- §4.3 (pp. 82–84): high-k insulators — at ε_i = 35.1 and t_i = 13.5 nm the model deviates.
- §4.4 (p. 85): "most of the deviation ... is due to the simplification taken for the potential at source/drain edge."

### R6. Compact threshold-voltage formulas [V] (A Ch. 3)
- **Definition.** Constant-current V_t, normalised by W/L (eqs. 3.6, 3.10; p. 44: V_g0 = V_t0 − E_g/4q).
- **Minimum-potential method (3.7)**, ΔV_t = −2√(b₁c₁)·e^{−πL/2λ₁}: error about 35 % at L/λ₁ = 1.5 (p. 52).
- **Why it fails (pp. 54–55).** It overestimates the x-integral by 10 % (symmetric) and the **y-integral by 150 % at L/λ₁ ≈ 3.5**. This is the end-region effect, recognised qualitatively.
- **Taylor method (3.15)**, a parabolic barrier giving erf terms (p. 58). In Fig. 3.14 (p. 59) it tracks the full integral at short L, but at L = 100 nm it gives ≈ −2 mV where the full integral and DESSIS give ≈ −10 mV (read off the figure).
- **Curve fits.** λ₁ is fitted by (3.31) and erf by (3.32).

### R7. Context [V]
- **Taur 2001 (E).** The exact 1-D Poisson–Boltzmann solution with mobile charge only, ψ(x) with the log-cosine form, eqs. (1)–(5). This is the long-channel 1-D reference for an undoped double gate, and "volume inversion only occurs under subthreshold conditions" (abstract). It has no 2-D end physics.
- **FTW (C).** The amplitudes b₂₁ and c₂₁ are named but not computed (p. 386). The body depletion edge is "oversimplified", and the reader is referred to Nguyen 1984 (p. 385).
- **Oh (D).** λ is extracted from the semilog slope of Δψ(V_DS) and does not depend on transverse position (Fig. 2b, p. 446). D also critiques the parabolic approximation (p. 445).

### Post-2006 literature located but not read in full
- Taur & Ning, 2nd ed. (2009), appendix "An analytical solution for the short-channel effect in subthreshold", pp. 575–581 [E3].
- Yu, Yuan, Song, Taur, *IEEE T-ED* 56(10), 2009: the nanowire analogue [E3, cited by Mahajan et al., engrXiv 10.31224/7486, 2026 [E2]].
- Xie, Xu, Taur, "Review and critique of analytic models of MOSFET short-channel effects in subthreshold," *IEEE T-ED* 59(6): 1569–1579, 2012 [E3].
- Xie, Xu, Ren, Taur, *Semicond. Sci. Technol.* 25(3) 035012, 2010: high-k [E3].
- Xie et al., *IEEE T-ED* 2013: ultrathin SOI, with "an empirical, channel length-dependent scale length" [E3].
- Lu, Liang, Wang, Taur, NSTI 2006, compact short-channel DG model [E2].
- J. Song, UCSD dissertation, 2010: FinFET compact model, short-channel part from 2-D Poisson [E2].
- Ray & Mahapatra, *IEEE T-ED* 56(2): 260–266, 2009: undoped DG; subthreshold quantities taken at the channel midpoint [E2].
- Taur et al., *IEEE EDL* 1995, effective length from 2-D sheet resistance [E2].
- Taur, *IEEE T-ED* 2000, shift-and-ratio channel length [E2].

---

## 3. Relationship to C-6 – C-11 (deliverable 2)

| Our construct | Where it sits in the Liang–Taur (LT) framework |
|---|---|
| C-6 exact zero-bias network; Rayleigh bounds G_rows ≤ G ≤ G_slice | **LT eq. (16) is the slice (Dirichlet) upper bound**, evaluated with the LT potential. In Slotboom variables this extends to **any V_ds** for a frozen subthreshold potential: I_rows ≤ I_2D(ψ) ≤ I_LT(16) [D, Notes N4]. For t_si = 10 nm the bracket is I_rows/I_slice ≥ 0.9948 for every L from 20 to 1000 nm at V_ds = 0.05 V and 1 V [N]. **The transport approximation of LT is therefore rigorously within 0.52 % for its own potential.** This puts a bound on the φ_n-constant assumption that A Fig. 2.8 only illustrated. |
| C-10 form λ·Ein(B/V_T) | With (13) in (16), V_ds → 0, L ≫ λ₁: **ΔL = (λ₁/π)·F(β, s)**, where β = b₁/V_T, F = ∫₀^β [1 − 1/M(z)] dz/z and M = film average of e^{z cos(πx/λ₁)} [D]. Thin-film limit: F → Ein(β); at s = 10⁻³ the difference is 5×10⁻⁸ [N]. For large β, F = ln β + C(s), with C(0) = γ = 0.5772 and C decreasing in s (C(0.886) = γ − 0.127; C(1.17) = γ − 0.216) [D, N]. **The C-10 bulk "descriptive identity" is the bulk counterpart of an exact corollary of LT.** |
| C-11 DG proposal: ΔL = λ_DG·Ein(B_proj/V_T) to ±5 % | This is the s → 0 limit of LT. Against LT's own finite-thickness formula, the pure-Ein form is high by **+2.1 to +2.4 %** at t_si = 5 nm, **+3.6 to +4.2 %** at 10 nm and **+6.2 to +7.2 %** at 25 nm (SiO₂, t_i = 1.5 nm, barrier 0.3–0.5 V) [N]. The proposed test would have measured LT's accuracy, not a new law. |
| C-11 H-lin (bulk) | The LT/FTW eigen-projection method transplanted to the doped bulk [D]. Its 13–24 % failure comes from the nonlinear junction corner. FTW already flags the depletion-edge oversimplification (C p. 385); the quantification is ours. |
| C-11 "main risk: screening at the n⁺ face" | In LT this is the source/drain-edge idealisation of A §4.2. **Its effect on ΔL is only logarithmic:** d ln ΔL / d ln b₁ = 1/F ≈ 0.27–0.33 [D, N]. A 10 % error in b₁ moves ΔL by about 3 %. A rough estimate from A Fig. 4.2(c) (t_si = 5 nm), with values read by eye, gives b₁ overestimated by ≈ 10 % and ΔL by ≈ 3 % [N, estimate]. |
| Higher-order modes | Their effect on ΔL relative to the lowest mode is ≤ 0.02 % for the SiO₂ cases and 0.18–0.36 % for the thick high-k case. That case violates LT's interpolation condition anyway [N]. Notes N5 explains why (O(β⁻²)). |
| Finite L | ΔL(L)/ΔL(∞) = 0.978–0.991 at L = 3λ₁, and 1 to within 10⁻⁴ at L ≥ 5λ₁ [N]. |

**Revisions to earlier conclusions:**
- C-11 §7.6 named a published λ·Ein formula as the thing that would remove the remaining novelty. For multi-gate bodies the formula is a corollary of LT eq. (16); a difference in notation does not create originality. **That novelty claim is withdrawn.**
- **The C-11 recommendation to run a DG test is withdrawn.**
- In bulk, C-10/C-11 still stand as verified but narrow results.
- A note has been appended to the C-11 report.

---

## 4. Ideas already covered by prior art (deliverable 3)

**Classification of the zero-drain-bias end correction, as the user asked:**

| Category | Verdict |
|---|---|
| Mathematically equivalent to an established expression | **Yes, in a limit.** The C-10/C-11 form equals the long-channel, thin-film, V_ds → 0 evaluation of LT eqs. (13) + (16) [D] |
| Special case of an existing analytical model | **Yes.** LT eqs. (13) + (16) hold for any L, V_ds and t_si (B; A Ch. 2). The nanowire analogue is Yu et al. 2009 [E3]. A textbook appendix also exists (Taur & Ning 2009) [E3] |
| Distinct but already understood quantity | **Yes, as a concept:** a subthreshold electrical length shorter than L because of source/drain fringing. Sources: FS93 (L − x_s − x_d) [E1]; A pp. 54–55 (y-integral overestimated by 150 %) [V]; Taur et al. 1995, L_eff ≠ L_met from 2-D sheet resistance [E2] |
| Unresolved, significant problem | **No** |

**Also covered:**
- scale length (C, D, B eq. 8);
- eigen-projection amplitude (B eqs. 11–12);
- cylindrical Bessel modes (D eq. 6);
- the current integral (B eq. 16);
- high-k and two-layer insulators (C; A §4.3);
- exact 1-D double-gate charge (E).

Applying any of these to a "new" geometry would not be original.

---

## 5. Uncertainties that remain open (deliverable 4)

### Open within the five documents
- **U1 — boundary data at the source/drain edge** (A §4.2, §4.4; B p. 1386). It is the stated dominant error; the dissertation offers no correction.
  - [D, order of magnitude] For an n⁺ doping of 10²⁰ cm⁻³ and a gate-to-edge difference of about 0.4 V, a gate-induced depletion depth W ≈ (2ε_siΔV/qN_D)^{1/2} ≈ 2.3 nm. That is comparable with t_si/2 for t_si ≤ 5 nm, consistent with the thickness trend in A Fig. 4.2.
  - [H] It is a nonlinear free-boundary corner problem on the source/drain side: the counterpart of the body-side corner nonlinearity C-11 found in bulk.
- **U2 — the insulator-gap interpolation** when t_i/t_si ≳ 0.3, and high-k insulators (B p. 1386; A §4.3). This overlaps the fringing-induced barrier-lowering literature, which was not reviewed.
- **U3 — the compact V_t formulas** (A eqs. 3.7, 3.15) leave out the long-channel end correction. See O1.
- **U4 — regime limits.** Subthreshold only; undoped body; Boltzmann statistics; no quantum effects; abrupt S/D. No graded S/D or underlap.

**Status after 2006: unknown to us.** Xie–Xu–Taur 2010/2012/2013 and the textbook appendix are the first things to read [E3].

### Resolved here
- **The transport approximation** is rigorously bounded (≤ 0.52 % at t_si = 10 nm, for a given potential).
- **Mode truncation** does not matter for the end correction.

---

## 6. Opportunities found (deliverable 5)

### O1 — Long-channel subthreshold short-channel effect = end correction (algebraic in L) [D, N]

**Proposition [D, Notes N6].** Within LT eqs. (13) + (16):

∫₀ᴸ dy N_∞/N(y) = L − ΔL_s − ΔL_d + R,  with |R| ≤ const·L·(b₁c₁/V_T²)·e^{−πL/λ₁},

where ΔL_{s,d} = (λ₁/π)·F(b₁/V_T or c₁/V_T, s) ≈ (λ₁/π)·[ln(b₁ or c₁ / V_T) + C(s)].

**Consequences.** The constant-current roll-off of A eq. (3.10) is

ΔV_t(L) = −V_T·ln[L / (L − ΔL_s − ΔL_d)] + O(V_T·(Δψ_min/2V_T)²),  with Δψ_min = 2√(b₁c₁)·e^{−πL/2λ₁}.

- The first-order effect of the central barrier lowering is already inside the overlapping end tails. Adding the minimum-potential term (3.7) to the end term therefore **double-counts** (checked numerically).
- Wherever Δψ_min ≪ kT/q, the roll-off is algebraic (≈ −V_T·(ΔL_s + ΔL_d)/L).
- The crossover, Δψ_min = kT/q, falls at L_× = (2λ₁/π)·ln(2√(b₁c₁)/V_T) ≈ ΔL_s + ΔL_d + O(λ₁/π): roughly where the two end regions meet. Here L_× = 40 nm (2.3 λ₁) at 50 mV and 47 nm (2.6 λ₁) at 1 V.
- FTW's exponential picture (C p. 386) and A's compact formulas (3.7) and (3.15) describe only the regime L ≲ L_×. Eq. (3.15) tends to zero faster than any power of 1/L (Notes N6).

**Numbers [N]** (`lt_vt_rolloff.py`): t_si = 10 nm, t_i = 1.5 nm, λ₁ = 17.73 nm, barrier 0.38 V, roll-off referred to L = 200 nm as in A Fig. 3.14.

| L (L/λ₁), V_ds | Full eq. (16) | End term alone | Min-potential (3.7) | Taylor (3.15) | Min-pot. + end (double-counts) |
|---|---|---|---|---|---|
| 50 nm (2.8), 1 V | −41.4 mV | −48.9 | −19.6 | −48.6 | −68.4 |
| 70 nm (3.9), 1 V | −19.2 mV | **−19.3** | −3.3 | −19.0 | −22.6 |
| 100 nm (5.6), 1 V | **−8.6 mV** | **−8.6** | −0.2 | **−3.2** | −8.8 |
| 140 nm (7.9), 1 V | −3.3 mV | −3.3 | −0.0 | −0.2 | −3.3 |
| 100 nm (5.6), 50 mV | −6.8 mV | −6.8 | −0.1 | −1.8 | −6.9 |

**Remainder scaling.** The remainder (full − end term) behaves as (Δψ_min/2V_T)²:
- 50 mV drain: +1.64 mV at 50 nm, +0.034 mV at 70 nm, < 0.001 mV at 100 nm;
- 1 V drain: +7.4 mV at 50 nm, +0.12 mV at 70 nm.

**Consistency with the dissertation's own data [V, read off the figure].** In A Fig. 3.14 (p. 59) at L = 100 nm:
- the full integral and the DESSIS simulation are both ≈ −10 mV;
- the Taylor formula is ≈ −2 mV;
- the minimum-potential formula is ≈ 0.

This is what the proposition predicts. Bias and threshold definition differ slightly from ours, so the agreement is qualitative.

**Further algebraic terms [D, estimates]:**
- **The "long-channel" reference itself is shifted:** about −5.4 to −6.4 mV at 200 nm and −1 mV at 1 µm, relative to an infinitely long device [N].
- **DIBL:** an algebraic component, ≈ −1.8 mV at 100 nm relative to the 200-nm reference [N]. It comes from ΔL_d growing as ln(c₁).
- **Subthreshold slope:** S increases by the fraction ≈ (λ₁/π)·V_T·[1/A₀ + 1/(A₀ + V_ds)]/(L − ΔL_s − ΔL_d), where A₀ is the source barrier E_g/2q + Δφ − V_g (0.38 V here). That is ≈ 1 % at 100 nm.
- **[H] Robustness:** the algebraic term should survive realistic source/drain profiles, because the boundary amplitude enters only through a logarithm.

**Risks:**

| | Assessment |
|---|---|
| Originality | **High risk.** FS93 already writes the bulk subthreshold current as ∝ 1/(L − x_s − x_d), the same algebraic structure. A pp. 54–55 recognises the end-region effect qualitatively. The double-gate ΔL is a corollary of LT eq. (16). What is absent from the five documents is only the closed form, the second-order remainder statement and the crossover. Compact-model literature has not been checked at full-text level. Ray–Mahapatra's midpoint evaluation [E2] suggests the term is often left out in practice. |
| Significance | Modest. Absolute effects are ≲ 10–20 mV and matter at L/λ ≈ 3–10. It bears on extracting the scale length from roll-off data and on long-channel references. |
| Feasibility | High. |

### O2 — A rigorous inner–outer theory of the source/drain boundary data (U1) [H]
- **Idea:** a nonlinear corner region on the source/drain side, matched to the evanescent mode, giving b₁,eff. This would affect V_t roll-off and DIBL **linearly**, unlike the end correction.
- **Originality:** high risk. Candidates already covering it: Xie–Xu–Taur 2010–2013, conformal-mapping models (Kloes and co-workers), and the fringing-induced barrier-lowering and underlap literature.
- **Significance:** moderate, for thin-body, high-k gate-all-around devices.
- **Feasibility:** medium.

### O3 — Two-sided transport bounds for analytic subthreshold models [D]
- **Idea:** combine Slotboom variables with the Dirichlet/Thomson principles. Eq. (16) is an upper bound; the rows integral is a lower bound.
- **Novelty:** classical tools; a small result.
- **Use:** a supporting element of O1 or O2, not a paper on its own.

### O4 — Not from these documents (from C-11)
- **Observation:** the coupled small-signal solve becomes ill-conditioned in deep subthreshold, whereas the decoupled equilibrium network stays accurate.
- **Scope:** a numerical-methods topic, suited to IEEE TCAD or the Journal of Computational Electronics rather than T-ED.
- **Prior art:** extended-precision TCAD practice.

---

## 7. Originality, significance and feasibility (deliverable 6)

| Candidate | Originality risk | Significance | Feasibility | My estimate of reaching T-ED level [S] |
|---|---|---|---|---|
| C-11 DG end-correction test | **Prior art** (special case of LT) | Low | High | ≈ 0 |
| O1 long-channel SCE = end correction | High (FS93 structure; LT corollary) | Modest | High (≈ 1–2 CPU-h of validation) | ≈ 10–20 %; higher for a brief or letter |
| O2 boundary-data theory | High | Moderate | Medium | ≈ 10–20 % until the 2010–2013 papers are read |
| O3 bounds | High (classical) | Low | High | Component only |

---

## 8. Recommendation and next action (deliverables 7–8)

**Recommendation: abandon the multi-gate end-correction direction as formulated, and keep one conditional reformulation (O1).**
- No DEVSIM double-gate deck.
- No multi-gate simulation campaign.
- The end-correction programme is closed as the main publication line.

**Most valuable next action: a decisive literature gate.** It costs no compute. Read at equation level:
1. Xie, Xu, Taur, *IEEE T-ED* 59(6) 1569–1579 (2012);
2. the Taur & Ning (2009) appendix, pp. 575–581;
3. Yu, Yuan, Song, Taur, *IEEE T-ED* 56(10) (2009);
4. Xie, Xu, Ren, Taur, *SST* 25 035012 (2010), and Xie et al., *IEEE T-ED* (2013);
5. the short-channel formulation of BSIM-CMG and of the Lu–Liang–Wang–Taur compact model (NSTI 2006);
6. one conformal-mapping current model (Kloes and co-workers, 2012–2014).

My tools cannot open IEEE or IOP full texts, so these PDFs would need to be attached. **If they cannot be obtained, I recommend closing O1 too**, rather than spending compute on an unverified novelty claim.

**Evidence required to proceed with O1:**
- (a) none of the six gives the logarithmic long-channel end correction or the algebraic term explicitly;
- (b) the user accepts a brief- or letter-level target.

**Validation design if O1 proceeds** (pre-registered, ≈ 1–2 CPU-h):
- 2-D drift–diffusion double gate and axisymmetric gate-all-around;
- L/λ₁ = 2–30, abrupt and graded source/drain;
- test the closed-form ΔL and the coefficient of the 1/L term in V_t roll-off and S;
- independent finite-volume cross-check (`fd2d_indep` adapted).

**If the gate fails,** formally close the end-correction programme and re-select the research problem outside classical subthreshold electrostatics. I would not manufacture a gap here.

**What would change these conclusions:**
- A published closed form of the algebraic term → O1 is not new.
- Evidence that real devices suppress the term (for example through series resistance or graded junctions) → O1 is not significant.
- A Xie/Taur result that already solves U1 → O2 is not new.

---

## 9. Record

| Item | Content |
|---|---|
| `phase3/c12/lt_endcorr.py` | K1–K7, K9: identities, closed form, ΔL integrals, bounds, full-series and finite-L checks, sensitivity, edge estimate → `results/phase3/C12/lt_endcorr.json` |
| `phase3/c12/lt_vt_rolloff.py` | K8: roll-off decomposition and rows/slice ratios → `results/phase3/C12/lt_vt_rolloff.json` |
| `phase3/tests/test_c12_lt_identities.py` | 18 tests |
| Corrected before recording | My first draft of O1 added the minimum-potential term (3.7) to the end term. That double-counts the first-order barrier lowering, so it was replaced by the N6 proposition. The script still prints the sum, labelled as double-counting, as evidence |
| Withdrawn | The C-11 novelty remark (§3) and the C-11 DG recommendation |
| Estimates, not measurements | K9: source-edge b₁ error read by eye from A Fig. 4.2(c). Values read off A Fig. 3.14 |
