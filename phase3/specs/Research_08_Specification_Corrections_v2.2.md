# Research 08 — Specification Corrections (R04/R05 v2.1 → v2.2)

**Date:** 2026-10-30. **Status:** correction note. The historical documents (Research_04 v2.1, Research_05 v2.1, Research_07) are **not modified**: they remain in the Project as issued. This note lists every changed passage with its old and new wording and the reason. Where this note and v2.1 differ, this note governs from 2026-10-30. Each correction is classed as *clarification* (no change to what is computed), *correction* (the old text was wrong or self-contradictory) or *new requirement*.

---

## 1. T0/T2 tolerance: factor-of-two inconsistency — *correction*

| Source | Old wording (v2.1) |
|---|---|
| R05 C-2 | "T0: abs(ψ_s,1COL,finest − ψ_s,M0) ≤ U_disc(ψ_s) + 1e-6 V … T2: abs(N_s,1COL,finest/N_s,M0 − 1) ≤ U_disc,rel + U_GCA,rel" |
| R04 §8.1 | "T0 … ≤ 2·U_disc(ψ_s) + 10⁻⁶ V … T2 … ≤ 2·(U_disc,rel + U_GCA,rel)" |

**New wording (both documents):**
- *T0 pass:* abs(ψ_s,1COL,finest − ψ_s,M0) ≤ 2·(U_disc(ψ_s) + U_iter(ψ_s) + U_GCA(ψ_s)) + 10⁻⁶ V.
- *T2 pass:* abs(N_s,1COL,finest/N_s,M0 − 1) ≤ 2·(U_disc,rel + U_iter,rel + U_GCA,rel).
- *Applicability:* both require the U_disc sequence to be classified **asymptotic** (R04 §7.4). If it is monotone non-asymptotic or oscillatory, the result is reported as **inconclusive**, not pass, regardless of the numbers.
- *Mandatory diagnostics:* the ratio abs(difference)/U_disc and the Richardson-extrapolated difference are reported. The latter is a diagnostic only.

**Reason.** The factor 2 is the same margin used by the claim rule (OC1) and the cross-code rule (R04 §7.7). A single factor across all comparison rules avoids testing consistency more strictly than claims. U_iter was missing from both versions. **Effect on C-2:** every T0/T2 point passed under both old forms and passes under the new one (classification asymptotic, U_iter ≤ 1e-16 relative).

## 2. D-1 grading across mesh levels — *correction*

| Source | Old wording (v2.1) |
|---|---|
| R05 D-1 | "Test: grading ratios are identical across levels (max deviation ≤ 1e-12), and breakpoints are nodes on every level." |
| R04 §7.3 | "Every level therefore has the same grading (geometric similarity), so h_k = h₀/r^k holds for the representative spacing" |

**New wording.**
- *Map:* each axis uses one fixed monotone map X(ξ) per segment between fixed breakpoints, identical on every level. Level k samples ξ uniformly with N_k = round(N₀·r^k) cells per segment, r = √2.
- *Tests:*
  - (i) breakpoints are nodes on every level;
  - (ii) the map is identical across levels, checked by evaluating X at common ξ;
  - (iii) monotone coordinates;
  - (iv) the actual per-segment refinement ratios N_{k+1}/N_k are recorded, together with h_rep,k/h_rep,k+1.
- *Grading:* adjacent-cell grading ratios depend on the level (Q^{1/N_k} for a geometric segment). They are **not** required to be identical across levels.
- *Uncertainty:* the observed order is computed with the **actual** h_rep ratio of each level pair, not with a nominal √2.

**Reason.** A fixed smooth map sampled more finely necessarily has adjacent ratios that approach 1. The v2.1 test was impossible to satisfy together with the fixed-map requirement. The D-1 replacement (`phase3/meshmap.py`, 36 tests) implements the new wording.

**Effect.**
- Per segment, rounding gives N_{k+1}/N_k = 1.400–1.429 instead of 1.414.
- For the whole depth axis, the actual h_rep ratios are 1.4124–1.4147.
- Checked on one point (ψ_s at V_G = 2.0 V, C-2 data): p₁₂₃ = 1.984 with nominal √2 and 1.991 with the actual ratios. The C-2 conclusions are unaffected.
- The C-3/C-4 analysis uses the actual ratios.

## 3. Numerical meaning of the T1 error target — *clarification*

**Old wording (R04 §5.2/§8.1, R05 C-3):** "Pass if abs(r̂₀) ≤ 3E and E ≤ E_req = 10⁻⁶ … E_req is a stated resolving-power design target, not a bound."

**New wording (added; the pass rule is unchanged):**
1. **E is a sum of estimates.** E_X = E_freq + E_eq + E_lin,X (L1/L2) or T + E_noise (L3). It is not a bound. The factor 3 in abs(r̂₀) ≤ 3E is a heuristic coverage margin.
2. **E_req = 10⁻⁶ fixes the resolving power.** It is the smallest inconsistency T1 is required to be able to detect. 3·E_req corresponds to a ~0.03 µV potential-reference error; see the R04 §5.2 sensitivity argument.
3. **What T1 verifies.** T1 verifies the *discrete* identity on each mesh level: the boundary-condition and linear-response implementation of the 2-D deck (D1–D4).
   - It does **not** verify mesh convergence. The convergence of G_h toward M0 is reported separately with R04 §7.4.
   - It does **not** verify the GCA at finite V_D.
   - It says **nothing about novelty**.
4. **E components are pre-registered** in the driver (`phase3/c3/t1.py`, module docstring) before execution:
   - *E_lin,L1* = max over three cross-sections of abs(G_sec − G_L1)/G, plus abs(G_L1 − G_NS)/G (cross-solver);
   - *E_lin,L2* = max abs(Re I_sec − Re Y)/G, plus abs(Re Y − G_L1)/G;
   - *E_eq* = abs(Δr̂₀) under τ_eq ÷ 100;
   - *E_freq* = abs(ReY(1 Hz) − ReY(0.1 Hz))/(0.99 G).
5. **Results stand as computed.** A result may not be re-evaluated with a different E definition after it is computed. Any later refinement of E is a new, versioned definition, applied to all configurations.

## 4. C-1(b) acceptance criteria — *correction* (replaces the R05 C-1(b) test figure)

**Old wording (R05 C-1(b)):** "one Newton step computed externally with this matrix must reproduce DEVSIM's own update to ≤ 1e-10 relative."

**New wording:** criteria B1–B6 as decided in `phase3/reports/C1b_Final_Acceptance_Review_2026-10-30.md` §2:
- B1: bitwise identity of matrix and RHS with the factored Newton system; the state is unchanged by export.
- B2: directional and column derivatives agree to finite-difference accuracy, with ≤ 1e-7 as the detection level.
- B3: equilibrated normwise backward error ≤ 1e-13.
- B4: cross-solver difference reported and compared with the *estimated* κ_Skeel·u as a consistency check, not a bound.
- B5: the linear-solve part of E_lin ≤ 1e-8 for T1.
- B6: damped steps are compared only through DEVSIM's update formulas.

The old figure is withdrawn and recorded as **not met** (1.4–1.9e-10, retained in the evidence).

## 5. Triple-point treatment for R-I(a) and R-I(b) — *new requirement*

**Old (R04 §4.2):** the lateral Si boundaries carry Dirichlet ψ and φ_n data; the lateral oxide boundaries are homogeneous Neumann; the interface carries continuity. Silent on the point where all three meet.

**New (added to R04 §4.2 and §4.3):**
- **R-I(a).** At each Si/SiO₂/lateral-boundary point the oxide node carries ψ = (the Si lateral data at that point). This is the discrete continuity condition with a Dirichlet twin. No other oxide node may carry a Dirichlet condition.
- **R-I(b).** The Si lateral boundary carries no ψ equation, so the interface continuity row must be active at the triple point. Verify this from the exported matrix.
- **Required checks (both models):**
  - potential jump at every coincident oxide/Si pair = 0;
  - discrete Gauss closure at round-off;
  - no single-entry oxide rows other than the gate and (R-I(a) only) the two triple points;
  - lateral uniformity of the V_D = 0 solution at round-off.

**Reason.** DEVSIM's mesher and interface assembly drop the interface pair at a node that is also a contact node for the same equation (Corner_Interface_Defect_Report_2026-10-30). Without the rule above the T1 identity fails by up to 94 %.

## 6. One-column reference formulation (R05 C-2) — *clarification*

**Old:** "electron continuity in Si with SG flux; … left/right Si edges: φ_n Dirichlet (= φ_n,slice)".

**New (added):** the reduced formulation is admissible — Poisson with a spatially uniform φ_n, no lateral contacts. Admissibility requires that equivalence with the full zero-current, frozen-hole formulation is **demonstrated under the actual boundary conditions** (equal lateral φ_n, zero normal electron flux at the interface and bottom, no generation or recombination). The demonstration must show all of the following:
- continuity residual ≈ 0;
- SG fluxes ≈ 0;
- the quasi-Fermi relations hold;
- ψ_s and N_s agree with the full formulation at round-off.

It is not a statement about the 2-D problem at V_D ≠ 0. C-2 provides this demonstration.

## 7. Reservoir sensitivity wording — *clarification*

S_res(a,b) is reported as the sensitivity between **two specified reservoir models** (R04 §4.3, Research_07 Issue 3). It is never described as a bound on all admissible reservoir conditions. This applies to every C-4/E-2 output.

---

## Erratum E1 (2026-10-30, C-3) — observed order with non-constant refinement ratios

§2 above states: "Checked on one point (ψ_s at V_G = 2.0 V, C-2 data): p₁₂₃ = 1.984 with nominal √2 and 1.991 with the actual ratios." That sentence is left unchanged as issued. **The value 1.991 is wrong.** It was computed as ln(d₂₃/d₁₂)/ln r₁₂ with the actual r₁₂. That formula is exact only when the refinement ratio is constant.

**Corrected rule (added to §2):**
- With actual ratios r₁₂ = h₂/h₁ and r₂₃ = h₃/h₂, the observed order solves p = [ln(d₂₃/d₁₂) + ln((r₁₂ᵖ − 1)/(r₂₃ᵖ − 1))]/ln r₁₂ (Celik et al., *J. Fluids Eng.* 130, 078001, 2008, monotone case).
- It is solved by fixed-point iteration.
- Implemented in `phase3/c3/t1_analysis.py::observed_order` (commit `9a88e17`) and unit-tested on exact power-law sequences (`phase3/tests/test_c3_analysis.py`).

**Corrected check value:**

| Point | Nominal √2 | v1 formula | Corrected |
|---|---|---|---|
| ψ_s, V_G = 2.0 V | 1.984 | 1.991 | **1.981** |

The C-2 conclusions remain unaffected: asymptotic classification and T0/T2 pass. The C-3 G-convergence results with the corrected rule are in `C3_T1_RIa_RIb_Report_2026-10-30.md` §4.1.
