# C-6 Research Phase Report — The Zero-Bias Source/Drain End Correction of the Long-Channel MOSFET

**Date:** 2026-10-31. **Branch:** `phase3/C6-endcorrection-2026-10-31`, from C-5 `2922f86`. This is an exploratory phase under the C-6 autonomy grant. All runs were limited zero-bias probes; no production physics or golden data were touched.

**Evidence labels:**
- **[V]** verified numerical result;
- **[P]** proved;
- **[S]** argument sketch;
- **[H]** working hypothesis;
- **[I]** interpretation.

Literature levels E1/E2/E3 are as in the C-5 audit.

## 1. Most important findings

1. **[P][V] The zero-bias GCA error is exactly an equilibrium-electrostatics problem.**
   - At V_D → 0 the Scharfetter–Gummel flux has ∂J/∂ψ = 0. The linear-response conductance therefore equals that of a resistor network whose edge conductances come from the equilibrium electron density (theory note, P4).
   - An independent network assembly reproduces DEVSIM's full coupled-Jacobian linear response at all 28 states where both were computed. The relative difference is ≤ 4.8e-6, always within DEVSIM's own adjoint error estimate for that state (2e-9 to 5.8e-6).
   - Hence Ward's §6.4 problem splits cleanly into two parts: (i) the equilibrium charge, and (ii) how current flows through it.

2. **[P][V] Two-sided rigorous bounds bracket the 2-D result:** G_rows ≤ G_2D ≤ G_slice, by Rayleigh monotonicity.
   - G_slice is Ward's quasi-1-D current hypothesis (6.33) evaluated with the **exact** 2-D charge.
   - The bounds hold at every state.

3. **[V] Affine length law.** 1/G(L) = (L − 2ΔL)/σ_∞ to within 1e-5 of ΔL for L ≥ 2 µm, comparing lengths at a fixed mesh level.

   | V_G (V) | ΔL per end (nm), L = 2 / 4 / 8 µm | Change from L = 1 → 2 µm |
   |---|---|---|
   | 1.0 | 369.864 / 369.873 / 369.874 | — |
   | 1.25 | 153.686 / 153.687 / 153.687 | 0.55 nm (end-layer overlap) |
   | 2.0 | 53.385 / 53.385 / 53.385 | — |

   The remainder after the ΔL term therefore decays beyond all algebraic orders, as P6 predicts. Ward's ε-ordering would in general allow an O(ε²) term.

   **Consequence:** ΔL(V_G) is a well-defined, L-independent, intrinsic length of the end structure.

4. **[V] The intrinsic end correction is large and strongly gate-dependent.** Abrupt n⁺ junction (10²⁰ cm⁻³), no overlap, constant mobility, F0 physics; W_dep = 304 nm. ΔL per end at L = 4 µm, mesh level 0:

   | V_G (V) | N_s (cm⁻²) | ΔL_2D (nm) | 1 − G_GCA/G_2D at L = 10 µm | L\* for 1 % error = 2ΔL/η |
   |---|---|---|---|---|
   | 1.0 (threshold) | 6.2e9 | **369.9** | 7.4 % | **74 µm** |
   | 1.25 | 1.6e11 | **153.7** | 3.1 % | 31 µm |
   | 1.5 | 3.9e11 | 92.9 | 1.9 % | 19 µm |
   | 2.0 | 9.0e11 | 53.4 | 1.1 % | 11 µm |
   | 3.0 | 1.9e12 | 29.9 | 0.6 % | 6 µm |

   Ward's quasi-1-D validity rule ε ≤ 0.01 means L ≥ 15 µm for these parameters ([W88] p. 16; [WOC90] p. 1104). It is **non-conservative by up to 5× near threshold** and conservative in strong inversion. Even at L = 8 µm, G_2D exceeds G_GCA by 10.2 % at threshold.

5. **[V] First quantitative test of Ward's §6.4 theory.** The a-priori predictions were fixed in C-5 before any 2-D run. The ratio ρ = ΔL_2D/ΔL_Ward is:

   | V_G (V) | 1.0 | 1.25 | 1.5 | 2.0 | 3.0 |
   |---|---|---|---|---|---|
   | ρ | 2.59 | 1.57 | 1.31 | 1.15 | 1.07 |

   Against the C-5 decision bands, the theory is **REFUTED at threshold, DEFICIENT at 1.25 V, and nominally CONFIRMED from 1.5 V up.** (C-5 left the 1.25 V and 1.5 V points of D-c only proposed; they are reported, not re-graded.)

6. **[V] Exact error budget of Ward's reduced model.** By construction,

   ΔL_2D = ΔL_Ward + [profile] + [charge model] − [spreading]

   where:
   - **profile** = ΔL_hybrid − ΔL_Ward, with ΔL_hybrid = Ward's own charge formula fed with the exact 2-D surface potential;
   - **charge model** = ΔL_slice − ΔL_hybrid;
   - **spreading** = ΔL_slice − ΔL_2D ≥ 0.

   All values in nm, at L = 4 µm:

   | V_G (V) | Ward | + profile | + charge model | − spreading | = 2-D |
   |---|---|---|---|---|---|
   | 1.0 | 142.7 | +226.7 | +0.9 | −0.4 | 369.9 |
   | 1.25 | 98.1 | +49.8 | +10.4 | −4.7 | 153.7 |
   | 1.5 | 70.9 | +14.1 | +15.1 | −7.2 | 92.9 |
   | 2.0 | 46.6 | −1.3 | +17.4 | −9.3 | 53.4 |
   | 3.0 | 27.9 | −4.8 | +17.2 | −10.3 | 29.9 |

   **Near threshold**, Ward's failure is almost entirely his phenomenological surface-potential profile. Given the true ψ_s(x), his corner-charge construction is accurate to 0.2 %: 369.4 vs 370.3 nm.

   **In strong inversion**, the apparent agreement (ρ = 1.07–1.15) is **partly a cancellation**:
   - the charge-model term (+17 nm; Ward's inversion-layer charge formula given ψ_s) is about 35–60 % of ΔL_Ward;
   - the neglected current spreading (−9 to −10 nm) is about 20–37 % of ΔL_Ward.

   Ward's quasi-1-D current hypothesis (6.33) overestimates ΔL by 0.1 % at threshold, 3 % at 1.25 V and 35 % at 3.0 V (spreading/ΔL_2D).

7. **[V][I] The missing physics is a depletion-scale surface-potential tail.**
   - The 2-D ψ_s(x) excess has a 2–6 nm layer next to the n⁺ edge in inversion (21 nm at threshold), plus a tail with a local decay length of 110–170 nm. That length is comparable to the quasi-2-D natural length ≈ 135 nm.
   - The tail amplitude falls from ~130 mV at threshold to ~1 mV at 3 V.
   - Ward assumes a pure oxide-fringe profile with decay length 2t_ox/π = 12.7 nm ([W88] p. 121).
   - [I] The tail is the long-channel limit of charge sharing. Near threshold the inversion charge responds exponentially to it.

8. **[V] Ward's domain truncation is accurate. The result is insensitive to well extent and weakly sensitive to well depth.**

   | Test | Effect on ΔL |
   |---|---|
   | Resolved abrupt n⁺ wells (R-II) vs Ward's equipotential-well boundary, same mesh family, L = 2 µm | +0.10 % (370.23 vs 369.86 nm, V_G = 1.0) and −0.17 % (1.25 V), both within mesh uncertainty (§5a); −1.9 % (2.0 V), resolved |
   | Doubling the well's lateral extent (0.3 → 0.6 µm) | −0.01 nm |
   | Doubling the well depth (0.5 → 1.0 µm) | +1.5–1.8 % |
   | Mesh levels 0 → 1 → 2 | ≤ 0.17 nm (≤ 5e-4 relative). Not asymptotic (§5a) |

   The well-depth sensitivity is the one Ward could only estimate in the section affected by the missing p. 142.

9. **[H rejected] Equivalent threshold-shift form (P7) in its simple version.**
   - ΔL·N_s is not constant in V_G.
   - A two-term law ΔL = ΔL_g + A/N_s fitted at two gate voltages mispredicts the other two by −13 % and +39 %.
   - The end term is algebraic (1/L) in L at every V_G, but its V_G dependence is not a one-parameter law.

## 2. Research paths considered

| Path | Status | Reason |
|---|---|---|
| D1: GCA error against idealized reservoirs | Rejected (C-5) | Physically subleading; reference-dependent |
| H-R / H-S / FB-1 | Not pursued | Same reason |
| C-W as registered (E-W1, finite matrix with R-II) | **Reformulated before running** | P4 showed only zero bias is needed and the network gives G cheaply, so a smaller design suffices: Ward's truncation first, then resolved wells as an ablation |
| Ward-truncation probe W1 (25 completed states + 1 development state) | **Done** | Findings 1–6, 9 |
| Resolved-well probe R1 (4 states + 1 development state) | **Done** | Finding 8 |
| Ward-hybrid and slice decomposition (analysis only) | **Done** | Finding 6. This is what turned "Ward is ~15 % off in strong inversion" into "Ward is right for partly wrong reasons" |
| P7 threshold-shift equivalence | **Tested, rejected** | Finding 9 |
| Local-ψ_s "fringe vs field" split | Reported only as a diagnostic | Ill-conditioned: an exponential sensitivity to the far-field reference made it L-dependent until it was referenced to the 2-D channel centre |
| Quasi-2-D tail as the basis for an improved reduced model | **Emerging, not developed** | §3 and §7 |
| Realistic junctions (overlap, grading, Fermi–Dirac, μ(N, E)) | Not done | Decides physical significance (§6) |

## 3. Strongest remaining candidate contribution

**C-E: "Source/drain end correction of the long-channel MOSFET conductance: an exact zero-bias decomposition and the first quantitative test of the asymptotic end-effect theory."** It has four parts:

1. **Exact framework [P].**
   - Zero-bias conductance = equilibrium-weighted conduction problem (P4).
   - Rigorous bounds (P5).
   - Affine length law with beyond-all-orders remainder (P6, [S] with numerical confirmation to 1e-5).
2. **Intrinsic ΔL(V_G) of an abrupt junction [V].**
   - L-independent, about 1.2 W_dep per end at threshold, falling to 0.1 W_dep at 3 V.
   - The linear-regime GCA validity length is L\* = 2ΔL/η, which replaces Ward's ε ≤ 0.01 rule.
3. **Test of Ward's §6.4 against pre-registered predictions [V].**
   - Refuted near threshold; nominally confirmed in strong inversion.
   - An exact error budget shows *why*:
     - the oxide-fringe profile misses a depletion-scale charge-sharing tail;
     - in strong inversion two errors of opposite sign partly cancel.
4. **[H] An improved reduced model** that keeps Ward's corner-charge construction (accurate to 0.2 % given ψ_s near threshold) and replaces his profile with a quasi-2-D charge-sharing tail. This is the natural next deliverable and is falsifiable against the verified 2-D data already in hand.

**Novelty status:**

| Part | Status | Basis |
|---|---|---|
| (3) Ward test and error budget | **`NOVELTY_POTENTIAL`** | Ward's §6.4 was never evaluated or published as a MOSFET result ([WOC90] p. 1108 "under investigation"; [W88] p. 138 "left to future investigators"; E1). [W92] (E2) covers buried-channel FETs; [RW91]/[WRO94] are E3 |
| (2) Intrinsic gate-dependent ΔL near threshold | **`NOVELTY_UNPROVEN`** | High prior-art risk. Channel-length extraction and its 2-D interpretation are crowded: Taur 2000 T-ED (E2); Narayanan et al. 1995 SSE (E2: 2-D numerical study, junction-based L definition "questionable"); the charge-sharing and quasi-2-D threshold literature (E3) |
| (1) Exact framework | **Not novel as mathematics** | Linear-response decoupling and Rayleigh bounds are classical. Its use as an exact error decomposition of a published reduced device theory is the methodological contribution |
| (4) Improved model | Not yet developed | — |

**Significance (candid):**
- The device class (long-channel, abrupt, no overlap, 20 nm oxide, 10¹⁶ cm⁻³) is classical.
- The result matters to anyone using long-channel GCA/charge-sheet conductance near threshold. Examples: mobility and threshold extraction on "long" test structures, and compact-model end corrections.
- Whether ΔL(V_G) survives realistic junctions decides relevance for T-ED (§6).
- Today's best venue estimate: *Solid-State Electronics* / *IEEE J-EDS* / *J. Comput. Electron.* for the device angle, or *SIAM J. Appl. Math.* for the asymptotic-verification angle. T-ED only if §6 succeeds and the improved model (4) works.

## 4. Mathematical arguments

`phase3/reports/C6_Theory_ZeroBias_EndCorrection_2026-10-31.md` contains:
- P4 (proved, SG-discrete and continuum);
- P5 (proved);
- P6 (sketch, with two stated gaps);
- P7 (hypothesis, rejected in simple form);
- the end-layer observation.

The error-budget identity (finding 6) is exact by construction. Its terms are defined operationally:

| Quantity | Definition |
|---|---|
| ΔL_slice | Network with every column shorted |
| ΔL_hybrid | Ward's p. 139 formula with v_s ← ψ_s,2D/(V_T lnλ) and α₀ₗ(ỹ) from §6.3 |
| ΔL_Ward | Ward's p. 139 formula with his p. 121 profile |

## 5. Numerical evidence (reproducible)

**Code:**

| File | Purpose |
|---|---|
| `phase3/c6/wprobe.py` | W1 deck, homotopy equilibrium, network G, bounds, DEVSIM cross-check |
| `phase3/c6/rprobe.py` | R1, resolved wells |
| `phase3/c6/w1_analysis.py` | Tables, Ward-hybrid, profiles |
| `phase3/c6/w1_figures.py` | Figure |
| `phase3/c5/ward64_quadrature.py` | Ward a priori. C-6 robustness fix for ỹ → 0 (ellipkm1); the C-5 predictions are reproduced to 1e-14 |

**Tests:** `phase3/tests/test_c6_network.py` covers network G against an analytic sheet and Rayleigh bounds on random 12-decade networks.

**Data:**
- `results/phase3/C6/w1/`: every JSON, logs, `RUNS.txt`, `RUNNER_END.txt`, `CANCELLED.txt`.
- `results/phase3/C6/analysis_w1/`: `w1_summary.json`, `w1_table.csv`, `C6_W1_figure.png`.

**Failures and cancellations (kept, not hidden):**
- **Out-of-memory kills** (7 GB cgroup; `dmesg` evidence in `CANCELLED.txt`):
  - L = 8 µm, level 0, V_G = 1.25, run concurrently. Re-run alone: OK.
  - L = 2 µm, level 2, V_G = 1.25, killed during the DEVSIM matrix-export cross-check. Level 2 was re-run with `--no-devsim` (network G only).
- **Runner B:** L = 8 µm level-1 jobs cancelled. The affine law is already established without them.
- **Runner D:** remaining level-2 jobs cancelled and re-queued as runner E. Runner D's earlier jobs were not run twice.

**Numerical quality:**

| Check | Result |
|---|---|
| Equilibrium | Newton updates ≤ 1e-16 V at τ = 1e-12 V |
| Network linear solves | Residual ≤ 1e-15 |
| Network vs DEVSIM | ≤ 4.8e-6, within the DEVSIM adjoint estimate at every state (28 states) |
| Mesh levels 0 → 1 → 2 | ΔL changes by ≤ 0.17 nm. Not asymptotic; working U ≈ 0.5 nm (§5a) |

### 5a. Level-2 convergence

| V_G (V) | ΔL_2D, levels 0 / 1 / 2 (nm) | Successive changes (nm) | Snapped well depth h_w, levels 0 / 1 / 2 (nm) |
|---|---|---|---|
| 1.0 | 369.864 / 369.827 / 369.656 | −0.037, −0.171 | 501.5 / 503.3 / 501.5 |
| 1.25 | 153.686 / 153.669 / 153.590 | −0.017, −0.079 | same |
| 2.0 | 53.385 / 53.358 / 53.320 | −0.027, −0.037 | same |

L = 2 µm. Level 2 was run with `--no-devsim`.

**The sequence is not in the asymptotic range.** The level 1 → 2 change exceeds the level 0 → 1 change. One identified confounder: the well depth is snapped to the nearest depth line, and that line differs at level 1 (+1.8 nm, against an observed sensitivity of +5.4 nm/µm). No observed order and no extrapolation are claimed.

**Conservative working uncertainty: U ≈ 0.5 nm per end** (≤ 0.2 % of ΔL). This is about 3× the largest observed change, an estimate, not a bound. The cause should be removed in the next probe by making h_w a fixed depth breakpoint.

Every conclusion in §1 involves differences ≥ 1 % except two: the resolved-well vs truncation differences at 1.0 V and 1.25 V (0.36 nm and 0.26 nm). Those two are **within U and therefore unresolved**; only the 2.0 V difference (1.0 nm, 1.9 %) is resolved. The affine-law statement (finding 3) compares lengths at a fixed mesh level and is a same-discretization result.

## 6. Unresolved scientific uncertainties

1. **Physical realism.**
   - Gate overlap (0 here) and lateral junction grading (abrupt here) will change the near-threshold tail, possibly strongly.
   - Fermi–Dirac statistics in the 10²⁰ well and μ(N, E_⊥) are absent. Constant μ means ΔL is a *charge-weighted*, not a *mobility-weighted*, quantity.
   - This decides whether finding 4 is physics or an artefact of the idealization.
2. **Prior art for ΔL(V_G) near threshold.** These full texts are not yet read:
   - Taur 2000;
   - Narayanan et al. 1995;
   - Terada–Muta 1979;
   - Chern et al. 1980;
   - Ward 1992;
   - Reyna–Ward 1991;
   - Liu et al. 1993;
   - Yau 1974.
3. **P6 rigour.** The two gaps stated in the theory note remain.
4. **Interpretation of the charge-model term.** The +17 nm in strong inversion is attributed [I] to the junction-adjacent bulk/spillover charge that Ward's inversion-layer formula excludes. Not demonstrated.
5. **Finite bias.** Everything here is V_D → 0. The O(V_D) behaviour of the end correction (Ward's drain-side formula) is untested.
6. **2-D only.** Width effects are excluded by construction.

## 7. Recommendation for the next investment (in order, each gated)

1. **Literature gate LG-1 (no computation).**
   - Read the eight works in §6.2.
   - Kill the "intrinsic ΔL(V_G)" claim if it is already quantified for abrupt junctions against 2-D drift-diffusion.
   - Several of these are paywalled. The user would need to supply PDFs (list in §6.2).
2. **Physical-significance gate (small, about 20–30 zero-bias states, about 2 CPU-h).** Repeat the V_G scan at L = 4 µm with:
   - gate overlap 0 / 10 / 30 nm;
   - a lateral Gaussian junction, σ = 10 and 30 nm;
   - Fermi–Dirac statistics in the well, if DEVSIM support is adequate.

   **Abandon C-E as device physics** if ΔL at threshold changes sign, or drops below 0.1 W_dep, for realistic overlap and grading. In that case keep only part (3), the Ward test, as an applied-math/verification paper.
3. **Improved reduced model (analysis only).**
   - Ward's corner charge + a quasi-2-D tail with natural length ℓ, giving a closed form for ΔL(V_G).
   - Pre-register its predictions against the existing 2-D data at the five V_G before fitting anything.
   - **Abandon** if it misses by more than 15 % at any V_G without refitting.
4. **Only then** a manuscript decision and venue choice.

**Not recommended:** large sweeps, finite-bias campaigns, or modern short-channel devices (ε ~ 1, no theory to test).
