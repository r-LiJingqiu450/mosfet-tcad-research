# WP3 pre-registration: resolved-well tier A-R and held-out states

**Date:** 2026-11-04 · **Author:** Claude · **Branch:** `benchmark/wp3-prereg-2026-11-04`, from the WP1 branch @ `e91e2b0`
**Status: FROZEN on commit.** No A-R or held-out state below has been simulated. Nothing in this document may be changed after the first WP2 run except by a dated, reasoned amendment appended at the end; the original text stays.
**Supporting estimates:** `benchmark/specs/support/` (`body_geometry.py`, `hole_geometry.py` and their outputs). These are planning estimates on simplified geometry, not device simulations.

---

## 1. Decisions in one paragraph

All A-R and held-out runs use a **standard ohmic body contact (configuration B1)** as the primary, externally scored configuration. Every state is also run with the **electron-blocking body used by all existing decks (configuration B0)** on the same mesh, as an internal companion. With B1 the device is a three-terminal electron network. The benchmark therefore reports three zero-bias conductances: drain self-conductance **G_DD**, source–drain transfer conductance **G_SD**, and drain–body leakage **G_DB**. Each has a cancellation-free reference form.

The scored current is the **electron** current only; holes are not part of the model.

Doping becomes fully analytic, with erfc grading laterally and vertically, so that no doping step falls on a mesh line.

---

## 2. Body contact: final configuration

### 2.1 Why it must be fixed explicitly

With no generation/recombination, the only drain–body path is conduction of minority electrons through the p substrate. Its conductivity is σ_b = q μ_n n_i²/N_A:
- 2.16e-14 S/cm at N_A = 1e18 cm⁻³;
- 2.16e-12 S/cm at N_A = 1e16 cm⁻³.

A Laplace solve on the actual well/body geometry (`support/body_geometry.py`) gives these factors in units of σ_b:

| Device | Grounding the body (B1 vs B0) raises G_DD by | and lowers G_SD by |
|---|---|---|
| dev1 | 0.60–0.66 | 0.39–0.40 |
| dev2 | 0.87–2.48 | 0.38–0.44 |

The ranges come from two variants: wells alone, and wells enlarged by the junction depletion width. Relative to the channel conductance, this gives the effects below. The channel G₀ is taken from the matching deck-W states as a proxy.

| State (ψ̂ = (ψ_s,∞ − ψ_b0)/V_T) | σ_b / G₀ | Estimated ΔG_DD/G₀ (B1 − B0) | Estimated ΔG_SD/G₀ |
|---|---:|---:|---:|
| dev1 (N_A 1e18, t_ox 2 nm, L 1 µm), ψ̂ = 12 | 6.05e-3 | +3.6e-3 … +4.0e-3 | −2.3e-3 … −2.4e-3 |
| dev1, ψ̂ = 16 | 1.31e-4 | +7.9e-5 … +8.7e-5 | −5.1e-5 … −5.2e-5 |
| dev1, ψ̂ = 21 | 1.00e-6 | +6.1e-7 … +6.6e-7 | −3.9e-7 … −4.0e-7 |
| dev2 (N_A 1e16, t_ox 20 nm, L 2 µm), ψ̂ = 12 | 8.21e-4 | +7.1e-4 … +2.0e-3 | −3.1e-4 … −3.6e-4 |
| dev2, ψ̂ = 16 | 1.74e-5 | +1.5e-5 … +4.3e-5 | −6.6e-6 … −7.6e-6 |
| dev2, ψ̂ = 21 | 1.31e-7 | +1.1e-7 … +3.2e-7 | −5.0e-8 … −5.7e-8 |

The A-W tolerances in WP1 are 7.8e-5 to 3.1e-4 relative. **At ψ̂ = 12 the body choice changes the answer by 10–50× the tolerance**, and at ψ̂ = 16 it is comparable to the tolerance. A benchmark that left the body boundary condition open would therefore be ill-posed.

### 2.2 B1: primary configuration (external scoring)

The boundary is the Si bottom y = 1.2 µm, over the full width −l_w ≤ x ≤ L + l_w:
- ψ = ψ_b0;
- φ_n = 0, i.e. n = n_i exp(ψ_b0/V_T) = n_i²/p₀;
- holes at equilibrium.

This is a standard ohmic contact and grounded like the source. The substrate depth of 1.2 µm is part of the device definition, because G_DB depends on it.

### 2.3 B0: companion configuration (internal; links to tiers A-W, B-W and C)

The boundary is the same Si bottom:
- ψ = ψ_b0;
- **zero normal electron flux.** No electron contact equation is imposed, so electrons are blocked.

These are exactly the boundary data of every existing deck in this repository. In a code whose contacts are always ohmic, B0 corresponds to a contact with zero electron surface recombination velocity and ψ fixed at ψ_b0.

B0 keeps the two-terminal structure, so G_SD = G_DD = G₀ and the same-mesh Rayleigh bracket remains available. That bracket is the only same-mesh bound we have.

### 2.4 Why B1 is primary

1. It is the default body contact in every drift–diffusion code, so external users need no special boundary condition. This serves H3 (external usefulness).
2. Every B1 terminal conductance has an exact, cancellation-free reference form (§3), so moving to three terminals costs no accuracy.
3. Running B0 on the same mesh keeps continuity with all existing tiers and measures the body effect directly. That measurement also tests the planning model in §2.1 (hypothesis H-B).

### 2.5 Holes: excluded from the observable

The model does not solve holes. A code that solves them adds a hole path through each n⁺ well. Its estimate is σ_p,well = q μ_p n_i²/N_D = 7.2e-17 S/cm times a geometric factor of 5.2, i.e. about 3.8e-16 S/cm per well (`support/hole_geometry.py`). Relative to G₀ this is:
- 1.05e-4 at dev1, ψ̂ = 12, which is **above** the A-W tolerance;
- 2.3e-6 at dev1, ψ̂ = 16;
- 1.4e-7 at dev2, ψ̂ = 12.

**Rule:** submissions report the **electron component** of each terminal current, or disable holes and hold them at equilibrium, p = n_i exp(−ψ/V_T). The body is ohmic for holes in both B0 and B1. This matters only to codes that solve holes, and those must still report electron current.

---

## 3. Observables and reference forms (three-terminal zero-bias network)

At equilibrium ∇φ_n = 0. The O(V_D) electron current therefore solves ∇·(σ₀∇δφ_n) = 0 with σ₀ = q μ_n n₀ (C6 theory, P4); its discrete form is the Scharfetter–Gummel resistor network with Laplacian 𝓛. Excite with φ = 1 at the drain and φ = 0 at the source and body (B1). The terminal currents are then I_D, I_S and I_B, with I_D + I_S + I_B = 0.

| Observable | Definition | Reference form (cancellation-free) | Second form (round-off check) |
|---|---|---|---|
| **G_DD** (primary) | ∂I_D/∂V_D | edge energy φᵀ𝓛φ (exact, because every fixed node except the drain has φ = 0) | −(I_S + I_B), both evaluated where φ = 0 |
| **G_SD** (primary) | −∂I_S/∂V_D | −I_S | reciprocity: excite the source with 1 (others 0) and take −I_D |
| G_DB | −∂I_B/∂V_D | −I_B | G_DD − G_SD |
| ΔL (derived) | (L/2)(1 − G_GCA/G_SD) | — | — |

The three-terminal energy identity was checked on the planning network: energy equals I_D to the printed precision (`support/body_geometry_output.txt`).

In B0 there is no body terminal, so G_DD = G_SD = G₀, and the WP1 bracket rule applies.

---

## 4. Device definition for A-R and the held-out states

The deck is resolved-well deck R as in `phase3/c7/endprobe.py`, **except where marked ★**.

- **Geometry:**
  - silicon −l_w ≤ x ≤ L + l_w with l_w = 0.3 µm, depth 0 ≤ y ≤ 1.2 µm;
  - SiO₂ thickness t_ox over the whole width;
  - gate on the oxide top over 0 ≤ x ≤ L (no overlap);
  - oxide top outside the gate is homogeneous Neumann.
- **Doping ★:** N_A uniform, plus donors N_D(x, y) = 1e20 · f_x(x) · f_y(y) cm⁻³, where
  - f_x = ½ erfc((x − x_c)/(√2 σ)) + ½ erfc(((L − x) − x_c)/(√2 σ)), with x_c = −√2 σ erfcinv(2N_A/N_D);
  - f_y = ½ erfc((y − h_w)/(√2 σ_y));
  - σ = σ_y = 15 nm, h_w = 0.5 µm.

  *Reason for ★:* the existing deck R uses a box step at y = h_w, and with box integration the node on the step owns half a cell on each side. That makes the effective well depth mesh-dependent. The vertical erfc removes the ambiguity for external codes.
- **Source and drain contacts:** ohmic, on the far edges x = −l_w and x = L + l_w for 0 ≤ y ≤ h_w:
  - ψ = V_T asinh(N_net/(2n_i)) + (applied bias, drain only);
  - φ_n = 0 (source) or V_D (drain).

  Below h_w the far edges are homogeneous Neumann. The oxide triple points at the far edges are tied to the adjacent contact ψ.
- **Body:** B1 or B0, as in §2.
- **Physics:**
  - Boltzmann statistics, constant μ_n = 1350 cm²/(V s), no generation/recombination;
  - T = 300 K, n_i = 1e10 cm⁻³, ε_Si = 11.7 ε₀, ε_ox = 3.9 ε₀, E_c − E_i = 0.56 eV;
  - ψ is referenced to the intrinsic level.
- **Gate:** ψ_gate = V_G + ψ_b0, i.e. V_FB = 0. In work-function terms, the gate work function equals the p-substrate work function, so bands are flat at V_G = 0. With these constants that is φ_m = χ_Si + 0.56 V − ψ_b0, which gives:
  - χ + 1.0362 V for N_A = 1e18 cm⁻³;
  - χ + 0.9172 V for N_A = 1e16 cm⁻³.
- **Meshes (network and DD paths):** the endprobe deck-R maps at levels 0, 1 and 2, with refinement factor √2 per level:
  - lateral `rprobe.lateral_nodes(L, level, 0.3)`;
  - depth `endprobe.depth_map(N_A)`, to 1.2 µm;
  - oxide `oxide_map(t_ox, 5)`.

## 5. Code paths and independence rules

| Path | What it computes | Independence |
|---|---|---|
| **P1 network** (new `benchmark/wp2/rdeck_probe.py`) | DEVSIM equilibrium; independent network assembly; G_DD, G_SD, G_DB by both forms; the B0 bracket; the DEVSIM coupled small-signal `G_E`/`G_M` with κ_G (level 1) | Shares the DEVSIM equilibrium with P2 and P3 |
| **P2 DD small-V_D** (new `benchmark/wp2/rdeck_dd.py`) | Full drift–diffusion at V_D = ±0.1 mV and ±1 mV; symmetric differences; DEVSIM default contact currents and the weighting-function current Σ_e J_e (w_i − w_j), with w from P1 | Same mesh as P1; different solution (finite bias) |
| **P3 coupled small-signal** | Inside P1, as in WP1 | Same discretization; ill-conditioned by design (atlas) |
| **P4 independent FV-R** (new `benchmark/wp2/fv_rdeck.py`) | Poisson–Boltzmann equilibrium and SG network on its own tensor mesh, for B1 and B0 | **May not import** from `phase3/`, `semiconductor_sim/` or DEVSIM. Re-declares its constants. Builds its own mesh (minimum spacing ≤ 0.5 nm at the interface and at x = 0, L; maximum ≤ 10 nm; written into its output). Its own Newton solve. |

Before any production run, each new script must pass a smoke test on the planning geometry. With uniform conductivity, P1's network and P4 must reproduce the `support/body_geometry.py` factors to 1e-3 relative at matched grid spacing. A failure is fixed and reported; the fix is not tuned against A-R results.

## 6. States

**A-R (6 states).** Each is run with B1 and B0, at levels 0, 1 and 2, plus P4.

| Id | Device | L (µm) | V_G (V) | Target ψ̂ |
|---|---|---:|---:|---:|
| AR-1 | dev1: N_A = 1e18, t_ox = 2 nm | 1 | 0.488174 | 12 |
| AR-2 | dev1 | 1 | 0.621433 | 16 |
| AR-3 | dev1 | 1 | 0.7835 | 21 |
| AR-4 | dev2: N_A = 1e16, t_ox = 20 nm | 2 | 0.488174 | 12 |
| AR-5 | dev2 | 2 | 0.621433 | 16 |
| AR-6 | dev2 | 2 | 0.7835 | 21 |

The V_G → ψ̂ mapping comes from the long-channel 1-D MOS problem, which is independent of the source/drain deck. These are the C10 values.

**Held-out (5 states).** Each is run with B1 and B0 at level 1 only, with P1 (including P3) and P2. These states are not used to set any threshold or tolerance.

| Id | Device | L (µm) | V_G (V) | ψ̂ |
|---|---|---:|---:|---:|
| HO-1 | dev1 | 1 | 0.588534 | 15 |
| HO-2 | dev1 | 1 | 0.686557 | 18 |
| HO-3 | dev2 | 2 | 0.588534 | 15 |
| HO-4 | dev2 | 2 | 0.686557 | 18 |
| HO-5 | dev1 | **2** | 0.488174 | 12 |

## 7. Reference values and tolerance (A-R, B1; applied separately to G_DD and G_SD)

- **Reference:** the P1 value at level 2.
- **u_disc** = max(|Q_L2/Q_L1 − 1|, |Q_P4/Q_L2 − 1|).
- **u_round** = |reference form / second form − 1| (§3).
- **u_solver** = the network relative residual.
- **Tolerance** = 2 (u_disc + u_round + u_solver), relative. ΔL tolerance follows as in WP1.

The observed convergence order from levels 0, 1 and 2 is reported as a diagnostic only; it does not enter the tolerance.

**Same-mesh scoring for B1** (no bracket exists): a same-mesh submission passes if it agrees with the reference within **τ_num = 1e-6** relative. *Basis:*
- the reference forms agree to ≤ 5.7e-11 (WP1, A-W);
- DEVSIM's coupled small-signal agrees with the network to 2.5–3.3e-7 in well-conditioned states (WP1, well-depth h_w = 1 µm and the C7 duplicate).

The value 1e-6 therefore passes legitimate same-mesh paths and fails ill-conditioned ones. B0 companions keep the bracket rule.

## 8. Hypotheses, thresholds and kill criteria

**H1: reference quality (A-R, B1).** A state enters tier A-R only if all of the following hold, for both G_DD and G_SD:

| Check | Condition |
|---|---|
| (i) cross-path | \|Q_P4/Q_L2 − 1\| ≤ 5e-3 |
| (ii) refinement | \|Q_L2/Q_L1 − 1\| ≤ 5e-3 |
| (iii) self-consistency at every level | reference vs second form ≤ 1e-9; KCL closure \|G_DD − G_SD − G_DB\|/G_DD ≤ 1e-9; the B0 companion source-side G₀ lies inside its own-mesh bracket |

**Kill:** fewer than 4 of 6 states pass. The A-R tier is then not released, the benchmark ships with internal tiers only, and the failure is reported.

**H-B: body effect (test of §2.1; never a kill).** Let Δ_DD = G_DD(B1)/G₀(B0) − 1 and Δ_SD = G_SD(B1)/G₀(B0) − 1, at level 2. The pre-registered bands are the §2.1 estimates widened by ×/÷3, to cover the deck-W→R proxy:

| State | Δ_DD band | Δ_SD band |
|---|---|---|
| AR-1 | [1.2e-3, 1.2e-2] | [−7.2e-3, −7.8e-4] |
| AR-2 | [2.6e-5, 2.6e-4] | [−1.6e-4, −1.7e-5] |
| AR-3 | [2.0e-7, 2.0e-6] | [−1.2e-6, −1.3e-7] |
| AR-4 | [2.4e-4, 6.1e-3] | [−1.1e-3, −1.0e-4] |
| AR-5 | [5.0e-6, 1.3e-4] | [−2.3e-5, −2.2e-6] |
| AR-6 | [3.8e-8, 9.7e-7] | [−1.7e-7, −1.7e-8] |

The planning model is supported if at least 5 of 6 states fall in band for each Δ, with the predicted signs. Otherwise the leakage model is revised and the miss is reported. Either way, the measured Δ values are published with the cases.

**H2: reference robustness on held-out states (kill for release).** All five states must meet all three conditions:
- (a) G_DD energy vs −(I_S + I_B) ≤ 1e-9;
- (b) KCL closure ≤ 1e-9;
- (c) the B0 source-side G₀ lies inside its own-mesh bracket.

**Kill:** any failure. The reference definition is then not robust; release stops and the failure is investigated.

**P-DD: full drift–diffusion resolution floor (atlas prediction; falsifiable; never a kill).**

The basis is the absolute current-error floor of DEVSIM's density-variable drift–diffusion on this deck family, as measured in the R-II finite-bias runs:
- terminal KCL residual 1e-13 to 6e-12 A/cm;
- cross-section spread 6e-11 to 1.0e-10 A/cm;
- both independent of the current level from 1e-5 to 1e-1 A/cm (WP1 cases `C__R2fb_*`).

**Prediction:** at V_D = ±1 mV, G_DD obtained from DEVSIM's default drain contact current misses the P1 reference (same level) by **> 1e-2** whenever G₀ · 1 mV ≤ 1e-10 A/cm. It is **within 1e-2** whenever G₀ · 1 mV ≥ 1e-8 A/cm. No prediction is made in between.

| Prediction | States |
|---|---|
| Predicted > 1e-2 | AR-1, AR-2, AR-3, AR-4, HO-1, HO-2, HO-3, HO-5 |
| Predicted ≤ 1e-2 | AR-6 |
| No prediction | AR-5, HO-4 |

The weighting-function variant and the ±0.1 mV results are measured and reported without a prediction.

**P-SS: coupled small-signal (atlas prediction; never a kill).** For dev1 states AR-1, AR-2, AR-3, HO-1, HO-2 and HO-5 (B1, level 1), `G_E` misses G_DD by > 1e-2. *Basis:* the A-W errors were 358 (ψ̂ = 12), 7.5 (ψ̂ = 16) and 4.6e-2 (ψ̂ = 21), and conditioning worsens as ψ̂ decreases. dev2 states are reported without a prediction.

## 9. Amendments to plan v2 (made before any run)

1. **H1(iv) removed from the reference criteria.** Plan v2 required DD at 0.1 mV to agree with the network within 0.5%. The measured DD floor (§8, P-DD: 6e-12 to 1e-10 A/cm) is at or above the 1 mV signal G₀ · 1 mV in 6 to 9 of the 11 planned states, depending on which end of the floor applies. This condition would therefore make H1 untestable rather than informative. It becomes the atlas prediction P-DD.
2. **No same-mesh bracket in B1.** Shorting the columns ties every column to the grounded body, so the "slice" bound collapses. B1 uses τ_num instead; B0 keeps the bracket.
3. **Doping made fully analytic**, with a vertical erfc (§4).
4. **Holes excluded** through the electron-current rule (§2.5).

## 10. Run budget, order and integrity

**Budget:**
- P1: 6 states × 2 body configurations × 3 levels = 36 runs; 10 held-out runs.
- P2: B1 at levels 1 and 2 plus B0 at level 1 for A-R, which is 18 runs; B1 at level 1 for held-out, which is 5 runs.
- P4: 12 runs.

That is about 81 runs in total. Level-0 runs take 15–100 s; level 2 is expected to take a few minutes.

**Order:**
1. Freeze this document (commit and push).
2. Write the WP2 scripts and pass their smoke tests (§5).
3. Run A-R.
4. Run held-out.
5. Evaluate H1, H-B and H2, and record P-DD and P-SS.
6. Export schema 1.1 cases and update the scorer: three-terminal observables, τ_num for B1, bracket only for B0.

**Integrity:**
- No threshold, state, level or tolerance rule may change after results are seen, except by appended amendment, and the original stays.
- Every run is logged with command, exit code and environment.
- Failed runs are kept.
- Production code, golden files, `main` and existing results remain untouched.

---

## Amendments (appended 2026-11-04 18:40 MST, before any production A-R or held-out run)

All evidence below comes from **non-production states only** (dev1/dev2 at V_G = 0.55 V and 1.0 V, labelled TRY/SMOKE in the logs). The original text above is unchanged. No threshold, state, level or kill criterion is changed.

**A1 — Smoke-test reference (§5).** The literal criterion **failed**: P1 and P4 differ from `support/body_geometry.py` by −2.5e-3 to −7.2e-3, against a limit of 1e-3. The cause is a difference in boundary weights. The planning script gives edges on the outer Neumann boundary unit weight, while box integration (DEVSIM and P4) gives them half weight. Against the same planning algorithm with box weights (`benchmark/wp2/smoke_geometry.py: planning_fv`), P1 and P4 agree to ≤ 9.2e-13 and agree with each other to ≤ 1.2e-12. **No change was made to P1 or P4.** For H-B, the planning factors move by < 0.75%, which is negligible against the ×/÷3 bands, so the bands are not changed. Evidence: `benchmark/wp2/logs/smoke_geometry.log` and `results/smoke_geometry.json`.

**A2 — P3 implementation (§5).** `phase3/c1/corner_study.eqnums` factors the raw, unscaled Newton matrix only to record the equation numbering, and on a dev2 level-1 trial state SuperLU reported that matrix "exactly singular". P3 therefore uses `rdeck_common.eqnums_capture`, which records the numbering without factoring. The G_E / G_M / κ_G computation is otherwise identical to `phase3/c6/wprobe.devsim_G`. Both routes gave identical results on a trial state (G_E relative −3.39e-4, κ_G = 4.0e13). Cost at level 1 on the largest geometry: 138 s and 3.0 GB.

**A3 — P4 production mesh (§5).** The mesh rules are unchanged: ≤ 0.5 nm at the interface and at x = 0, L, and ≤ 10 nm everywhere. Production P4 (`--mesh refined`) uses:
- first cell 0.125 nm laterally and 0.025 nm in depth;
- growth 1.04 laterally and 1.025 in depth;
- cap 5 nm, with 1.25 nm in the well-bottom band.

Every state is also run with the base settings (0.25 nm, 0.05 nm, 1.08, 1.05, 10 nm, 2.5 nm) as a P4 self-convergence companion. *Basis* (trial state dev1, V_G = 1.0 V, ψ̂ = 27.9):
- refined vs base differs by −4.75e-4;
- P4 with a 0.25 nm interface cell moves +1.7%, onto P1 level 0;
- P1 level 0 differs from base P4 by −1.19%.

P1's interface cells (0.25 / 0.177 / 0.125 nm at levels 0 / 1 / 2) are therefore expected to dominate P1's discretization error. This is recorded as an expectation; the H1 thresholds stay as registered.

**A4 — Contact extent (§4).** P1 and P2 inherit the existing deck-R convention: the far-edge contacts end on the depth line nearest h_w, which lies at 0.498259–0.503363 µm depending on device and level. P4 places the contact at exactly 0.5 µm. `benchmark/wp2/contact_extent_check.py` (P4, non-production states) bounds the effect at ≤ 6.5e-9 relative for G_DD and G_SD and ≤ 1.5e-11 for G_DB. The device definition stays at 0.5 µm.

**A5 — P2 details and the symmetric-difference term (§5, §8).**
- P2 Newton settings: τ = 1e-12 V on both ψ and φ_n, with step halving on failure; s_couple = 1.
- The symmetric difference carries a physical O(v²) term, with v = V_D/V_T. In both smoke states G_sym/G_net − 1 = 2.49e-4 at 1 mV and 2.49e-6 at 0.1 mV, i.e. v²/6.
- P-DD is evaluated exactly as registered (±1 mV, threshold 1e-2 ≫ 2.5e-4).
- A Richardson value G_R = (100·G(0.1 mV) − G(1 mV))/99 is reported as a diagnostic only.

**A6 — Exploratory observation and an informed expectation (not a pre-registered prediction).** The smoke states showed two things.
- (a) DEVSIM's post-solve contact currents (`get_contact_current`) at both drain and source carry an absolute error of about 1e-11 A/cm at level 0. At dev1, V_G = 0.55 V the drain current is off by +7.8 relative at 1 mV and the source current has the wrong sign. This is consistent with Newton updates of n ≈ 1e20 cm⁻³ in the wells at machine precision, which leave flux residuals of about 1e-13 A/cm per well node.
- (b) The circuit current −V1.I and the weighting-function current agree with each other, and their Richardson values agree with the same-run network G_DD to ≤ 2.2e-10.

**Expectation E-CIRC** (informed by these smoke states, so weak evidence): in every production P2 run, the Richardson value from the circuit current agrees with the same-level network G_DD within 1e-6. It is recorded so that the production runs test it, and it is not used for any release decision.

**A7 — Execution.** `benchmark/wp2/run_queue.py` runs the 81 registered jobs with two workers under a 7 GB memory budget and single-threaded BLAS. Every run's command, exit code, wall time and peak memory is appended to `benchmark/wp2/results/ledger.jsonl`, with per-run logs in `benchmark/wp2/logs/production/`.

---

## Amendment A8 (appended 2026-11-05, during production; changes no state, level, threshold, tolerance rule or kill criterion)

**A8 — Execution correction after an out-of-memory kill.** Some production results existed when this was written; nothing below depends on them.

- **What happened.** At 01:59 UTC the run `P2_AR-4_B1_lev2` was killed by the kernel (exit −9 after 630 s, peak RSS 5.57 GB). The kernel report shows a **memory-cgroup limit of 5.83 GiB** for every process started from the session shell, not the 8 GB of the machine that A7 assumed. A post-hoc diagnostic (`benchmark/wp2/forms_diagnostic.py`, below 1 GB) that I had started outside the queue was running at the same time. Whether the job alone would have exceeded the limit is not known.
- **Session 1 stopped cleanly.** The queue process was paused (SIGSTOP) by a watcher the moment the job running at that time (`P1_AR-5_B0_lev2`) wrote its ledger line, and was then terminated. No job was interrupted.
- **Corrections to `run_queue.py` for session 2 on:**
  - memory budget 5.4 GB (was 7.0);
  - the dev2 level-2 P2 runs are estimated at 5.9 GB, and a job above the budget may run when nothing else is running, so these run alone;
  - the log of every earlier attempt is kept as `<job>.attemptN.log`;
  - no other process is started while the queue runs.
- The failed attempt stays in the ledger, and its log is kept (`P2_AR-4_B1_lev2.attempt1.log`).
- **Contingency, fixed now.** If a dev2 level-2 P2 run exceeds the limit even when running alone, it is reported as not executable in this environment. Those runs feed only the atlas prediction P-DD and the expectation E-CIRC (never a kill), which then rest on the level-1 P2 runs of the same states. This will be stated with the results.
