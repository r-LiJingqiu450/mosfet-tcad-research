# Low-current MOSFET drift–diffusion benchmark

A code-agnostic reference set for the **zero-bias (linear-response) terminal conductances of 2-D MOSFETs at low current**. It comes with documented failure signatures of common ways of evaluating that current. It complements the open drift–diffusion discretization benchmark of Farrell & Peschka (2019), which covers forward-biased p-i-n diodes. No claim of methodological novelty is made.

**Releases**
- **WP1 (2026-11-04)** exported existing project data, with no new simulation: internal tiers A-W and B-W, plus stress tier C. Schema 1.0.0; these 72 case files are unchanged since.
- **WP2 (2026-11-04 MST; runs 2026-11-05 01:32–04:59 UTC)** adds new production runs for the externally usable resolved-well deck, in schema 1.1.0:
  - **A-R**: three-terminal reference with an established tolerance;
  - **HO**: held-out states;
  - **B-R**: electron-blocking-body companions.

  All of it was pre-registered in `specs/WP3_Preregistration_AR_and_HeldOut_2026-11-04.md` before the first production run. Amendments A1–A9 are appended there with dates; none changes a state, threshold or tolerance rule. Results and evaluation: `reports/WP2_Report_2026-11-04.md`.
- **WP4 (2026-11-05; runs 06:17–11:02 UTC)** scores DEVSIM's own user-facing output paths as submissions. No case file changes. The paths are:
  - the default post-solve contact current in double precision (H3): fails 17 of 17 scorable cases;
  - the same current with DEVSIM's extended precision: passes 6 of 6 A-R cases on the level-1 mesh;
  - DEVSIM's built-in AC analysis, in double and in extended precision: fails 6 of 6 each.

  Pre-registered in `specs/WP4_Preregistration_External_Paths_2026-11-04.md`, with amendments W4-A1 to W4-A3 (execution and release form only). Results: `reports/WP4_Report_2026-11-05.md`. NiTCAD was not built, because Eigen 5.0.1 is unavailable; its evaluation is on hold.
- **WP6 (2026-11-05; runs 11:53–18:06 UTC)** is an AI-agent study, not a benchmark change; no case file changes.
  - Fresh headless Claude Code sessions on Claude Opus 5.5 computed G_DD for the 11 A-R and HO three-terminal states with the kit's DEVSIM deck. Arm A had no checker: 33 runs, 3 per task.
  - **Result: 0 silent errors among 33 valid runs** (95% upper bound 0.106). 32 runs reported a value within 1.1e-4 of the reference; one ran out of time.
  - Every run called DEVSIM's default contact current and 30 ran its AC analysis. All set them aside in favour of bias-step differences of the circuit current or the conductance network, each with an independent check.
  - The registered stop rule therefore applied, and arm B (with the scorer as checker) was not run.
  - Pre-registered in `specs/WP6_Preregistration_Agent_Study_2026-11-05.md` (amendments W6-A1, W6-A2). Results, sensitivity and limits: `reports/WP6_Report_Agent_Study_2026-11-05.md`.
- **WP7 (2026-11-05)** prepares the public release and the manuscript. No case file changes.
  - `paper/`: a LaTeX manuscript whose every figure, table and quoted result is generated from the released files.
  - `benchmark/tools/build_release.py`: builds a benchmark-only release snapshot (MIT for code, CC BY 4.0 for data and text).
  - Report: `reports/WP7_Report_Release_and_Manuscript_2026-11-05.md`.

| Tier | Cases | Content | Scored by |
|---|---:|---|---|
| **A-R** | 6 | Resolved-well deck R, ohmic body (B1), three-terminal G_DD / G_SD / G_DB at ψ̂ ≈ 12, 16, 21 for two devices. The reference is P1 at level 2; all 6 passed the pre-registered H1. | Tolerance for any mesh or code; τ_num = 1e-6 on a published mesh |
| HO | 5 | Held-out deck-R states, level 1 only, never used to set a threshold. All 5 passed H2. | τ_num on the published level-1 mesh only |
| B-R | 23 | Electron-blocking-body (B0) companions of every A-R and HO state, one case per mesh | Same-mesh bracket only |
| A-W | 3 | Ward-truncated deck W (internal) | Tolerance + same-mesh bracket |
| B-W | 28 | Deck-W deep-subthreshold states (internal) | Same-mesh bracket only |
| C | 41 | Deck-R states at V_G ≥ 1.0 V, well-depth states and finite-bias stress groups | Bracket / tolerance where established / not scorable |

## Scope envelope (every case)

| Setting | Value |
|---|---|
| Statistics | Boltzmann |
| Mobility | constant, μ_n = 1350 cm²/(V s) |
| Holes | not part of the observable. Report the **electron** current, or hold p = n_i exp(−ψ/V_T). |
| Generation/recombination | none |
| Temperature and intrinsic density | T = 300 K, n_i = 1e10 cm⁻³ |
| Permittivities | ε_Si = 11.7 ε₀, ε_ox = 3.9 ε₀ |
| Gate | ψ_gate = V_G + ψ_b0 (V_FB = 0); ψ is referenced to the intrinsic level |
| Dimensionality | 2-D: currents in A/cm, conductances in S/cm per unit width |

Each case file contains the full device specification under `device`: geometry, doping, every boundary condition and the gate work function. A case can therefore be set up without reading the project's code.

## Observables and reference method

At equilibrium ∇φ_n = 0, so the O(V_D) electron current is that of the medium σ₀ = q μ_n n₀(x, y). It does not depend on the electrostatic response (`phase3/reports/C6_Theory_ZeroBias_EndCorrection_2026-10-31.md`, P4). Every code can therefore be scored on the V_D → 0 limit of its terminal currents.

The reference values are terminal currents of the discrete Scharfetter–Gummel network of this linear-response problem, evaluated in cancellation-free form. This is the established weighting-function / characteristic-potential evaluation of contact current:
- E. Palm and F. Van de Wiele, *IEEE T-ED* ED-32 (1985) 2052–2059;
- H. Nah et al., IEEE 2003, pp. 75–78;
- P. Farrell and D. Peschka, *CAMWA* 78 (2019) 3731–3747.

The loss of accuracy of direct flux sums at low current in heavily doped regions has been known since Polsky & Rimshans (1981), as cited by Palm & Van de Wiele.

**Three-terminal cases** (A-R, HO). The body is ohmic (B1): ψ = ψ_b0 and φ_n = 0 at y = 1.2 µm. The drain is excited with φ_n = 1, and source and body are held at 0.

| Observable | Definition | Reference form | Scored |
|---|---|---|---|
| G_DD | ∂I_D/∂V_D | edge energy φᵀ𝓛φ (second form −(I_S + I_B)) | yes |
| G_SD | −∂I_S/∂V_D | −I_S (second form: reciprocity) | yes |
| G_DB | −∂I_B/∂V_D | −I_B | reported only |
| ΔL | (L/2)(1 − G_GCA/G_SD), nm per end | — | yes (tolerance) |

**Two-terminal cases** (B-R, A-W, B-W, C) have G₀ = lim I_D/V_D with the electron-blocking body (B0), and ΔL = (L/2)(1 − G_GCA/G₀).

## Resolved-well device (A-R, HO, B-R)

| Item | Specification |
|---|---|
| Silicon | −l_w ≤ x ≤ L + l_w, with l_w = 0.3 µm; depth 0 ≤ y ≤ 1.2 µm |
| Oxide and gate | SiO₂ of thickness t_ox; gate over 0 ≤ x ≤ L |
| Doping | N_A uniform plus analytic n⁺ wells N_D = 1e20 f_x(x) f_y(y) cm⁻³: lateral erfc with the metallurgical junction at x = 0 and L (σ = 15 nm); vertical erfc at h_w = 0.5 µm (σ_y = 15 nm) |
| Source/drain | ohmic, on the far edges for 0 ≤ y ≤ h_w |
| dev1 | N_A = 1e18 cm⁻³, t_ox = 2 nm, L = 1 µm; gate work function χ + 1.0362 V |
| dev2 | N_A = 1e16 cm⁻³, t_ox = 20 nm, L = 2 µm; gate work function χ + 0.9172 V |

The body choice matters at the level of the tolerances: at ψ̂ = 12, B1 raises G_DD by 0.4 % (dev1) relative to B0.

**Meshes.** `meshes/mesh_<device>_L<L>_lev<n>.json` gives the exact node lines of every reference mesh, together with the triangulation, contact-extent and gate rules. A box-integration Scharfetter–Gummel code can rebuild the identical network and declare `same_mesh`.

## Scoring rules (`tools/score.py`)

1. **Sign.** A conductance ≤ 0 fails.
2. **Three-terminal, same mesh.** Declare `same_mesh: true` and `mesh_id`, which must be a key of `observables.same_mesh_references` (default: the case mesh). G_DD and G_SD must agree with the reference-path value on that mesh within **τ_num = 1e-6** relative. G_DB is compared and reported but does not decide.
3. **Three-terminal, any mesh** (A-R only). Two checks apply:
   - |G/G_ref − 1| ≤ tolerance, with tolerance = 2 (u_disc + u_round + u_solver), where u_disc = max(|P1_L2/P1_L1 − 1|, |P4/P1_L2 − 1|) and P4 is an independent finite-volume solver;
   - |ΔL − ΔL_ref| ≤ the ΔL tolerance.

   A same-mesh submission on a coarser published mesh is scored by rule 2 only. Its distance from the level-2 reference is reported.
4. **Two-terminal, same mesh.** The Rayleigh bracket G_rows (1 − 1e-12) ≤ G ≤ G_slice (1 + 1e-12). It bounds **this discrete network on this mesh** and is not a bound on the continuum value. B-R cases also accept G_DD or G_SD as G₀.
5. **Two-terminal tolerance** (A-W, C where established). This is the WP1 rule, k = 2 times (discretization + round-off + solver).
6. If no check applies, the case is `not_scorable`. A case passes when every scored check passes.

Because each cross-path spread is part of its own tolerance, the independent solvers pass by construction. That is a consistency check, not independent validation.

| Tier | G_DD tolerance (relative) | ΔL tolerance (nm) |
|---|---|---|
| A-R dev1 (ψ̂ = 12 / 16 / 21) | 2.8e-3 / 3.7e-3 / 4.7e-3 | 1.33 / 1.68 / 2.14 |
| A-R dev2 (ψ̂ = 12 / 16 / 21) | 2.6e-4 / 2.5e-4 / 2.4e-4 | 0.166 / 0.154 / 0.143 |

## Submitting

```json
{"submission": {"name": "my code", "current_evaluation": "weighting function"},
 "results": {
   "AR__R_L1_NA1e18_tox2_VG0.488174_B1": {"G_DD_S_per_cm": 3.53e-12, "G_SD_S_per_cm": 3.50e-12, "same_mesh": false},
   "HO__R_L1_NA1e18_tox2_VG0.588534_B1": {"G_DD_S_per_cm": 6.18e-11, "same_mesh": true, "mesh_id": "dev1_L1_lev1"}}}
```

```bash
python benchmark/tools/score.py my_submission.json --out my_report.json
```

**Worked examples** (`examples/`) re-submit the project's own output paths and, from WP4, DEVSIM's user-facing ones:

| Example | Result |
|---|---|
| Three-terminal reference values (A-R, HO), same mesh | 11/11 pass |
| DEVSIM default drain contact current (`get_contact_current`), Richardson value, same mesh | 11/11 fail |
| DEVSIM circuit-node current (−x[V1.I]), Richardson value, same mesh | 11/11 pass |
| DEVSIM coupled small-signal `G_E`, level-1 mesh | 11/11 fail |
| Independent finite volume (`wp2/fv_rdeck.py`), other mesh | 6/6 pass |
| Two-terminal source-side reference (A-W, B-W, C, B-R), same mesh | 86/86 pass |
| Two-terminal drain-side flux sum, same mesh | 38 fail, 48 pass |
| WP4: DEVSIM drain and source contact currents with extended precision, Richardson values, level-1 mesh (A-R) | 6/6 pass |
| WP4: DEVSIM built-in AC analysis (`solve(type='ac')`, 1e-6 Hz), double precision, level-1 mesh (A-R) | 6/6 fail |
| WP4: the same in extended precision | 6/6 fail |

## Failure signatures in the data (atlas)

- **Direct drain-side flux sum.** On every A-R and HO mesh it misses the same-mesh reference by more than τ_num. At ψ̂ = 12 the miss is up to about 4e3 relative, with the wrong sign. The source side of the same solution is exact to round-off. The looser B0 bracket misses it on 4 of 23 B-R meshes.
- **DEVSIM's post-solve contact current** (`get_contact_current`) has an absolute error of 3e-12 to 2.7e-11 A/cm, set by round-off in the n⁺ contact regions. On the dev1 meshes the raw values at ±0.1 mV are bit-identical across gate voltages whose G_DD spans four decades. Below about 1e-11 A/cm the value carries no channel information. **Use the circuit-node current or a weighting function instead.** Their Richardson values agree with the reference to ≤ 1.2e-9 in all 23 drift–diffusion runs.

  **DEVSIM's extended precision repairs this current** (WP4), when `extended_solver`, `extended_model` and `extended_equation` are set:
  - the current then equals the circuit-node current to float64 resolution;
  - its Richardson value matches the reference to ≤ 6.2e-10 in all 6 A-R states on the level-1 mesh;
  - the cost is ×11–16 wall time, with the same Newton iteration count.
- **Coupled small-signal output** (DEVSIM full Jacobian). The backward residual is near machine precision, yet the error is up to −7e3 relative at ψ̂ = 12, with κ_G ≈ 1e17. A residual check alone accepts it.

  **DEVSIM's built-in AC analysis is this same solve** (WP4). At 1e-6 Hz, `solve(type='ac')` reproduces `G_M` bit for bit in 5 of 6 A-R states. Its error reaches 1.7e3 relative, and 3.3e3 in extended precision. Extended precision cannot help here: the AC system is factored and solved once, in double precision.

## Reproduce

```bash
python benchmark/wp2/run_queue.py --list            # the 81 registered WP2 jobs (ledger: wp2/results/ledger.jsonl)
python benchmark/wp2/evaluate_wp2.py                # registered evaluation -> wp2/results/evaluation.json
python benchmark/tools/export_cases.py              # deterministic; rewrites cases/, meshes/, INDEX.json, MANIFEST.sha256
python benchmark/tools/make_examples.py
python -m pytest -q benchmark/tests
python benchmark/wp4/run_queue_wp4.py --list        # the 18 registered WP4 jobs (ledger: wp4/results/ledger.jsonl)
python benchmark/wp4/evaluate_wp4.py                # registered WP4 evaluation -> wp4/results/evaluation_wp4.json
python benchmark/wp4/posthoc_wp4.py                 # WP4 post-hoc diagnostics -> wp4/results/posthoc_wp4.json
```

Source files are checked before export:
- WP1 sources must be byte-identical to commit `c05139b4ae037ece34257ff0382e84a91b059773`;
- WP2 sources must be byte-identical to `454db70e71660a896b23ca378d2f5c3413808204`, the commit that holds the production results.

Each case lists its sources with their SHA-256. The WP2 and WP4 runs use DEVSIM 2.11.2, NumPy 2.3.5 and SciPy 1.17.0. A DEVSIM equilibrium needs `DEVSIM_MATH_LIBS` pointing at BLAS/LAPACK. DEVSIM's extended precision needs a build that reports `extended_precision: True` in `get_parameter(name='info')`.

## Limits

- Everything here is agreement between internal code paths: the DEVSIM network, an independent finite-volume solver and drift–diffusion. **Nothing is experimental validation.**
- WP4 scores DEVSIM's own user-facing output paths, but the reference equilibrium also comes from DEVSIM. No code independent of DEVSIM has been scored; NiTCAD is on hold.
- G_DB and the HO and B-R cases have no established tolerance.
- WP6 tested one agent (Claude Opus 5.5, custom system prompt, a stated 1% requirement, a ready deck). Its null result does not transfer to other agents, prompts or tasks.
- The physics envelope is narrow: constant mobility, Boltzmann statistics, no holes or generation/recombination, and zero or small V_D.
- The licence for public release is not yet decided.
