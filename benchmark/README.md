# Low-current MOSFET drift–diffusion benchmark (WP1 freeze)

A code-agnostic reference set for the **zero-bias (linear-response) drain conductance of 2-D MOSFETs at low current**, together with documented failure signatures of common ways of evaluating that current. It complements the open drift–diffusion discretization benchmark of Farrell & Peschka (2019), which covers forward-biased p-i-n diodes.

**Status: WP1 (2026-11-04).** This release exports only *existing* project data, with no new simulation. The external-facing tier **A-R** (resolved source/drain wells with ohmic contacts, deep subthreshold) is reserved and will be added in WP2 after pre-registration (WP3). The tiers present now are internal (A-W, B-W) or stress cases (C). No claim of methodological novelty is made.

## Scope envelope

Every case shares these physics settings:

| Setting | Value |
|---|---|
| Statistics | Boltzmann |
| Mobility | constant, μ_n = 1350 cm²/(V s) |
| Holes | not solved; p = n_i exp(−ψ/V_T) |
| Generation/recombination | none |
| Temperature and intrinsic density | T = 300 K, n_i = 1e10 cm⁻³ |
| Permittivities | ε_Si = 11.7 ε₀, ε_ox = 3.9 ε₀ |
| Gate | ψ_gate = V_G + ψ_b0 (V_FB = 0), with ψ referenced to the intrinsic level |
| Dimensionality | 2-D: currents in A/cm, conductances in S/cm per unit width |

The full device specification (geometry, doping, every boundary condition) is written into each case file under `device`, so that it can be set up without reading the project's code.

## Observable and reference method

**G₀ = lim_{V_D→0} I_D / V_D.** At equilibrium ∇φ_n = 0, so the O(V_D) current is that of the medium σ₀ = q μ_n n₀(x, y) and does not depend on the electrostatic response. See `phase3/reports/C6_Theory_ZeroBias_EndCorrection_2026-10-31.md`, P4.

The reference value is the terminal current of this linear-response problem, with φ_n = 0 at the source and 1 at the drain. It is evaluated on the **source side**, where it is free of cancellation, and it equals the edge-energy form φᵀLφ to about 1e-11. This is the established weighting-function / characteristic-potential evaluation of contact current:
- E. Palm and F. Van de Wiele, *IEEE T-ED* ED-32 (1985) 2052–2059;
- H. Nah et al., IEEE 2003 conference, pp. 75–78;
- P. Farrell and D. Peschka, *CAMWA* 78 (2019) 3731–3747.

The loss of accuracy in direct flux sums at low current in highly doped regions has been known since Polsky & Rimshans (1981), as cited by Palm & Van de Wiele.

The derived quantity **ΔL** (nm, per end) is defined as ΔL = (L/2)(1 − G_GCA/G₀), where G_GCA is the long-channel 1-column sheet conductance divided by L, given in each case.

## Tiers in this release

| Tier | Cases | Content | Scored by |
|---|---:|---|---|
| A-W | 3 | Ward-truncated deck W, N_A = 1e18 cm⁻³, t_ox = 2 nm, L = 1 µm, V_G = 0.488174 / 0.621433 / 0.7835 V. Reference cross-checked by an independent finite-volume solver (≤ 0.16 % in ΔL); a level-1 refinement exists at 0.7835 V. | tolerance (established) + same-mesh bracket |
| B-W | 28 | C10 deep-subthreshold deck-W states (N_A = 1e16 / 1e18, t_ox = 1–20 nm). In 16 of them the drain-side flux sum is documented as wrong. | same-mesh bracket only |
| C | 41 | 28 resolved-well (deck R) zero-bias states at V_G ≥ 1.0 V (2 of them with a level 0/1 pair, which gives a tolerance); 4 well-depth states; 9 finite-bias drift–diffusion stress groups (V_D ≤ 20 mV) | bracket (same mesh) / tolerance where established / not scorable |

**Deck W** puts Dirichlet ψ = ψ_bi on vertical segments at the junction faces. Standard TCAD cannot set this up, which is why the external-facing reference tier must come from deck R (WP2).

## Scoring rules (`tools/score.py`)

1. **Sign:** G₀ ≤ 0 fails.
2. **Same-mesh bracket**, applied *only* when the submission declares `same_mesh: true`. The check is G_rows (1 − 1e-12) ≤ G ≤ G_slice (1 + 1e-12). The bracket comes from Rayleigh monotonicity and bounds **this discrete network on this mesh**. It is not a bound on the continuum value. Example: the independent finite-volume ΔL lies outside the DEVSIM-mesh bracket at 2 of the 3 A-W states, yet agrees within 0.16 %.
3. **Tolerance**, applied only where `uncertainty.status = established`. The tolerance is k × (discretization + round-off + solver), with k = 2. This is a stated convention, not a statistical coverage. The discretization term is max(|ΔL_independentFV − ΔL_ref|, |ΔL_level1 − ΔL_level0|) where these exist. The G₀ tolerance follows from d ln G = 2 dΔL / (L − 2ΔL).
4. If no check applies, the case is `not_scorable`. A case passes when every applicable check passes.

Because the cross-path spread is part of the uncertainty estimate, the independent finite-volume result passes by construction. That is a consistency check, not independent validation.

## Submitting

```json
{"submission": {"name": "my code", "current_evaluation": "weighting function"},
 "results": {"<case id>": {"G0_S_per_cm": 1.2e-10, "same_mesh": false}}}
```

```bash
python benchmark/tools/score.py my_submission.json --out my_report.json
```

Worked examples (re-submissions of the project's own output paths) are in `examples/`:

| Example | Result |
|---|---|
| Source-side reference | 63/63 pass |
| Drain-side flux sum | 19 fail (16 B-W + 3 A-W) |
| Coupled small-signal G_E | 3/3 fail |
| Independent finite volume | 3/3 pass |

## Failure signatures already in the data (atlas seed)

- **Drain-side flux sum.** In 16 of 28 B-W states it lies outside the same-mesh bracket. At V_G = 0.488174 V its relative error in G is 119, while the source side and the edge energy agree to about 1e-11.
- **Coupled small-signal output** (DEVSIM full Jacobian, `G_E` equilibrated SuperLU, `G_M` UMFPACK). At 0.488174 V the errors are −358 and +361 relative, with a backward residual near 1e-16 and κ_G ≈ 3.8e16. At 0.7835 V they are −4.6 % and +6.7 %. A residual check alone accepts all of these.

## Reproduce

```bash
python benchmark/tools/export_cases.py      # deterministic; rewrites cases/, INDEX.json, MANIFEST.sha256
python benchmark/tools/make_examples.py
python -m pytest -q benchmark/tests
```

Every source file is checked to be byte-identical to commit `c05139b4ae037ece34257ff0382e84a91b059773`. Each case lists its sources with their SHA-256. The exporter only copies numbers, or applies the closed-form ΔL relation used by the producing script.

## Limits

- Only internal-path agreement is available so far. Nothing here is experimental validation.
- B-W and most C cases have no established uncertainty.
- The finite-bias groups use `s_couple = 0` (drain ψ held at equilibrium) except the three `s1` groups. `s_couple = 1` is the standard ohmic contact.
- The licence for public release is not yet decided.
