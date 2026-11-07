# Release notes

## 1.0.0 (2026-11-07)

First public release. Release candidate rebuilt on 2026-11-07 with the WP8 results; the date is set to the publication
date when the release is published. The snapshot was built by `benchmark/tools/build_release.py` from the research
repository; the commit message of the release commit records the source commit and the SHA-256 of the tarball.

**Benchmark (schema 1.0.0 and 1.1.0, 106 case files).**
- **A-R** (6 three-terminal states, two devices, ψ̂ ≈ 12, 16, 21): reference `G_DD`, `G_SD` with a tolerance usable on any
  mesh (2.4e-4 to 4.7e-3 relative) and same-mesh references (1e-6) on three published meshes. All 6 passed the
  pre-registered acceptance (H1).
- **HO** (5 held-out states): same-mesh references on the level-1 mesh; all 5 passed H2.
- **B-R** (23 electron-blocking-body companions), **A-W** (3), **B-W** (28) and **C** (41) tiers from earlier project data.
- Scorer (`benchmark/tools/score.py`), schemas, published meshes and worked examples.

**Failure-mode atlas** (DEVSIM 2.11.2, these devices and meshes).
- The drain-side flux sum misses the same-mesh reference on every three-terminal mesh, by up to about 4e3 relative.
- DEVSIM's default contact current has an absolute round-off floor of about 2e-12 to 3e-11 A/cm and fails all 17
  scorable cases. Extended precision repairs it (6 of 6 A-R states) at about 11 to 16 times the wall time.
- DEVSIM's AC analysis returns the coupled small-signal solve, whose backward error is near machine precision while the
  output is wrong by up to about 3e3 relative. It fails on all 6 A-R states in double and in extended precision.
- Why (WP8, pre-registered): an exact solve of the same double-precision system is still wrong (4.2e-4 to 1.6e3
  relative), and one-ulp rounding of the Jacobian entries moves the output by 20 to 93% of the error. An exploratory
  analysis places the sensitivity in the neutral drain, where the potential derivatives of the fluxes cancel; their
  rounding leaves a floor of a few 1e-9 S/cm that does not depend on the gate bias.

**Grid convergence** (WP8, pre-registered). The A-R references converge at second order (observed order 1.98 to 2.22);
every registered tolerance contains the Richardson-extrapolated value. The WP8 code, ledger and records are in
`benchmark/wp8/`; report: `benchmark/reports/WP8_Report_SmallSignal_Mechanism_and_Grid_Convergence_2026-11-06.md`.

**AI-agent study** (pre-registered; Claude Opus 5.5 in Claude Code, arm A without the scorer).
- 33 runs: 0 silent errors (95% upper bound 0.106), 32 correct and claimed, 1 without a report (time limit).
- The registered stop rule applied, so arm B was not run.

**Manuscript** (`paper/`): LaTeX source, the asset script and the built PDF, updated with the WP8 results. The
author's approved AI-use disclosure, the corresponding e-mail, the acknowledgments, the repository URL and the DOI are
to be added.

**Known limits.** Agreement between computational paths only (no experimental validation); narrow physics envelope;
DEVSIM is the only external code scored (the NiTCAD evaluation is on hold); one agent configuration tested.
