# Active research handoff — 2026-11-01

- Branch: `phase3/direction-reset-2026-11-01`, based on `fde679b561bddf35663933890d7ace617bf9b928`.
- Active question: test whether the weak-inversion end correction has an electron-free dimensionless limit `F(eta,chi)=Delta L/W_dep`, rather than continuing the already closed extraction-impact narrative.
- Closed extraction direction: C6–C9 established gate-dependent `Delta L`, affine conductance, extraction consequences, halo sign reversal and their prior-art overlap. The latest fit remains verification only.
- Well-depth correction: the 0.50→1.00 um interval crosses the proposed threshold and changes the bottom-boundary margin from 0.70 to 0.20 um; it is not a strict refutation of a `3 W_dep` semi-infinite criterion. See `phase3/reports/Direction_Reset_Correction_and_Next_Investment_2026-11-01.md`.
- New decision evidence: C8 E8-1 six-point data were reduced to `F=DeltaL/W_dep`, `eta`, `chi`, and `N_s/(N_A W_dep)`. Fitting `F=F_inf+A/chi` on alpha=15,50 predicts alpha=150 within 0.45% and 0.40% on the two existing structures. This is an exploratory mechanism test, not independent validation.
- Next gate: one third structure with substantially different `eta`, a preregistered prediction and an independent discretization. Proceed only if the prediction is within 5%; close the direction if it exceeds 15% or is boundary-dominated.
- Validation scope: the combined command `pytest -q semiconductor_sim phase3/tests` reports 171 passed, 4 skipped, 4 warnings. Production code and golden files remain unchanged.
