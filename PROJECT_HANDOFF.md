# Active research handoff — 2026-11-01

- Branch: `phase3/direction-reset-2026-11-01`, based on `fde679b561bddf35663933890d7ace617bf9b928`.
- Active question: quantify how the already observed bias-dependent intrinsic `Delta L(V_G)` changes practical multi-length and gate-dependent extraction results.
- Well-depth correction: the 0.50→1.00 um interval changes `Delta L` by 1.47115%, but it crosses the proposed threshold and also changes the bottom-boundary margin from 0.70 to 0.20 um. It is retained as a boundary-sensitivity check, not as a strict refutation of a `3 W_dep` semi-infinite criterion. See `phase3/reports/Direction_Reset_Correction_and_Next_Investment_2026-11-01.md`.
- New evidence: existing C6 `L=2,4,8 um` states fit the affine `1/G` law with `Delta L_fit` = 369.865, 153.686, 53.385 nm at `V_G` = 1.0, 1.25, 2.0 V; maximum relative fit residuals are `3.71e-6`, `9.41e-8`, `5.77e-9`.
- DEVSIM comparison correction: the machine-readable network-vs-E relative difference is `3.1334709456e-7`; the second is `3.3407079818e-7`; the earlier `1.57e-7` was a half-difference reporting error.
- Next specification: `phase3/next_investment/Extraction_Impact_Specification_2026-11-01.md`; current post-processing outputs are under `results/phase3/next_investment/extraction_impact/`.
- Validation scope: the combined command `pytest -q semiconductor_sim phase3/tests` reports 171 passed, 4 skipped, 4 warnings. Production code and golden files remain unchanged.
