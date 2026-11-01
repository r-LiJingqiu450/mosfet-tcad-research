# Codex research handover at C12

**Verified base:** `phase3/C12-multigate-literature-2026-11-01` at `0d2390c8d49741eddd4d605a60cf25b3b49c0139`.
This document records the state reconstructed from the restored Git bundle and the independently fetched GitHub repository. It adds no simulation data and changes no production model.

## Evidence state

The C12 branch contains the complete Phase 3 source, specifications, reports, result JSON/CSV/NPZ files, figures and logs. Bundle restoration passed `git bundle verify` and `git fsck --full --strict`; the subsequent GitHub fetch passed an independent `git fsck`. Every 1,297 C12 tracked files matches the supplied research ZIP byte-for-byte. Four ZIP export records remain outside Git: `ENVIRONMENT.txt`, `GIT_INFO.txt`, `MANIFEST_SHA256.txt`, and `REPRODUCE.md`.

The protected S3 golden is unchanged at SHA-256 `17ca32745b69c34355113347eb0f51f95772cf8a76f76b1decdd85925ef67828`. Phase 2B candidate S1/S2 hashes are `93244365320f005eb61ae242f5416cbb4cf9a9af92b1e63917370e336adc0bd0` and `7883dfdd80fe5fc3b768ed3e21088f6b212682dd74048c033a935dd1162afeac`.

## Scientific state reconstructed from C1–C12

- C1 established partial DEVSIM capabilities: AC, matrix export and current extraction work, with a triple-point Si/SiO₂ interface defect recorded.
- C3–C4 closed the numerical reference scope as `REFERENCE_READY` with caveats; novelty remained `NOVELTY_UNPROVEN`. T1 v2 passed its registered coverage, while the conservative L=10 µm level-3 conditioning case remains non-resolving.
- C5 rejected the original idealized GCA error-decomposition direction after the Ward full-text audit and recommended a gated source/drain end-correction direction.
- C6–C9 found a realistic near-threshold end correction and an affine-in-length behavior over the tested range. The scaling law was supported for the registered Addendum A and partial for Addendum B; strong-inversion and halo limits were qualified.
- C10 found several pre-registered similarity and asymptotic hypotheses rejected, and identified deep-subthreshold drain-side network-current cancellation. Independent finite-volume checks retained the numerical values.
- C11 rejected both parameter-free predictors tested: H-lin errors −13.3% to −24.2%, H-q2D errors −12.0% to −35.1%. The nonlinear pre-asymptotic transition is the current physical explanation; it remains a finding, not a universal theorem.
- C12 established that the proposed double-gate end correction is a special case of Liang–Taur (2004) in the stated limits. The multi-gate end-correction campaign is therefore closed as the primary publication line. O1 remains conditional on a literature gate; O2/O3 remain hypotheses or supporting mathematics.

## Recommended next decision

Do not launch a new multi-gate simulation campaign before reading the unresolved Xie/Taur, Taur–Ning, Yu/Taur, BSIM-CMG and conformal-mapping sources listed in `phase3/reports/C12_MultiGate_Literature_Assessment_2026-11-01.md`. If the literature gate finds the logarithmic/algebraic correction already established, close O1. If it does not, run only the pre-registered small decisive test with explicit falsification thresholds.

Current status is `MORE_EVIDENCE_NEEDED`; no SCI novelty claim is supported. The existing numerical limitations remain: deep-subthreshold coupled small-signal conditioning, unresolved predictive theory for the nonlinear boundary amplitude, and the narrow scope of the tested realistic-device scaling.

## Reproduction

The authoritative C12 branch and artifacts are in the GitHub repository. The restoration and backup audit, archive hashes, remote branch verification and object-reachability checks are in the separate G-0 evidence package. A future agent should begin by checking this branch against remote C12, then create a new research branch for any registered experiment. The `main` branch remains protected.
