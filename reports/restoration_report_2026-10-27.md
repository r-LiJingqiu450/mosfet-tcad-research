# MOSFET Phase 1 v1.2 restoration report

## Recovery source and repository state

The restoration branch is `restore/phase1-v1.2`, created from the freshly
fetched `origin/main` at commit `c64c463ba35fb3c2f177c9558217a20387046795`.
The previous `work` branch was not used as the restoration base. The revision
archive was used as the source of restored files; `胡小从仿真项目.zip` was
inspected only for cross-checking and was never copied over the revision.

Archive SHA256 values are recorded in `verification/archive_sha256.txt`:

- `MOSFET_Phase1_Revision_v1_2 (1).zip`:
  `881ebe97cede3a19b08112453b404ef87052531bed1ec1456730126392f20618`
- `胡小从仿真项目.zip`:
  `cb03e8162363b983db32c9fe4c343ea90f22fabcf8de7dbd1e47371efeaa29ca`

The revision archive contains 112 files. Its release manifest has 111 entries;
all 111 file hashes match. The source manifest has 46 entries and matches the
embedded `baseline/recovered_v1_source.zip` exactly. Eleven of those entries
intentionally differ in the revised working tree because they are the files
changed by `changes.patch`; this is evidence of the revision, not a silent
overwrite. The embedded historical delivery archive is retained under
`history/`.

The only pre-existing repository file that conflicted with the archive was
`README.md`. Its original project description was merged into the revision
README. `.gitignore` was preserved. `requirements-dev.txt` now includes the
revision's `requirements.txt` and the verified environment pins:
NumPy 2.3.5, SciPy 1.17.0, Matplotlib 3.10.8, and pytest 9.1.1.

## Current reproducibility evidence

The following commands ran in the restored checkout:

```text
.venv/bin/python -m pip install --requirement requirements-dev.txt  exit 0
.venv/bin/python -m pytest semiconductor_sim -q                    exit 0
.venv/bin/python scripts/reproduce_phase1.py --fresh                exit 0
```

The direct suite and the reproduction suite both collected 95 tests and
reported 91 passed, 4 skipped, 0 failures, and 0 errors. The four skips are
the stage1 and stage2 golden regression/key-stability tests. The archive has
no `stage1.json` or `stage2.json`; only the S3 production golden exists at
`semiconductor_sim/results/golden/stage3.json`. The skipped tests therefore
remain unverified historical baselines, not passing tests.

The reproduction workflow completed all four stages: numerical study and
attribution, full pytest with JUnit output, report generation, and artifact
verification. It rewrote ordinary result and report outputs as part of the
documented workflow but did not regenerate or modify the frozen golden. The
S3 golden SHA256 was
`17ca32745b69c34355113347eb0f51f95772cf8a76f76b1decdd85925ef67828` both
before and after the run. The complete logs are `study_log.txt`,
`tests_phase1.txt`, `verification/pytest.xml`,
`verification/report_generation.txt`, `verification/artifact_check.txt`, and
`verification/run_status.json`.

The current run emitted four S2 warnings during the equilibrium degeneration
closure test: an Armijo line-search floor, a Newton iteration budget warning,
and two Gummel iteration budget warnings. The test still passed because the
closure assertion is configured to tolerate this difficult diagnostic case.
These warnings are retained as a numerical risk for the next S2 convergence
study.

## Defect audit

The historical v1 CLM path was reproduced with a minimal fake table. Calling
the recovered v1 solver with `clm=False` in a saturated case left
`dL=0` and `y_p=L`, then evaluated the pinch profile with denominator
`L-y_p=0`. The actual result was one `invalid value encountered in divide`
warning and non-finite profile values. The revision's guard
`if saturated and dL > 0.0` was exercised by the same case: the profile was
finite, no warning was emitted, and `y_p/L=1.0`. The focused regression
result is in `verification/clm_zero_length_audit.txt`, and the production
test `test_clm_disabled_profile_finite_at_high_drain` passes.

The threshold inverse has explicit linear and saturation branches. Both branch
tests pass. The production study reports the saturation branch at the chosen
criterion, with `vt_cc=1.0432740076 V`, `vt_long=0.9963499971 V`, and a
criterion offset of `29.2925 mV`; the branch is checked by reconstructing the
requested current rather than by a label alone.

The revised Poisson solver requires both a correction criterion and a scaled
residual/backward-error criterion. A one-iteration budget is explicitly tested
to remain non-converged. The full stage3 solver tests pass, including table
rejection of non-converged slices.

The independent continuum reference tests pass at four bias points. The
reported maximum current relative error is `0.0019098835` (0.19099%). The
spatial refinement study reports source sheet-density errors of 3.3123%,
0.8962%, 0.2357%, and 0.06065% as the initial spacing changes from 0.5 nm to
0.0625 nm. This establishes numerical refinement evidence for the current
reduced model; it is not experimental validation or a full 2D drift-diffusion
reference.

## Remaining material and next research step

Stage1 and stage2 historical golden files are still absent. The original
project ZIP contains document/PDF artifacts rather than a second executable
source tree, so it cannot supply those missing machine-readable baselines.
The restored package also explicitly limits itself to a classical GCA model;
it does not claim 2D drift-diffusion, quantum transport, mobility/AC effects,
or a finished paper.

The next defensible step is to add independently generated S1/S2 golden files
only after defining and reviewing their configurations, then rerun the same
freeze-protected workflow. In parallel, the S2 equilibrium Newton/Gummel
warnings should be isolated with iteration histories and residual thresholds.
Only after those baselines are independently checked should S4 or any new
physical model be considered.
