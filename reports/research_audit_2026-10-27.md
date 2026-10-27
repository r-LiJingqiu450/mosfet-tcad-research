# MOSFET/TCAD research audit and baseline delivery

## Repository and environment evidence

The inspected checkout is `r-LiJingqiu450/mosfet-tcad-research`, branch
`work`. Before development, `work` and `origin/main` both pointed to
`c64c463ba35fb3c2f177c9558217a20387046795`. The remote read succeeded through
the platform HTTPS Git proxy. The remote contains `main` and
`r-LiJingqiu450-patch-1`; no other branch, tag, or recoverable object contains
the historical research implementation.

That snapshot contained only `README.md`, `.gitignore`, and
`requirements-dev.txt`. There was no `S1`, `S2`, `S3`, `common` numerical
kernel, test suite, golden JSON, raw data, configuration, or prior report.
Consequently, the historical CLM divide-by-zero, threshold branching,
Newton stopping, and golden-file issues cannot be reproduced or attributed to
the current repository.

The development virtual environment uses Python 3.12 with pinned NumPy 2.3.5,
SciPy 1.17.0, Matplotlib 3.10.8, and pytest 9.1.1. The initial environment
only had pytest in `.venv`; adding the numerical packages to
`requirements-dev.txt` was required for the numerical tests.

## Implemented and verified

This branch adds a deliberately scoped S1 baseline: an intrinsic one
dimensional Poisson–Boltzmann slab with a prescribed surface potential and a
neutral far boundary. The solver uses a damped Newton iteration, reports a
dimensionless residual, refuses non-finite input, and estimates surface charge
with a second-order one-sided derivative.

The test suite contains three checks:

1. zero surface potential gives the exact neutral solution;
2. positive and negative potentials have odd potential and charge symmetry;
3. grid refinement reduces the change in surface charge.

The command `.venv/bin/python -m pytest -q` passes all 3 tests. The experiment
command `MPLCONFIGDIR=/tmp/mplconfig .venv/bin/python
experiments/run_s1_baseline.py` evaluates 9 surface potentials on 3 grids,
for 27 rows total. All rows converged, with maximum dimensionless residual
`1.67e-12`. The outputs are stored under `results/s1_baseline/` as JSON, CSV,
NPZ, PNG, a log, a report, and a manifest.

## Scientific scope and limitations

The experiment provides a numerical baseline and a grid study. It is not an
independent physical validation because no external reference solution or
historical golden data is present. The model omits oxide capacitance, doping,
fixed charge, mobility, channel current, CLM, threshold extraction, and
MOSFET/GCA behavior. S2, S3, and S4 remain unimplemented pending source code,
parameter definitions, and reference data.

No original data or trusted golden file was overwritten. The generated baseline
is identified as new work and must not be used as a replacement golden result
for the missing historical project.

## Delivery state

The implementation is committed on the `work` branch and `main` remains
unchanged. The branch was pushed successfully to GitHub. Pull request creation
through `gh pr create` was attempted but GitHub GraphQL returned
`Forbidden`; no PR number is claimed. The saved cloud environment draft now
uses the expanded dependency installation and startup validation, and records
the exact branch HEAD for reproducible restoration. The user must review and
publish that draft in environment settings before relying on a fresh task to
restore the prepared filesystem.
