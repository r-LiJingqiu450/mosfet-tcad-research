# mosfet-tcad-research

MOSFET numerical simulation, model verification and reproducible research.

## Current baseline

The repository snapshot at `c64c463ba35fb3c2f177c9558217a20387046795` did not
contain the historical S1/S2/S3 implementations, tests, golden files, or raw
data. A small, explicit S1 baseline is now provided in
`common/poisson_boltzmann.py`: an intrinsic one-dimensional
Poisson–Boltzmann slab solved with damped Newton iteration. It is intended to
exercise numerical residuals, polarity symmetry, and grid refinement before
MOSFET/GCA physics is added.

Install the pinned development environment and run the tests with:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install --requirement requirements-dev.txt
.venv/bin/python -m pytest -q
```

The reproducible S1 experiment is run with:

```bash
MPLCONFIGDIR=/tmp/mplconfig .venv/bin/python experiments/run_s1_baseline.py
```

It writes machine-readable results and a report under
`results/s1_baseline/`. S2, S3, the GCA model, and any comparison to a
historical golden result remain pending because their source and reference
materials are absent from the repository.
