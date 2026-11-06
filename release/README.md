# Low-current MOSFET drift–diffusion verification benchmark

An open, code-agnostic reference set for the **zero-bias (linear-response) terminal conductances of two-dimensional
n-channel MOSFETs at low current**. It comes with a scorer, published meshes, documented failure signatures of common
ways of evaluating the terminal current (a failure-mode atlas), and a pre-registered AI-agent study.

The accompanying manuscript is in `paper/` (`paper/main.pdf`).

**Version 1.0.0.** Release notes: `RELEASE_NOTES.md`.

## What is in it

| Path | Content |
|---|---|
| `benchmark/cases/` | 106 case files (schema 1.0.0 / 1.1.0), each with the complete device specification, reference values, uncertainty and provenance (SHA-256 of every source file) |
| `benchmark/meshes/` | the exact node lines and triangulation rules of every reference mesh |
| `benchmark/schema/` | JSON schemas of a case and of a submission |
| `benchmark/tools/score.py` | the scorer (Python standard library only) |
| `benchmark/examples/` | worked submissions and their score reports |
| `benchmark/README.md` | the benchmark definition: scope, observables, tiers, scoring rules, atlas |
| `benchmark/specs/` | the pre-registrations (with their dated amendments) |
| `benchmark/reports/` | the reports of every work package |
| `benchmark/wp2/`, `benchmark/wp4/` | the production runs of the reference values and of the external paths: code, ledgers, records, evaluations |
| `benchmark/wp6/` | the AI-agent study: harness, task kit, schedule, every run record and transcript, audits, evaluations |
| `results/phase3/…` | earlier project results that the internal-tier case files cite as sources (unchanged, at their original paths) |
| `phase3/…`, `semiconductor_sim/common/…` | the project modules that the benchmark code imports or that case files name as producing scripts (unchanged, at their original paths) |
| `paper/` | manuscript source, the script that makes every figure, table and quoted number from the files above, and the built PDF |
| `SHA256SUMS` | SHA-256 of every file in this release |

Every file sits at the path it had in the research repository, so every SHA-256 recorded in a case file and every
import in the benchmark code resolve inside this release. The research repository's full history is not part of the
release; the release snapshot records the commit it was built from (see `RELEASE_NOTES.md`).

## Quick start: score your code

1. Compute `G_DD` (and `G_SD`) in S/cm per unit width as the `V_D → 0` limit of your terminal currents for the A-R
   states (any mesh), or on a published mesh (`benchmark/meshes/`) and declare `same_mesh`.
2. Write a submission (format: `benchmark/README.md`, *Submitting*; schema: `benchmark/schema/submission.schema.json`).
3. Score it:

```bash
python benchmark/tools/score.py my_submission.json --out my_report.json
```

The scope envelope (Boltzmann statistics, constant mobility, holes at equilibrium, no recombination, T = 300 K, V_FB = 0)
and the scoring rules are in `benchmark/README.md`.

## Reproduce

```bash
python -m pytest -q benchmark/tests                      # scorer, examples, harness and release checks
python benchmark/wp2/evaluate_wp2.py                     # registered evaluation of the reference runs
python benchmark/wp4/evaluate_wp4.py                     # registered evaluation of DEVSIM's output paths
python benchmark/wp6/harness/evaluate_wp6.py --stage 1   # registered evaluation of the agent study
python paper/scripts/make_assets.py                      # every figure, table and number of the manuscript
cd paper && latexmk -pdf main.tex                        # the manuscript
```

Rerunning the simulations themselves needs DEVSIM 2.11.2 with `DEVSIM_MATH_LIBS` pointing at BLAS/LAPACK
(e.g. `/usr/lib/x86_64-linux-gnu/libblas.so.3:/usr/lib/x86_64-linux-gnu/liblapack.so.3`); extended precision needs a DEVSIM
build that reports `extended_precision: True`. The runs used Python 3.12, NumPy 2.3.5, SciPy 1.17.0 and, for the figures,
Matplotlib 3.10.8. The job lists and commands are in `benchmark/README.md` and the work-package reports.

## Licences

- **Code** (every `*.py` file, including the scorer, the harnesses and the project modules): MIT License, `LICENSE`.
- **Data, documentation and text** (case files, meshes, schemas, examples, run records, transcripts, reports,
  pre-registrations, the manuscript and its figures): Creative Commons Attribution 4.0 International, `LICENSE-DATA.md`.
- **Third-party files** keep their own licences: `benchmark/wp6/kit/docs/` holds DEVSIM's manual and examples
  (Apache License 2.0; see `DEVSIM_LICENSE` and `DEVSIM_NOTICE` there).

## Citing

See `CITATION.cff`. A DOI will be added when the release is archived.

## Limits

Everything here is agreement between computational paths; nothing is experimental validation. The physics envelope is
narrow. The only external code scored is DEVSIM. The agent study tested one model in one configuration. See
`benchmark/README.md` (*Limits*) and the manuscript (*Discussion*).
