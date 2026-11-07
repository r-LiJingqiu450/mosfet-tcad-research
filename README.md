# mosfet-tcad-research

MOSFET numerical simulation, model verification and reproducible research.

[![DOI](https://zenodo.org/badge/1435421516.svg)](https://zenodo.org/badge/latestdoi/1435421516)

## Published: low-current MOSFET drift–diffusion verification benchmark, v1.0.0

Version 1.0.0 (2026-11-07) is the first public release of this project. It is an open, code-agnostic reference set
for the zero-bias (linear-response) terminal conductances `G_DD` and `G_SD` of two-dimensional n-channel MOSFETs. The
devices have resolved n⁺ wells and operate in deep subthreshold (ψ̂ ≈ 12–21; drain current at 1 mV from 1.7e-15 to
1.6e-8 A/cm).

The release contains:
- a scorer and published meshes;
- a failure-mode atlas of terminal-current evaluation in DEVSIM 2.11.2;
- pre-registered numerical studies and a pre-registered AI-agent study;
- the accompanying manuscript, a preprint that has not been peer reviewed.

| | |
|---|---|
| **GitHub release** | [v1.0.0](https://github.com/r-LiJingqiu450/mosfet-tcad-research/releases/tag/v1.0.0) |
| **Zenodo, this version** | [10.5281/zenodo.24040430](https://doi.org/10.5281/zenodo.24040430) |
| **Zenodo, all versions** | [10.5281/zenodo.24040429](https://doi.org/10.5281/zenodo.24040429) (resolves to the latest version) |
| **Release snapshot** | tag `v1.0.0` = commit `7d450799ced424bcb716cc2884587e378225c397` (branch [`release/v1.0.0-candidate`](https://github.com/r-LiJingqiu450/mosfet-tcad-research/tree/release/v1.0.0-candidate)) |
| **Built from** | research commit `88f26a64995662e77c83a270ed84b5a567da70fd` (branch [`benchmark/wp8-mechanism-gci-2026-11-06`](https://github.com/r-LiJingqiu450/mosfet-tcad-research/tree/benchmark/wp8-mechanism-gci-2026-11-06)) |
| **Manuscript** | [`paper/main.pdf`](https://github.com/r-LiJingqiu450/mosfet-tcad-research/blob/v1.0.0/paper/main.pdf) in the release (preprint) |
| **Licences** | code: MIT ([`LICENSE`](https://github.com/r-LiJingqiu450/mosfet-tcad-research/blob/v1.0.0/LICENSE)); data, documentation and text: CC BY 4.0 ([`LICENSE-DATA.md`](https://github.com/r-LiJingqiu450/mosfet-tcad-research/blob/v1.0.0/LICENSE-DATA.md)) |

### Key findings

- **Reference values.**
  - The 6 primary three-terminal states (A-R) passed the acceptance criteria fixed before the runs. All 5 held-out states passed as well.
  - The tolerance is 2.4e-4 to 4.7e-3 relative on any mesh, and 1e-6 on the published meshes.
  - Under mesh refinement the references converge at second order (observed order 1.98 to 2.22).
- **Default contact current in DEVSIM 2.11.2.**
  - It has an absolute round-off floor of 2.4e-12 to 2.7e-11 A/cm and fails all 17 scorable cases.
  - Extended precision repairs it, at 10.9 to 15.7 times the wall time.
- **DEVSIM's AC analysis** (the coupled small-signal solve).
  - Its backward error is near machine precision, yet its output is wrong by up to 3.3e3 relative.
  - Solving the same double-precision system exactly does not repair it.
  - Rounding the Jacobian entries by one unit in the last place moves the output by 20 to 93% of the error.
  - An exploratory analysis traces the error floor (3.2e-9 to 8.6e-9 S/cm) to the neutral drain. There, the potential derivatives of the fluxes cancel in exact arithmetic.
- **AI-agent study** (pre-registered): 33 runs, 0 silent errors (95% upper bound 0.106).
- **Limits.**
  - Everything is agreement between computational paths, not experimental validation.
  - The physics envelope is narrow.
  - DEVSIM is the only external code scored.

### Cite

> Hu, X. (2026). *Low-current MOSFET drift-diffusion verification benchmark* (Version 1.0.0) [Data set]. Zenodo.
> https://doi.org/10.5281/zenodo.24040430

```bibtex
@misc{hu2026lowcurrent,
  author    = {Hu, Xiaocong},
  title     = {Low-current {MOSFET} drift-diffusion verification benchmark},
  version   = {1.0.0},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.24040430},
  url       = {https://doi.org/10.5281/zenodo.24040430}
}
```

**Author:** Xiaocong Hu, School of Physics, Guizhou University (ORCID
[0009-0001-5813-7638](https://orcid.org/0009-0001-5813-7638)).

**Machine-readable metadata:**
[`CITATION.cff`](https://github.com/r-LiJingqiu450/mosfet-tcad-research/blob/v1.0.0/CITATION.cff) in the release. To cite
every version at once, use the concept DOI 10.5281/zenodo.24040429.

### Reproduce

1. **Get the release.** Either:
   - from Zenodo: download the archived zip from https://doi.org/10.5281/zenodo.24040430 and unpack it. It unpacks to
     `r-LiJingqiu450-mosfet-tcad-research-7d45079/`. It is GitHub's source archive of the tag.
   - with git: `git clone --branch v1.0.0 --depth 1 https://github.com/r-LiJingqiu450/mosfet-tcad-research.git benchmark-v1.0.0`.

   Both hold the same 1,121 files.
2. **Check the files.** In the release directory, run `sha256sum -c SHA256SUMS` (on macOS: `shasum -a 256 -c SHA256SUMS`).
3. **Install the software.** Use Python 3.12 with the package versions of the runs:
   `python -m pip install numpy==2.3.5 scipy==1.17.0 matplotlib==3.10.8 pytest jsonschema`.
   - **Simulations:** rerunning them also needs `devsim==2.11.2`, with `DEVSIM_MATH_LIBS` pointing at BLAS/LAPACK
     (for example `/usr/lib/x86_64-linux-gnu/libblas.so.3:/usr/lib/x86_64-linux-gnu/liblapack.so.3`). Without DEVSIM,
     the tests that need it are skipped.
   - **Manuscript:** it needs LaTeX with latexmk, siunitx, booktabs, natbib and cleveref.
4. **Rerun** the checks, the registered evaluations, the manuscript's figures, tables and numbers, and the manuscript:

   ```bash
   python -m pytest -q benchmark/tests                      # scorer, examples, harness and release checks
   python benchmark/wp2/evaluate_wp2.py                     # registered evaluation of the reference runs
   python benchmark/wp4/evaluate_wp4.py                     # registered evaluation of DEVSIM's output paths
   python benchmark/wp6/harness/evaluate_wp6.py --stage 1   # registered evaluation of the agent study
   python benchmark/wp8/evaluate_wp8.py                     # registered evaluation of WP8 (mechanism, grid convergence)
   python benchmark/wp8/xsens_summary.py                    # summary of the exploratory WP8 sensitivity analysis
   python paper/scripts/make_assets.py                      # every figure, table and number of the manuscript
   cd paper && latexmk -pdf main.tex                        # the manuscript
   ```

   The release's acceptance run reran all of these. The regenerated outputs were byte-identical to the released files,
   and the rebuilt manuscript's text was identical to the released PDF.
5. **Score your own code.**
   1. Compute `G_DD` and `G_SD` for the A-R states, on any mesh or on a published mesh (declared `same_mesh`).
   2. Write a submission as described in
      [`benchmark/README.md`](https://github.com/r-LiJingqiu450/mosfet-tcad-research/blob/v1.0.0/benchmark/README.md).
   3. Run `python benchmark/tools/score.py my_submission.json --out my_report.json`.

   The scorer needs only the Python standard library. If `jsonschema` is installed, it also validates the submission.
6. **Rebuild the reproducible tarball** (optional).
   1. In a full clone, run `git checkout 88f26a64995662e77c83a270ed84b5a567da70fd`.
   2. Run `python benchmark/tools/build_release.py --version 1.0.0`.

   This writes `dist/mosfet-lowcurrent-dd-benchmark-1.0.0.tar.gz`. The release commit message records its SHA-256,
   `1ebaef951eba59661a79855873f53e3d5a169e8098c22921ee404c383a126f4c`.

More detail is in the release's
[`README.md`](https://github.com/r-LiJingqiu450/mosfet-tcad-research/blob/v1.0.0/README.md) and
[`RELEASE_NOTES.md`](https://github.com/r-LiJingqiu450/mosfet-tcad-research/blob/v1.0.0/RELEASE_NOTES.md).

### How this release was made

This release was produced with AI tools under the author's direction.

- **Codex (OpenAI):** used in the earlier project phases (auditing, numerical-method development, tests, reproducibility
  workflows, MOS-capacitor experiments, reports), whose code and results the case files cite.
- **Claude Opus 5.5 (Anthropic):**
  - designed the benchmark studies and their pre-registrations;
  - wrote the benchmark and analysis code;
  - ran the simulations;
  - produced the figures;
  - drafted the manuscript.

The pre-registrations, ledgers, run records and transcripts are included so that every step can be checked. The
manuscript's section *Use of generative AI* gives the full disclosure.

## Repository layout

| Branch | Content |
|---|---|
| `main` | This page and the development requirements. |
| `release/v1.0.0-candidate` | The v1.0.0 release snapshot (tag `v1.0.0`). It holds only the release files. |
| `benchmark/*` | The work packages of the benchmark (WP1–WP8, from 2026-11-04). The latest research branch is `benchmark/wp8-mechanism-gci-2026-11-06`. |
| `phase2a/*`, `phase2b/*`, `phase2c/*`, `phase3/*`, `codex/*`, `restore/*`, `work` | Earlier project phases: source recovery, verification audits, reference studies and research-direction reports. |
