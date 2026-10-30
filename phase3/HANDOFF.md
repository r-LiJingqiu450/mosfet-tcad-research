# Phase 3 Research Handoff — 2026-10-30

**Final status: `C1_PARTIAL`.** C-1 (a), (b) and (c) are all feasible. Two items remain open: the C-1(b) acceptance wording, and the triple-point deck defect that blocks C-3. **C-2 may begin** under the conditions in `reports/C1_Feasibility_Report_2026-10-30.md` §7. **Package E was not started.**

## Branches and commits (this handoff)

| Branch | HEAD | Parent | Content |
|---|---|---|---|
| `phase3/B-rebuild-2026-10-30` | `43264ff` | `7073f20` (Phase 2C, remote) | Package B reconstruction |
| `phase3/C1-2026-10-30` | see `git log` (handoff commit on top of `f4ca0bf`) | `43264ff` | C-1 (a)(b)(c) feasibility, reports, this file |

`main` (`c64c463`) and all remote branches are untouched. No force operation was used. The branches were **not pushed**: the session could not obtain GitHub write access. They are preserved in the git bundle `mosfet-tcad-phase3-B-C1-2026-10-30.bundle`, which includes the Phase 2C ancestry, and in a source ZIP.

## Restore

```bash
git bundle verify mosfet-tcad-phase3-B-C1-2026-10-30.bundle
git clone mosfet-tcad-phase3-B-C1-2026-10-30.bundle repo && cd repo
git branch -a        # phase3/B-rebuild-2026-10-30, phase3/C1-2026-10-30, phase2c/...
# to publish (from a machine with write access; never main):
git remote add gh https://github.com/r-LiJingqiu450/mosfet-tcad-research
git push gh phase3/B-rebuild-2026-10-30 phase3/C1-2026-10-30
```

## Reproduce

```bash
python3.12 -m venv venv312 && venv312/bin/pip install devsim==2.11.2 numpy==2.3.5 scipy==1.17.0 matplotlib==3.10.8 pytest==9.1.1
export DEVSIM_MATH_LIBS=/lib/x86_64-linux-gnu/libblas.so.3:/lib/x86_64-linux-gnu/liblapack.so.3
venv312/bin/python -m pytest semiconductor_sim phase3/tests -q -p no:cacheprovider    # 106 passed, 4 skipped
venv312/bin/python phase3/package_b.py --out <new dir>
venv312/bin/python phase3/c1/c1b_jacobian.py --out <new dir>
venv312/bin/python phase3/c1/c1a_ac.py --lateral onecol --oxide-corner-pin 1 --out <new dir>
venv312/bin/python phase3/c1/c1c_currents.py --out <new dir>
```

All three C-1 scripts reproduced bit-for-bit from commit `f4ca0bf` on a clean tree (`results/phase3/C1/repro_f4ca0bf/PROVENANCE.txt`).

## File index

| Path | What |
|---|---|
| `phase3/m0.py`, `phase3/package_b.py`, `phase3/tests/test_m0.py`, `phase3/runlog.py` | Package B (M0 reference) |
| `results/phase3/B_rebuild_2026-10-30/{current,pinned}/` | B-1…B-5 outputs, run JSON, figures, logs |
| `phase3/c1/deck.py`, `onecol.py`, `measure.py` | C-1 DEVSIM deck, independent 1COL FV solver, measurement helpers |
| `phase3/c1/c1a_ac.py`, `c1b_jacobian.py`, `c1c_currents.py` | C-1 (a), (b), (c) test scripts |
| `results/phase3/C1/c1{a,b,c}/attempt*/` | Every attempt, including failures, with `EXIT.txt` notes and full `run.log` |
| `results/phase3/C1/repro_f4ca0bf/` | Clean-tree reproduction runs (final evidence) |
| `results/phase3/C1/regression/` | pytest logs in both environments, golden SHA before/after, scratch-file isolation record |
| `phase3/reports/Package_B_Reconstruction_Report_2026-10-30.md` | Package B report |
| `phase3/reports/C1_Feasibility_Report_2026-10-30.md` | C-1 report: versions, solver parameters, results, capability matrix, C-2 recommendation, limitations |

## Protected items (verified)

- Golden `semiconductor_sim/results/golden/stage3.json` SHA-256 `17ca32745b69c34355113347eb0f51f95772cf8a76f76b1decdd85925ef67828` is unchanged.
- `git diff 7073f20 -- semiconductor_sim` is empty: no production solver, test or dataset was modified.
- The test-suite scratch file `results/golden/_test_common_npz.npz` is untracked and excluded from all commits and archives (documented in `results/phase3/C1/regression/scratch_file_isolation.txt`).
- No credentials, tokens, virtual environments or downloaded third-party sources are committed.
