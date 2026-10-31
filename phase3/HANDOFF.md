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


---

## Update 2026-10-30 (C-2 package) — commit reconciliation

The table above was written inside commit `5b26809`, which cannot contain its own SHA. Full SHAs, verified with `git log` and `git bundle list-heads`:

| Ref | Full SHA | Parent |
|---|---|---|
| Phase 2C baseline | `7073f208003ab32eb02bb3f75a75010e47338cc3` | — |
| `phase3/B-rebuild-2026-10-30` | `43264ff86f71c4d382540b250ba4fa668046459c` | 7073f20 |
| C-1 code | `f4ca0bfaed56b099448cad4d22b15636fceee30c` | 43264ff |
| `phase3/C1-2026-10-30` (final HEAD; archived as `…_source_5b26809.zip` and in the first bundle) | `5b26809a2c02e53cb4fb23187b4876c0558595e3` | f4ca0bf |
| C-2 code and results | `fe44d689597150da64c805183760fd5e696eeb99` | 5b26809 |
| `phase3/C2-2026-10-30` final HEAD | given in `SHA256SUMS`/`bundle_verification.txt` of the C-2 archive (a commit cannot record its own SHA) | fe44d68 |

C-2 reports: `phase3/reports/Corner_Interface_Defect_Report_2026-10-30.md`, `C1b_Numerical_Acceptance_Assessment_2026-10-30.md`, `C2_OneColumn_T0T2_Report_2026-10-30.md`.

Independent re-execution: `results/phase3/C2/independent_verification_2026-10-30T14/VERIFICATION.md`. Remote preservation: **not pushed**; the archive is a verified bundle plus ZIP.

---

## Update 2026-10-30 (C-3 package) — `REFERENCE_PARTIAL`, `NOVELTY_UNPROVEN`

Branch `phase3/C3-2026-10-30` from `c909bf11a0c6b4a41dcef613be7a4bb0d4a29b61` (C-2 final HEAD, verified). The final HEAD is in the archive's `bundle_verification.txt`. **Not pushed** (no write access; not retried). `main` and the remote branches are untouched. Package E is not started.

**Decision:** `phase3/reports/C3_Research_Decision_Report_2026-10-30.md`.

**Key commits:**

| Commit | Content |
|---|---|
| `9dc3307` | T1 driver with pre-registered E definitions (`phase3/c3/t1.py`); 1COL φ_n |
| `dc3f069` | `body_corners` option (R-I(b) only; default off); `fb_preview`; T1 analysis |
| `896d9d9` | T1 matrix complete (12 runs, 72 points); analysis v1 |
| `b8ad92d`, `352daf7`, `f2f1c5f`, `87752c9` | Post-hoc L1 linear-solve diagnostic (runs 1–3) |
| `9a88e17`, `318ca54` | Observed-order correction (non-constant ratios), unit tests, analysis v2, Research_08 Erratum E1 |
| `34b0b59` | Final regression: 146 passed, 4 skipped, both environments; golden unchanged |

**Reproduce (pinned venv, as above):**

```bash
export DEVSIM_MATH_LIBS=/lib/x86_64-linux-gnu/libblas.so.3:/lib/x86_64-linux-gnu/liblapack.so.3
venv312/bin/python phase3/c3/t1.py --L 1 --level 0 --modes a,b --vgs 0.8,1.25,2.0 --l3 1 --out <new dir>   # ~95 s
bash phase3/c3/run_t1_matrix.sh <out_dir> "1 0 1" "1 1 0" "1 2 0" "1 3 0" "3 0 1" "3 1 0" "3 2 0" "3 3 0" "10 0 1" "10 1 0" "10 2 0" "10 3 0"   # full matrix, "L level l3"; ~45 min on 2 cores split in two streams; the script hard-codes PY=/home/claude/venv312/bin/python
venv312/bin/python phase3/c3/t1_analysis.py --in results/phase3/C3/t1 --out <new dir>
venv312/bin/python phase3/c3/rib_inspect.py --L 1 --level 0 --vg 2.0 --body-corners 0 --out <new dir>
venv312/bin/python phase3/c3/l1_refine_diag.py --L 10 --level 3 --mode a --vg 2.0 --out <new dir>
venv312/bin/python phase3/c3/fb_preview.py --L 1 --level 0 --out <new dir>
```

**Reports:**
- `C3_T1_RIa_RIb_Report_2026-10-30.md`
- `C1b_Final_Acceptance_Review_2026-10-30.md` (+ Addendum A)
- `../specs/Research_08_Specification_Corrections_v2.2.md` (+ Erratum E1)
- `Finite_Bias_Discrepancy_Experiment_Design_2026-10-30.md`
- `Novelty_Challenge_Matrix_2026-10-30.md` (+ C-3 addendum)
- `C3_Research_Decision_Report_2026-10-30.md`

**User action needed for the literature gate:** download and attach the two Caltech PDFs, Ward, Odeh & Cohen 1990 (CaltechAUTHORS record 31400) and Ward 1988 PhD thesis (CaltechTHESIS 440). This environment's fetch proxy refused them.

---

## Update 2026-10-31 (C-4 package) — `REFERENCE_READY` (scoped), `NOVELTY_UNPROVEN`

Branch `phase3/C4-2026-10-30` from `e4e4528` (C-3 head, verified). The final HEAD is in the archive's `bundle_verification.txt`. **Not pushed** (no write access; not retried). `main`, golden files and production S1–S3 are untouched.

**Decision:** `phase3/reports/C4_Research_Decision_Report_2026-10-30.md`.
- The registered 48-state FB-1 stop/go is **not** recommended in its current form.
- Package E is not started.

**Key commits:**

| Commit | Content |
|---|---|
| `409e521` | **Pre-registration** of T1 v2, corrected R-I(b) and R-I(c) (`phase3/specs/C4_T1v2_Preregistration_2026-10-30.md`); `phase3/c4/{deck4,t1v2}.py` |
| `eb8309b` | Perturbation analysis (Props 1–3) and registered TC-1 rules |
| `3a18f37`, `d7846e1`, `bfcca33`, `1ede271` | T1 v2 sets B, C, D, D′ (64 runs, 192 points) and TC-1 (6 runs) |
| `daa183e` | `phase3/c4/jobrunner.py` (checkpointing, bounded runner; replaced the open-ended wait chains) |
| `b6cee5a`, `11068c7` | Analyses (T1 v2; ℓ_res and TC-1) |

**Reproduce (pinned venv):**

```bash
export DEVSIM_MATH_LIBS=/lib/x86_64-linux-gnu/libblas.so.3:/lib/x86_64-linux-gnu/liblapack.so.3
venv312/bin/python phase3/c4/t1v2.py --model b --L 1 --level 0 --out <dir>             # ~5 s
venv312/bin/python phase3/c4/t1v2.py --model c --lres 0.3 --L 1 --level 0 --out <dir>
venv312/bin/python phase3/c4/fb4.py --model c --lres 1 --L 1 --level 0 --out <dir>
venv312/bin/python phase3/c4/tc1.py --L 1 --level 0 --out <dir>                        # ~40 s
venv312/bin/python phase3/c4/jobrunner.py --jobs results/phase3/C4/jobs_fb4_c.json --out <dir> --limit-s 2400
venv312/bin/python phase3/c4/t1v2_analysis.py --in results/phase3/C4/t1v2_B results/phase3/C4/t1v2_C results/phase3/C4/t1v2_D results/phase3/C4/t1v2_Dprime --v1 results/phase3/C3/t1 --out <dir>
venv312/bin/python phase3/c4/c4_post_analysis.py --fb4 results/phase3/C4/fb4 --tc1 results/phase3/C4/tc1 --out <dir>
```

**User action needed:** attach the two Ward PDFs, Ward–Odeh–Cohen 1990 (CaltechAUTHORS 31400) and the Ward 1988 thesis (CaltechTHESIS 440). The Task E audit is blocked without them.

---

## C-5 (literature audit and direction decision), 2026-10-31

**Branch:** `phase3/C5-literature-2026-10-31`, from `d3cc880`. No simulation was run.

**Reports:**
- `phase3/reports/C5_Ward_FullText_Audit_2026-10-31.md` — equation-level audit of Ward 1988 and Ward–Odeh–Cohen 1990. Printed p. 142 of the thesis is missing; the affected items (§6.5 only) are marked.
- `phase3/reports/C5_Research_Direction_Decision_2026-10-31.md` — answers to the five questions, the E-W1 design and the abandonment criteria.

**Analysis code:** `phase3/c5/ward64_quadrature.py` → `results/phase3/C5/ward64/ward64_predictions.json`. This is the a-priori evaluation of Ward's published §6.4 formulas and the pre-registered predictions for E-W1.

**Decisions:**
- D1: rejected.
- C-W (R-II, verification and decomposition of Ward's O(ε) end correction): recommended, gated by literature gate LG-1 and then E-W1.
- E-W1 needs authorization and a frozen pre-registration.

**Reproduce:** `python3 -m phase3.c5.ward64_quadrature results/phase3/C5/ward64` (scipy only; ~10 s).

---

## C-6 (exploratory research phase: zero-bias end correction), 2026-10-31

**Branch:** `phase3/C6-endcorrection-2026-10-31`, from `2922f86`.

**Reports:**
- `phase3/reports/C6_Research_Phase_Report_2026-10-31.md` — findings, paths, candidate, uncertainties, recommendation.
- `phase3/reports/C6_Theory_ZeroBias_EndCorrection_2026-10-31.md` — P4–P7 and the end-layer observation.

**Code:**
- `phase3/c6/wprobe.py` — Ward-truncated n⁺ ends.
- `phase3/c6/rprobe.py` — resolved n⁺ wells.
- `phase3/c6/w1_analysis.py`, `phase3/c6/w1_figures.py`.
- Test: `phase3/tests/test_c6_network.py`.

**Data:** `results/phase3/C6/{w1, dev, analysis_w1}`, plus the job lists `jobs_w1_*.json`, the runner outputs and `CANCELLED.txt`.

**Reproduce one state** (about 1–3 min at level 0; level ≥ 2 needs `--no-devsim` within 7 GB):
`DEVSIM_MATH_LIBS=... venv312/bin/python phase3/c6/wprobe.py --L 4 --level 0 --vg 1.25 --out <dir>`

**Analysis:** `python3 phase3/c6/w1_analysis.py && venv312/bin/python phase3/c6/w1_figures.py`

---

## C-7 (realism, model test, scaling law), 2026-10-31

**Branch:** `phase3/C7-realism-model-2026-10-31`, from `08d61cb`.

**Report:** `phase3/reports/C7_Research_Phase_Report_2026-10-31.md`.

**Registrations:** `phase3/specs/C7_Realism_Test_Registration_2026-10-31.md` and `phase3/specs/C7_M7_Model_Preregistration_2026-10-31.md` (with Addenda A and B).

**Code:** `phase3/c7/`.
- `endprobe.py`: decks W/R, N_A, t_ox, overlap, grading, ψ_bi shift, Lombardi post-processing.

**Reproduce one state** (scaled device, about 15 s):
`DEVSIM_MATH_LIBS=... venv312/bin/python phase3/c7/endprobe.py --deck W --L 1 --level 0 --vg 1.2914 --tox 2 --na 1e18 --out <dir>`

**Analyses:**
`venv312/bin/python phase3/c7/c7_analysis.py; venv312/bin/python phase3/c7/c7_scaling.py; python3 phase3/c7/m7_failure_diag.py`

---

## C-8 (theory of the end correction; extraction consequences), 2026-10-31

**Branch:** `phase3/C8-theory-extraction-2026-10-31`, from `2cfa8e0`.

**Report:** `phase3/reports/C8_Research_Phase_Report_2026-10-31.md`.

**Registration:** `phase3/specs/C8_Preregistration_2026-10-31.md` (commit `8650c76`, before the runs).

**Reproduce:**
- E8-1: `python3 phase3/c8/e8_subth_eval.py`
- E8-2: `venv312/bin/python phase3/c8/e8_extraction.py`

The states come from `phase3/c7/endprobe.py`; job lists are in `results/phase3/C8/jobs_*.json`. Run them with `phase3/c4/jobrunner.py`.

---

## C-9 (independent audit of C-8; halo realism; originality reassessment), 2026-10-31

**Branch:** `phase3/C9-audit-realism-2026-10-31`, from `536e1ac`.

**Report:** `phase3/reports/C9_Research_Phase_Report_2026-10-31.md`.

**Registration:** `phase3/specs/C9_Preregistration_2026-10-31.md` (commit `b83892f`).

**Probe:** `phase3/c9/probe9.py`.
- Wraps `phase3/c7/endprobe.py` and adds a Gaussian halo, terminal-charge outputs, and short-L meshing.

**Reproduce one state** (about 15 s):
`venv312/bin/python phase3/c9/probe9.py --deck R --L 1 --vg 1.32 --out <dir>`

**Batches:** `phase3/c4/jobrunner.py --jobs results/phase3/C9/jobs_A.json` (and `jobs_B.json`).

**Analysis:** `venv312/bin/python phase3/c9/e9_analysis.py`
