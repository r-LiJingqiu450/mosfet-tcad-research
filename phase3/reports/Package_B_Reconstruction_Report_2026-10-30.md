# Package B — Reconstruction Report (2026-10-30)

**Specification:** Research_05 v2.1, Package B (B-1…B-5).
**Branch / commit:** `phase3/B-rebuild-2026-10-30` @ `43264ff` (parent: Phase 2C `7073f20`).
**Verdict:** **PASS in both environments.** Pinned: Python 3.12.3, NumPy 2.3.5, SciPy 1.17.0. Current: Python 3.13.18, NumPy 2.5.3, SciPy 1.18.1. The outputs of the two environments are bitwise identical.

## 1. Why a reconstruction

The original 2026-10-29 Package B files were never pushed:
- source: `phase3/m0.py`, `package_b.py`, `tests/test_m0.py`;
- data: `results/phase3/B`.

They were lost with the session container. A search found neither `mosfet-tcad-phase3-stage0-B-Dpre.gitbundle` nor `…_source.zip`. The places searched were:
- Project documents (`project_search`);
- conversation attachments and `/mnt/user-data/uploads`;
- the whole container filesystem (`find /` for `*gitbundle*`, `*.bundle`, `*stage0-B*`, `*Dpre*`).

The previous report (`claude/Phase3_Package_B_Report.md`) was used **only as an acceptance reference**. Every number below was recomputed.

## 2. What was built (production files read-only)

| File | Content |
|---|---|
| `phase3/m0.py` | Parameter-injected copy of `stage3_mosfet/reference.py`: `_em1mx` (l.16–19), `reference_slice` (l.22–44), `reference_current` (l.47–56). Parameters come from a frozen `M0Params` and tolerances from `Tolerances`, with a `.tightened(f)` variant. Helpers: `conductance0`, `ith_per_um`, `vt_cc_direct`, `vt_cc_tabulated` (PCHIP on ln I′), `delta_vg`, and `psi_b0` (−0.3571586 V). |
| `phase3/tests/test_m0.py` | 12 tests: bitwise slice equality with production; B-2 ≤ 1e-13; injection sensitivity; G0 as the V_D → 0 limit; agreement between the two threshold finders; finder status codes |
| `phase3/package_b.py` | B-1…B-5 driver. It records input SHA-256s, versions, timestamps and the golden SHA before and after the run |
| `phase3/runlog.py` | Run-record helper (R05 rule 3) |

## 3. Commands and exit codes

| Command | Exit | Wall time | Output |
|---|---|---|---|
| `python3 phase3/package_b.py --out results/phase3/B_rebuild_2026-10-30/current` | 0 | 9.06 s | `current/`, `current_stdout.log` |
| `venv312/bin/python phase3/package_b.py --out …/pinned` | 0 | ~9 s | `pinned/`, `pinned_stdout.log` |
| `pytest phase3/tests` (both environments) | 0 | 12 passed | `phase3_tests_{current,pinned}.log` |

## 4. Results vs the acceptance reference

| Check | Criterion | 2026-10-30 (both envs) | 2026-10-29 report | Status |
|---|---|---|---|---|
| B-1, 42 points of `independent_reference.csv`: max e_rep (production and m0) | ≤ 1e-12 pinned / ≤ 1e-9 current | **0.0 / 0.0** | 0.0 (current only) | PASS (pinned now run) |
| B-2: m0 vs production | ≤ 1e-13 | **0.0** | 0.0 | PASS |
| B-3: U_GCA,rel over all 42 points, tolerances ÷ 100 (estimate) | ≤ 1e-6 | **1.70e-8** (`quad` estimate dominates); tight − default ≤ 1.35e-14 | 1.70e-8 | PASS |
| B-5: 18 Phase 1 background shares | ≤ 1e-6 | **0.0** | 0.0 | PASS |
| B-5: R05 quoted values (6 significant digits) | within rounding | 3.2e-7, 8.5e-7, 2.5e-7 | — | PASS |
| R-II comparison window (f_bg < 1e-3, V_D = 50 mV) | — | V_G ≥ **0.50 V** | 0.50 V | reproduced |

**B-4 derived quantities** (`m0_F0.csv`: 40 V_G × 8 V_D):

| Quantity | Value |
|---|---|
| V_T,cc (direct root on M0, I′_th = 1e-7 A/µm, V_D = 50 mV) | **1.0428740 V** (reproduced) |
| PCHIP finder, 0.05 V grid | −0.114 mV (reproduced) |
| PCHIP finder, 0.01 V local grid (new) | **+0.0012 mV**, which confirms the ≤ 0.01 V recommendation for 2-D curves |
| M0 − Phase 1 production table | −0.400 mV (reproduced; a same-model cross-check, not an acceptance item) |
| max abs(I′/(V_D·G0) − 1) for V_D ≤ 1 mV | 1.909 % (vs 1 mV/(2V_T) = 1.934 %; reproduced) |

The CSV also provides ψ_s in both references: φ_s bulk-referenced, and ψ_s = φ_s + ψ_b0 intrinsic-referenced. These are the inputs needed for C-2 T0.

## 5. Not reconstructed

- `results/phase3/D_pre/`, the pre-DEVSIM GCI calibration with the in-house 1-D FV slice, and `phase3/meshmap.py` (D-1). These were not part of this authorization. **D-1 is needed before C-2.**
- Package A remains open: Phase 2D `01abe99` and G-0 `8e61dc76` are not on the remote.

Golden `stage3.json` SHA-256 `17ca3274…7828` is unchanged before and after every run.
