# C-9 post-hoc diagnostics (not pre-registered)

These runs were made **after** the C-9 report (commit `1edc8ad`), during the pre-delivery audit of the package.
They are disclosed in the C-9 report: erratum E9-1 and addendum A9. None of them changes a pre-registered verdict.

| Directory | Code version | Content |
|---|---|---|
| `aend0.175/` | probe9 v1 (`probe9-2026-10-31`) via `phase3/c9/diag9_mesh.py` | D1: L = 1 µm solved on the short-device end mesh (lateral end segment a = 0.175 µm), 2e18 halo, V_G = 1.32 V |
| `level1/` | probe9 v1 | D2: L = 0.5 µm, level 1; D3 + sweep: L = 1 µm, level 1, 2e18 halo, 11 V_G; controls at 1.32 V (5e17 halo, no halo) |
| `level1_v2/` | probe9 v2 (`probe9-2026-10-31b`, Lombardi impurity fix) | the same level-1 states for the halo devices |
| `aend0.175_v2/` | probe9 v2 via `diag9_mesh.py` | D1 repeated with v2 |

- **Valid in all directories:** constant-mobility quantities (`const_mu`, `dL_net_nm`). The Lombardi post-processing does not enter them, and v1/v2 give identical constant-mobility values.
- **Invalid for halo devices in v1 directories:** the Lombardi quantities (`lombardi_mu`). See erratum E9-1 in `phase3/reports/C9_Research_Phase_Report_2026-10-31.md`.
- **Run records:**
  - `DIAG_RUNS.txt` covers D1–D3. Its first three lines are a tooling error (`/usr/bin/time` not installed; no simulation started) and are kept as a record.
  - `RUNS.txt` and `RUNNER_END.txt` in each directory come from `phase3/c4/jobrunner.py`.
  - Job lists: `jobs_D_level1.json`, `jobs_v2_level1.json`, `jobs_v2_aend.json`.
- **Analysis:** `phase3/c9/e9_posthoc.py`, which writes `results/phase3/C9/analysis/e9_posthoc.json`.
