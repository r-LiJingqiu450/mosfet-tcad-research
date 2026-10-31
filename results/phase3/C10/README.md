# C-10 data

| Directory | Contents |
|---|---|
| `runs/` | 22 pre-registered states (endprobe, deck W), with job-runner records `RUNS.txt` and `RUNNER_END.txt` |
| `netcheck/` | Two states re-solved with `phase3/c10/probe10.py` (source-side and swapped-boundary currents, E10-1 check) and their `netcheck_*.json` diagnostics |
| `ni_image/` | Two post-hoc exact-similarity states (`probe10ni.py`; n_i = 1e12 cm⁻³; a mathematical test, not a device) |
| `tail/` | Far-field tail fits of the existing states (`tail_existing.json`) and a preview of the C-10 states |
| `analysis/` | `c10_result.json` (all verdicts, the corrected ΔL from the source-side current, the registered drain-side values `dL_drain_prereg`); `c10_figure.png/.pdf` |

Two design files sit at the top level:
- `design.json`: design numbers computed before any run;
- `jobs_*.json`: the job lists.

**Measured quantity.** ΔL is computed from the source-side network current, which is cancellation-free (E10-1). In 16 of the 22 runs, the endprobe field `result.const_mu.dL_net_nm` (from the drain-side current) violates the Rayleigh bounds and must not be used.

**Software environments:**
- `analysis/c10_result.json` and `c10_figure.*` were produced with the documented venv (Python 3.12, numpy 2.3.5).
- The system Python 3.13 (numpy 2.5.3) reproduces `c10_result.json` except for ≤ 2e-15 relative differences in five tail-fit fields (λ_fit, β). These come from LAPACK least squares.
- `tail/*.json` were produced with the system Python before the pre-registration commit. `c10_analysis.py` does not read them.
