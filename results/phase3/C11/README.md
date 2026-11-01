# C-11 data

| File or directory | Contents |
|---|---|
| `c11_eval.json` | Pre-registered evaluation of H-lin and H-q2D on the 24 C-10 deep-subthreshold states: predictions, errors, verdicts, B_lin, B_fit, and far-tail B_∞ (device-1 far-tail fits are invalid; see the report) |
| `c11_eval_stdout.txt` | Console output of the evaluation |
| `indep/fd2d_VG*.json` | Independent finite-volume solver (`phase3/c11/fd2d_indep.py`) for device 2 at ψ̂ = 21 and 12 |
| `devsim_check/` | DEVSIM coupled small-signal check (`probe10.py --devsim`): state JSONs (key `devsim`), netcheck diagnostics, logs |
| `c11_figure.png/.pdf` | C-11 figure (`phase3/c11/c11_figures.py`) |

**Environments:**
- `c11_eval.json` and the figure were produced with the venv (Python 3.12, numpy 2.3.5).
- `indep/*.json` were produced with the system Python 3.13 (numpy 2.5.3). That solver does not need DEVSIM.
