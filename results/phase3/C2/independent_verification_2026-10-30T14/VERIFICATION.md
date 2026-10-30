# Independent re-execution of the C-2 package (2026-10-30, ~14:15 UTC)

**Why this record exists.** Commit `fe44d68` and the three untracked reports were produced in this container between 12:55 and 13:27 UTC. They came from an earlier run of the same authorization whose conversation was not available to the verifying run after a session restart. Before committing and archiving anything, the decisive computations were re-executed from the code at `fe44d68` into this directory.

| Re-executed (DEVSIM 2.11.2, venv 3.12.3, NumPy 2.3.5, SciPy 1.17.0) | Exit | Numbers compared with the recorded file | Differences |
|---|---|---|---|
| `corner_study.py --config A0 --level 0` | 0 | 482 | 0 |
| `corner_study.py --config A1 --level 0` | 0 | 244 | 0 |
| `corner_study.py --config X0 --level 0` | 0 | 301 | 0 |
| `corner_study.py --config X1 --level 0` | 0 | 223 | 0 |
| `corner_study.py --config X1 --level 1` | 0 | 223 | 0 |
| `corner_study.py --config X2 --level 0` | 0 | 235 | 0 |
| `c2/onecol_deck.py --form P --level 0` | 0 | 2354 | 0 |
| `c2/onecol_deck.py --form P --level 3` | 0 | 4730 | 0 |
| `c2/onecol_deck.py --form F --level 0` | 0 | 2354 | 0 |
| `c2/onecol_deck.py --form F --level 3` | 0 | 4730 | 0 |
| `c2/onecol_deck.py --form P --level 3 --points [[V_T,long, 0]] --tag vtlong` | 0 | 72 | 0 |

Timestamps, paths and git fields were excluded from the comparison.

**Direct cross-check against M0 (`phase3/m0.py`), outside the analysis script.**

| Point | ψ_s difference, L3 − M0 | N_s relative difference |
|---|---|---|
| V_G = 2.0 V, φ_n = 0 | −1.219e-5 V | +1.391e-5 |
| V_G = 1.25 V, φ_n = 0.5 V | −8.26e-7 V | +2.428e-5 |
| V_G = 0.8 V, φ_n = 0.05 V | −6.37e-7 V | +1.490e-5 |

These are identical to `t0t2.csv`.

**Code review of the corner repair** (`deck.build_mesh_explicit(corner='tie')`, `set_lateral_dirichlet`):
- The one-node oxide contact sits on the triple-point node.
- Its row is ψ − (PsiLat + s_couple·bias), with the same expression and data as the Dirichlet Si twin. It is therefore the discrete continuity condition, and it remains consistent when the lateral data follow V_D.
- Its displacement flux is assigned to the same terminal (`oxc_d` on circuit node `drain_bias`).
- No other oxide node is constrained.

**Regression:** see `../regression_final/`. `pytest semiconductor_sim phase3/tests` gives 142 passed, 4 skipped, exit 0 in both environments. Golden SHA is unchanged, and the production diff vs `7073f20` is empty.
