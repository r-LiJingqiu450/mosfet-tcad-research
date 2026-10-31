# C-9 halo states re-run with probe9 v2 (erratum E9-1)

**What is here:** the 55 halo states of the C-9 matrix, re-solved after the delivery audit found an error in the Lombardi impurity term for halo devices:
- 5e17 halo: L = 0.5 and 1 µm;
- 2e18 halo: L = 0.5, 1 and 2 µm;
- 11 V_G values each, level 0.

**Version:** `probe9-2026-10-31b`. Job list: `results/phase3/C9/jobs_v2_halo.json`.

**The error (probe9 v1):**
- The Lombardi term received N_A + N_D − N_halo. `endprobe.lombardi` clips it at 1 cm⁻³.
- Wherever N_halo > N_A + N_D, the C·N^λ term collapsed.

**The fix (v2):** the term receives the total ionised impurity, N_A + N_D + N_halo.

**What did not change:**
- The equilibrium solution, the terminal charges and every constant-mobility quantity are unchanged. The Lombardi model enters only as post-processing (C-6 P4).
- The constant-mobility values reproduce the v1 values, verified in `results/phase3/C9/analysis/e9_posthoc.json` and in the report.

**What is kept:** the original v1 states in `results/phase3/C9/runs/` and the original analysis `e9_result.json`, both unchanged as historical evidence.

**Corrected analysis:**
```
python phase3/c9/e9_analysis.py --halo-runs results/phase3/C9/runs_v2 --tag v2
```
This writes `e9_result_v2.json`.
