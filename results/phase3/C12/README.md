# C-12 data

Lightweight checks inside the Liang–Taur double-gate model. **No device simulation.**

| File | Contents |
|---|---|
| `lt_endcorr.json` / `lt_endcorr_stdout.txt` | `phase3/c12/lt_endcorr.py`, checks K1–K7 and K9 for 5 geometries × 3 barriers (see list below) |
| `lt_vt_rolloff.json` / `lt_vt_rolloff_stdout.txt` | `phase3/c12/lt_vt_rolloff.py`, check K8: constant-current V_t roll-off of the full LT integral against L06 eqs. (3.7) and (3.15) and against the end-correction term, plus rows/slice ratios. t_si = 10 nm, t_i = 1.5 nm, barrier 0.38 V, V_ds = 0.05 and 1 V, L = 20–1000 nm |

**Checks in `lt_endcorr` (K1–K7, K9):**
- λ₁ against the quoted values;
- conjugate-function identity;
- b₁ closed form against the projection and against quadrature;
- ΔL (slice, rows, Ein), C(s) and the thin-film limit;
- finite-L convergence;
- full-series truncation;
- sensitivity to b₁;
- a rough estimate of the source-edge error, read from a figure (K9).

**Environment.** Python 3.12.3 venv: numpy 2.3.5, scipy 1.17.0. No DEVSIM needed. Runtime about 3 s in total.

**Reproduce:**
```bash
python -m phase3.c12.lt_endcorr
python -m phase3.c12.lt_vt_rolloff
python -m pytest phase3/tests/test_c12_lt_identities.py -q
```

**Source documents read** (user attachments, not committed; SHA-256):

| Document | SHA-256 |
|---|---|
| Liang & Taur, IEEE T-ED 51(8) 2004 | `6bd74a36bff32d38de2c11ad87b03dd44dc99390a925baa940e72cdd66132329` |
| Liang, Ph.D. dissertation, UCSD 2006 (eScholarship qt49f030x3) | `b8e99dd8f83b75489caaa9efb47f6598fce36a4adb77a63b6217f77bf2cf767e` |
| Frank, Taur & Wong, IEEE EDL 19(10) 1998 | `04fa4d780f7ab3d2bbccf89db9d818821a66e8cfbbcdfc6292b60816e38f0d66` |
| Oh, Monroe & Hergenrother, IEEE EDL 21(9) 2000 | `3b42a897fefd16fd6ccfb2fafee6b1117190776571e6a777d48c90bd1bdc28eb` |
| Taur, IEEE T-ED 48(12) 2001 | `a79f05f6666e8ee53dc9acb7887a65830a6225523c73f8345d813201d41ef4f5` |

**Caveats:**
- Everything here is a statement about the LT model: the analytic potential, eq. (16), and its boundary idealisations. **It is not physical validation.**
- K9 and the figure readings are estimates by eye.
