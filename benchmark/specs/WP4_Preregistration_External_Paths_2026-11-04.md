# WP4 pre-registration: external paths, DEVSIM extended precision, NiTCAD feasibility

**Date:** 2026-11-04 (MST) · **Author:** Claude · **Branch:** `benchmark/wp4-external-2026-11-04`, from the WP2 release `2fe8859`
**Status: FROZEN on commit.** No WP4 run on a benchmark state has been made. Smoke runs on the non-production state SMOKE-sub (dev1, V_G = 0.55 V, level 0) are reported in §3.3 and §4.3. Every prediction informed by them is labelled as such. Changes after the first production run are allowed only as dated amendments appended at the end.

## 1. Context

**Plan v2, WP4:** "External paths: DEVSIM default contact current; NiTCAD Linux build feasibility."

**Tool self-check (2026-11-04), Branch A.** Across WP1, WP3 and WP2:
- no convergence failure, floating-point exception or mesh failure occurred;
- the only blocker was memory, resolved in amendment A9.

So no external tool is needed to proceed. The self-check also noted that every failure mode measured in WP2 is a float64 dynamic-range effect: the n⁺ contact edge conductances exceed G by up to 5e17. DEVSIM's documented answer to that class of problem is extended precision (`extended_solver`, `extended_model`, `extended_equation`). The installed DEVSIM 2.11.2 reports `extended_precision: True`, and every WP2 run used double precision.

WP4 therefore has four parts:
- **A.** Formal scoring of the DEVSIM default path (H3), from existing data.
- **B1.** Whether extended precision repairs DEVSIM's post-solve contact current.
- **B2.** Whether it repairs DEVSIM's built-in small-signal AC analysis.
- **C.** The NiTCAD verdict.

## 2. WP4-A: H3 with the DEVSIM default contact current (no new simulation)

**H3 (plan v1 §5):** "At least one external code path ... is scored on Tier A and Tier B. **Kill:** every external path passes every case trivially."

- **Path:** DEVSIM `get_contact_current` at the drain. The Richardson value of the ±0.1 and ±1 mV symmetric differences is taken from the WP2 P2 records (commit `454db70`).
- **Scored on:**
  - the 6 A-R and 5 HO three-terminal cases (same-mesh τ_num rule);
  - the B-R level-1 companions where a B0 P2 run exists: AR-1 to AR-6 (bracket rule).

**The outcome is already known** from the WP2 worked examples: 11/11 three-terminal cases fail. WP4-A only records it as the formal H3 assessment, together with the B-R counts. It is not a test.

## 3. WP4-B1: DEVSIM extended precision and the post-solve contact current

### 3.1 Runs

- **Code:** `benchmark/wp4/rdeck_dd_ext.py`. It sets the three DEVSIM flags before anything is built, then runs the **unchanged** WP2 P2 code (`benchmark/wp2/rdeck_dd.py`). Mesh, equilibrium procedure, Newton rule (τ = 1e-12 V, 2 extra iterations), biases (±0.1, ±1 mV) and current evaluations are the same as in WP2.
- **What extended precision changes** (DEVSIM source, commit `e3a25f0`, `src/math/ExternalPreconditioner.cc`): Newton residuals, models and equations are evaluated in float128, while every LU factorization and solve stays in double (UMFPACK). Newton thus works as mixed-precision iterative refinement.
- **States:** AR-1 to AR-6, body B1, **level 1**. Each run is paired with the WP2 double-precision P2 run on the same state and mesh. That makes 6 runs.

### 3.2 Observables

The circuit-node current I_C (−x[V1.I]) is the finite-bias reference current of the same run. Per bias:
- **drain error** e_D = |I_D − I_C|;
- **source error** e_S = |I_S + I_C + I_B|. This uses KCL, with DEVSIM's sign convention I_D = +I, I_S = −I.

Per run:
- the Richardson value G_R of each current, compared with the same-run network G_DD;
- the wall time.

### 3.3 Smoke evidence (SMOKE-sub, level 0; not a production state)

| Quantity | Double (WP2 smoke) | Extended |
|---|---|---|
| e_D, max over biases | 1.03e-11 A/cm | 6.7e-30 A/cm |
| e_S | 1.26e-11 A/cm | 6.4e-30 A/cm |
| DEVSIM drain-current G vs network at ±1 mV | +7.8 | +2.491e-4 (= the physical v²/6) |
| same, at ±0.1 mV | +242 | +2.491e-6 |
| final Newton updates | about 1e-16 V | about 1e-34 V |
| wall time | 18 s | 165 s |

### 3.4 Predictions

All three are informed by the single smoke state above, so they are weak evidence, as E-CIRC was in WP2.

| Id | Prediction | WP2 double precision, same runs |
|---|---|---|
| **X1** | In all 6 runs, max over the four biases of e_D ≤ 1e-20 A/cm, and of e_S ≤ 1e-20 A/cm | e_D = 3.3e-12 to 1.8e-11 A/cm; e_S = 2.4e-12 to 1.7e-11 A/cm |
| **X2** | In all 6 runs, the Richardson G_DD from the DEVSIM default drain current agrees with the same-run network G_DD within τ_num = 1e-6 | it fails τ_num in 6/6 |
| **X3** | E-CIRC holds: the circuit-current Richardson G agrees with the network within 1e-6 in all 6 runs | — |

X1 and X2 are falsifiable. The competing hypothesis is that a double-precision LU cannot refine a Jacobian whose κ_G reaches about 1e17 on these meshes (WP2, P3). The level-1 states are deeper in subthreshold (AR-1, ψ̂ = 12) and have larger meshes than the smoke state. **Each prediction is evaluated per run and reported as n/6. None is a kill.**

## 4. WP4-B2: DEVSIM's built-in small-signal AC analysis, double vs extended

### 4.1 Runs

- **Code:** `benchmark/wp4/ac_probe.py`. Same device, mesh and equilibrium; the drain source V1 carries acreal = 1 V.
- **Solve:** `devsim.solve(type='ac', frequency=1e-6 Hz)`. G_AC = −Re(V1.I), compared with the same-run network G_DD. At 1e-6 Hz the capacitive correction is below 1e-12 relative.
- **States:** AR-1 to AR-6, body B1, level 1, each in double and in extended precision. That makes 12 runs.

### 4.2 Why this path matters

It is what a DEVSIM user would call for a small-signal conductance. In extended mode the AC system is assembled in float128 but **factorized and solved once in double, with no refinement** (same source file). So, unlike Newton, it cannot profit from float128 residuals.

### 4.3 Smoke evidence (SMOKE-sub, level 0)

G_AC vs network: **+15.3 in double, +219 in extended.**

### 4.4 Predictions

| Id | Prediction | Basis |
|---|---|---|
| **Y1** | Double precision: \|G_AC/G_DD − 1\| > 1e-2 for AR-1, AR-2 and AR-3 | P-SS (`G_M` errors 1.7e3, 42, 0.31 at level 1) and the smoke state |
| **Y2** | Extended precision: \|G_AC/G_DD − 1\| > 1e-2 for AR-1, AR-2 and AR-3 as well | the code path (single double-precision solve) and the smoke state |

dev2 states (AR-4 to AR-6) are reported without a prediction; their `G_M` errors in WP2 were 25, 0.55 and 4.6e-3. None of this is a kill.

## 5. WP4-C: NiTCAD feasibility (verdict recorded before any build attempt)

**Source.** `github.com/nihaltodurkar-cms/nitcad` @ `95c2254` (2026-10-31), MIT licence. Its README lists among its capabilities "Accurate low currents (weighted terminal current with an error estimate)". That is the remedy this benchmark recommends, which makes NiTCAD the most relevant external code.

**Findings by inspection and syntax-only compilation (GCC 13.3, `-std=c++23`):**

| Item | Finding |
|---|---|
| Build system | Hard gates: `CMakeLists.txt` stops with "NiTCAD is Windows-only." and "NiTCAD requires MSVC"; uses vcpkg at a pinned baseline, `/W4 /WX`, CMake ≥ 4.3. The local CMake is 3.28. |
| Core code | 74 of 82 non-GUI, non-renderer translation units compile syntax-only. Of the other 8: 1 is missing `#include <algorithm>`; 5 use Windows APIs (`<intrin.h>`, `LoadLibrary` for the optional MKL, `RtlGetVersion` for provenance, a console handler, one more in `gate_tunnels.cpp`); 2 fail an Eigen version assertion. The renderer (Direct3D 12 / Direct2D, used by `nitcad plot`) and the GUI are Windows-only and not needed. |
| **Blocker** | `eigen_sparse_lu.cpp` and `eigen_ilut.cpp` `static_assert` **Eigen 5.0.1**. The comment states that the failure handling depends on Eigen 5.0.1 internals and that "a different Eigen must be re-checked against these before this assertion is changed". Eigen 5.0.1 is distributed from `gitlab.com`, which this environment cannot reach. PyPI offers only Eigen 3.4.1 (`cmeel-eigen`). |

**Decision.** Building NiTCAD against Eigen 3.4.1 would override the author's own verification gate. A result from such a build could not be presented as NiTCAD's, so **NiTCAD is not built in WP4** unless the official source archive is supplied:
- `eigen-5.0.1.tar.gz` from `https://gitlab.com/libeigen/eigen/-/archive/5.0.1/eigen-5.0.1.tar.gz`;
- SHA-512 `b337d3bc38440db190a8f1fbc4eabc0098e69fcc95bfba195fe039ffb942cae2a7f0153f3094f35fa26325750d1c62e20cccaf916a41f5c7f248ec5e5d30a942`, as recorded by vcpkg's `eigen3` port, version 5.0.1.

If it is supplied, the plan is:
1. Port the CLI to Linux with platform-glue patches only, listed file by file. No numerical code is changed.
2. Validate the port with NiTCAD's own Catch2 test suite.
3. Write a **separate dated pre-registration** of the deck mapping before any NiTCAD result on a benchmark case is computed.

## 6. Execution and integrity

- **Queue:** `benchmark/wp4/run_queue_wp4.py`, with the WP2 queue design:
  - two workers, memory budget 5.4 GB, single-threaded BLAS, `MALLOC_MMAP_THRESHOLD_=131072` (A9);
  - every attempt is appended to `benchmark/wp4/results/ledger.jsonl`, and the logs of earlier attempts are kept.
- **Time limit:** a run still going after 4 h is stopped and reported as not completed.
- **Evaluation:** X1–X3 and Y1–Y2 are evaluated exactly as written, by `benchmark/wp4/evaluate_wp4.py`.
- **Unchanged:** the WP2 code paths, the WP2 results, and WP1 and WP2 case files. New benchmark artefacts (for example the atlas entries) are added only after the evaluation, as a new release.
