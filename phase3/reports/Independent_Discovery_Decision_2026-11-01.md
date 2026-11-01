# Independent discovery decision: two-length differential extraction

## Decision

The strongest remaining candidate is a **two-length differential estimator for
zero-bias MOSFET extraction**. It is more promising than extending the
conditioned small-signal audit and more falsifiable than fitting another
semi-empirical κ(α∞) curve. The present result is a validated model-level
demonstration, not yet a claim of a new publication algorithm. Status:
`MORE_EVIDENCE_NEEDED`.

## Opportunity comparison

| Candidate | Closest limitation | Decisive test | Current decision |
|---|---|---|---|
| κ(α∞) end-layer law | Existing κ curve is tabulated and its proposed logarithmic branch was rejected in C8 | Derive and predict unseen bias/geometry states without refitting | Scientifically interesting, but theory and prior-art risk remain high |
| Conditioned small-signal gate | Compared outputs are not automatically the same physical quantity; the numerical QA protocol is already established | Equivalence-controlled finite-bias and precision study | Retain as QA only, not the main contribution |
| Two-length differential extraction | Single-length conductance/charge extraction is biased by gate-dependent end terms; multi-length practice exists but the exact cancellation for this paired conductance/charge estimator has not been demonstrated in the project | Predict mobility from two unseen lengths in a laterally graded junction and compare with the known constant-mobility parameter | **Chosen for focused investigation** |

The candidate is worth testing because it makes a sharp prediction: if the
zero-bias affine laws hold, the estimator must recover the bulk mobility while
the per-end correction may change sign and magnitude. A failure in a new
grading configuration would directly bound its scope.

## Mathematical statement

For each fixed gate bias under the zero-drain-bias equilibrium model, write

\[
G(L)=\frac{\sigma_\infty}{L-2\Delta L},\qquad
Q_e(L)=q\left(N_\infty L+2E\right).
\]

Here `L` is the metallurgical channel length in cm, `Q_e` is the integrated
electron sheet charge per unit width over `[0,L]`, and `E` is a per-end excess
charge term. Finite differences over two lengths remove the end terms:

\[
\frac{d(1/G)}{dL}=\frac{1}{\sigma_\infty},\qquad
\frac{dQ_e}{dL}=qN_\infty,
\]

so

\[
\mu_{2L}=\frac{[d(1/G)/dL]^{-1}}{dQ_e/dL}.
\]

This is an algebraic consequence of the model assumptions, not an independent
physical law. It does not establish validity at finite drain bias, under
transient conditions, or with quantum/Fermi–Dirac corrections.

## New falsifiable test

Before evaluating the new files, the prediction was fixed: for a resolved
source/drain deck with lateral erfc grading (`σ = 15 nm`), `N_A = 10^18 cm⁻³`,
`t_ox = 2 nm`, level 0, constant `μ = 1350 cm²/(V s)`, and `L = 1, 2 µm`,
the two-length estimator should recover 1350 cm²/(V s) across `V_G = 1.32,
1.60, 2.00 V`. The gate-dependent `ΔL` was not assumed to be positive or
length-independent in advance.

Six new DEVSIM runs completed with exit code 0. The direct per-end corrections
were `32.7806`, `−4.7846`, and `−33.9541 nm` at the three gate voltages. The
corrections agree between the two lengths to at most `8.71×10⁻⁶ nm`. The
two-length mobility errors were:

| `V_G` | `μ₂L` (cm²/V s) | relative error |
|---:|---:|---:|
| 1.32 V | 1349.9999765 | −1.74×10⁻⁸ |
| 1.60 V | 1349.9999929 | −5.27×10⁻⁹ |
| 2.00 V | 1349.9999977 | −1.67×10⁻⁹ |

This supports the estimator **within this zero-bias discrete model** and shows
that it remains valid when the end correction changes sign. It is not an
independent verification of the underlying device physics: the conductance and
charge come from the same DEVSIM equilibrium solution, and the reference
mobility is the configured constant parameter.

## Closest prior work and remaining gap

The project’s literature audit identifies multi-length split-C–V and effective
channel-length extraction as established practice: Guo–Chung–Hsu (IEEE T-ED
1994, DOI `10.1109/16.324592`), Romanjek et al. (SSE 2005, DOI
`10.1016/j.sse.2004.08.019`), and Taur (IEEE T-ED 2000, DOI
`10.1109/16.817582`). The available records are abstract/metadata level for
the first two, so no novelty claim is made from their absence of an identical
formula. The present result may become a contribution only if it is shown to
give a useful, independently measurable correction under models and operating
conditions where existing extraction procedures fail.

The strongest invalidating tests are finite drain bias, nonconstant mobility,
Fermi–Dirac/quantum statistics, and a separately implemented terminal-charge
calculation. The current phase does not claim those cases.

## Reproducibility and limits

- Script: `phase3/independent_discovery/two_length_differential.py`.
- Raw data: `results/phase3/independent_discovery/graded_probe/`.
- Machine-readable output: `results/phase3/independent_discovery/analysis/`.
- The attempted `L = 0.5 µm` run is preserved as a failed mesh assertion; no
  result is inferred from it.
- No production S1/S2/S3 solver, official golden, or historical result was
  modified.

Further investment is justified only as a bounded validation study: add one
independent terminal-charge path and one finite-drain-bias test before treating
this as a publication direction. If either breaks the cancellation materially,
the contribution should be reported as a zero-bias diagnostic rather than
expanded into a manuscript.
