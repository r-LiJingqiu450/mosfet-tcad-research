# Phase 2B literature screening

Date: 2026-10-28  
Scope: gradual-channel approximation (GCA), drift-diffusion numerical verification, discretization uncertainty, and independent 2D/TCAD comparison.

This is a screening report, not a novelty claim. The repository contains a one-dimensional Python S1/S2/S3 workflow and no independent 2D device solver or measured data. The sources below define comparison methods and established model limits; they do not establish that this project has a new model or a publication-ready result.

## Sources suitable for verification

| ID | Source | Why it matters | Link |
|---|---|---|---|
| L1 | A. A. Gummel, “A self-consistent iterative scheme for one-dimensional steady state transistor calculations,” *IEEE Transactions on Electron Devices* 11(10), 455–465 (1964). DOI: `10.1109/T-ED.1964.1536099` | Original self-consistent Gummel iteration. It is the historical basis for treating the deprecated S2 Gummel path as a numerical method comparison rather than a production convergence claim. | [IEEE DOI](https://doi.org/10.1109/T-ED.1964.1536099) |
| L2 | D. L. Scharfetter and H. K. Gummel, “Large-signal analysis of a silicon Read diode oscillator,” *IEEE Transactions on Electron Devices* 16(1), 64–77 (1969). DOI: `10.1109/T-ED.1969.16754` | Source of the Scharfetter–Gummel current discretization used by S2. It supports a discretization reference, while leaving nonlinear iteration and mesh error to be tested separately. | [IEEE DOI](https://doi.org/10.1109/T-ED.1969.16754) |
| L3 | H. C. Pao and C. T. Sah, “Effects of diffusion current on characteristics of metal-oxide-semiconductor transistors,” *Solid-State Electronics* 9(10), 927–937 (1966). DOI: `10.1016/0038-1101(66)90013-6` | Foundational MOSFET drift/diffusion and gradual-channel analysis reference. It is an analytical/model baseline for separating GCA/model error from numerical error. | [Elsevier DOI](https://doi.org/10.1016/0038-1101(66)90013-6) |
| L4 | Y. Tsividis and C. McAndrew, *Operation and Modeling of the MOS Transistor*, 3rd ed., Oxford University Press (2011). ISBN 978-0-19-517015-9 | Standard derivations for charge-sheet/GCA assumptions, threshold and channel-charge conventions. It is appropriate for checking units, branch definitions, and applicability limits. | [Oxford University Press](https://global.oup.com/academic/product/operation-and-modeling-of-the-mos-transistor-9780195170159) |
| L5 | C. J. Roy, “Review of code and solution verification procedures for computational simulation,” *Journal of Computational Physics* 205(1), 131–156 (2005). DOI: `10.1016/j.jcp.2004.10.036` | Separates code verification, solution verification, and validation. This is the framework used here to keep residual, mesh, and independent-reference evidence distinct. | [Elsevier DOI](https://doi.org/10.1016/j.jcp.2004.10.036) |
| L6 | W. L. Oberkampf and T. G. Trucano, “Verification and validation in computational fluid dynamics,” *Progress in Aerospace Sciences* 38(3), 209–272 (2002). DOI: `10.1016/S0376-0421(02)00005-2` | General V&V terminology and evidence hierarchy. The CFD setting differs from semiconductor devices, so this is a methodological reference rather than a device model reference. | [Elsevier DOI](https://doi.org/10.1016/S0376-0421(02)00005-2) |
| L7 | I. B. Celik et al., “Procedure for estimation and reporting of uncertainty due to discretization in CFD applications,” *Journal of Fluids Engineering* 130(7), 078001 (2008). DOI: `10.1115/1.2960953` | Grid-convergence-index style reporting. The project can use the same discipline for S1/S2 grids, while reporting that the governing equations are different. | [ASME DOI](https://doi.org/10.1115/1.2960953) |
| L8 | S. Selberherr, *Analysis and Simulation of Semiconductor Devices*, Springer (1984). DOI: `10.1007/978-3-7091-8752-4` | Semiconductor device discretization and nonlinear-solver reference, useful for positioning 1D drift-diffusion verification and a future 2D comparison. | [Springer](https://doi.org/10.1007/978-3-7091-8752-4) |
| L9 | P. A. Markowich, C. A. Ringhofer, and C. Schmeiser, *Semiconductor Equations*, Springer (1990). DOI: `10.1007/978-3-7091-0650-2` | Mathematical reference for drift-diffusion/Poisson equations, boundary conditions, and equilibrium limits. | [Springer](https://doi.org/10.1007/978-3-7091-0650-2) |
| L10 | ASME, *Standard for Verification and Validation in Computational Fluid Dynamics and Heat Transfer*, V&V 20-2009 | Defines a formal validation-uncertainty framing. It is a standard, not a claim that a CFD acceptance threshold transfers directly to TCAD. | [ASME standard page](https://www.asme.org/codes-standards/find-codes-standards/v-v-20-standard-verification-validation-computational-fluid-dynamics-heat-transfer) |
| L11 | Berkeley BSIM Group, BSIM4 model documentation and release material | A practical compact-model comparison target. It is useful for checking whether a 1D GCA result is being compared with a compact model outside its documented scope. | [BSIM model page](https://bsim.berkeley.edu/models/bsim4/) |
| L12 | Synopsys, Sentaurus Device product/documentation portal | A possible independent commercial TCAD comparison route; access, licensing, mesh, and model settings must be documented before any comparison is called independent. | [Synopsys Sentaurus Device](https://www.synopsys.com/silicon/tcad/device-simulation/sentaurus-device.html) |

The execution environment could not query Crossref directly: the API request was rejected by the configured proxy with HTTP 403. DOI and publisher links are therefore recorded for human verification rather than presented as a fresh API retrieval. No reference was added solely from an unverified search snippet.

## Questions that can be tested with the current project

1. **GCA error decomposition.** For a fixed geometry and doping, sweep (V_G), (V_D), channel length, and inversion level. Compare the production GCA result with the existing independent integral/reference path, then report model error, spatial error, nonlinear residual error, and interpolation error separately. A useful result is an error surface with explicit validity boundaries, not a single aggregate percentage.

2. **Threshold and branch sensitivity.** Compare threshold extraction based on the documented branch criterion with the Pao–Sah/Tsividis charge-sheet convention. Report changes in (V_T), (I_D), and extracted subthreshold/saturation boundaries under grid refinement. This tests convention sensitivity and does not establish a new threshold theory.

3. **CLM attribution.** Hold the GCA equation and grid fixed, then compare no-CLM, legacy-CLM, and revised production branches. The existing Phase 1 ablation can be extended by checking whether CLM corrections are being combined with threshold or discretization errors that partly cancel. The result must show signed component errors, not only the total error.

4. **1D-to-2D transferability.** Build or obtain a separately configured 2D finite-volume/drift-diffusion reference for one geometry and a small bias set, or use a licensed TCAD tool with a frozen model deck. Compare potential, carrier density, current, and terminal charge profiles. A 2D comparison is required before making claims about short-channel, edge, fringing, or nonuniform-junction behavior.

5. **Numerical uncertainty budget.** Combine mesh refinement, nonlinear stopping residual, independent-reference discrepancy, and parameter perturbation into a reproducible uncertainty table. Roy/Oberkampf/Celik provide the reporting vocabulary; the numerical thresholds must be derived from this project's scales and not imported unchanged from CFD.

## Gaps and decision

The repository currently has no measured device data, no independent 2D solver, no commercial TCAD deck, and no restored historical S1/S2 golden files. The literature supports the proposed verification questions and the use of Gummel/Scharfetter–Gummel terminology; it does not support a claim of algorithmic originality or experimental accuracy.

Recommended Phase 2C entry condition: proceed only with a verification-focused study after candidate S1/S2 baselines are versioned and independently checked. Defer any claim about 2D/short-channel applicability until a separately configured 2D reference is available. Treat the deprecated Gummel path as a documented nonconvergence diagnostic throughout.

