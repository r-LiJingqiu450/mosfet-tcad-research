# Literature check B2: small-signal methods and device-simulation numerics

**Date:** 2026-11-07 · **Author:** Claude · **Step:** B of the pre-submission plan (literature), second part. It follows B1 (`Literature_Check_B1_SmallSignal_and_GridConvergence_2026-11-07.md`).

**Sources.** Four papers supplied by the author on 2026-11-07, each read in full or searched in full. The PDFs are not stored in the repository.

**Result.**
- The decisive source, Laux (1985), does not describe the WP8 finding. Laux establishes sinusoidal steady-state analysis, which is DEVSIM's AC analysis and, at zero frequency, WP8's coupled solve, as the accurate small-signal method. He discusses finite word length only for the transient method.
- The closest prior art is Polak et al. (1987), on the cancellation in the discrete current density where neighbouring densities are large and nearly equal. That is the nonlinear-residual counterpart of the WP8 Jacobian mechanism, not the mechanism itself.
- The novelty statement stays "not found in the literature we reviewed". The reviewed literature now includes the canonical small-signal paper.

## 1. Questions

- **B2-1.** Does any source describe a round-off limit on small-signal outputs computed from the linearized (ψ, n, p) system? The specific cases: at low current, through cancellation in the Jacobian rows, or as an error floor that refinement cannot remove.
- **B2-2.** What is the closest prior art for the cancellation theme, and how does WP8 relate to it?
- **B2-3.** How do the sources characterize the accuracy of sinusoidal steady-state analysis?

## 2. Sources

### S7. Laux, "Techniques for small-signal analysis of semiconductor devices"
*IEEE Trans. CAD* CAD-4(4), 472–481, Oct. 1985, doi:10.1109/TCAD.1985.1270145. All ten pages were read as page images.

**The three methods compared:**
- transient excitation with Fourier decomposition (FD);
- incremental charge partitioning (CP);
- sinusoidal steady-state analysis (S³A).

**How S³A works:**
- It linearizes the (V, n, p) equations about the dc operating point and solves [J + jωD] X̃ = B (eqs. 10–11).
- In FIELDAY it solves this by block SOR that reuses the LU factors of the Newton Jacobian, or by direct factorization at high frequency.
- At ω → 0 it is the coupled solve J δx = δb examined in WP2, WP4 and WP8.

**Accuracy statements:**
- S³A provides "accurate, rigorously correct results" and is "extremely accurate".
- "High accuracy is obtained because the device is formally linearized and calculations are done directly in the frequency domain."
- Finite word length appears only for FD. The current spike must not cause "a loss of current continuity due to the finite computer word length", and "a typical dc solution using FIELDAY maintains current continuity to ten decimal digits".

**Examples:**
- p-n diodes;
- a 1-µm n-MOSFET at V_D = 3–4 V and V_G = 1–4 V, mostly above threshold. Below threshold Laux notes only that the capacitance matrix is reciprocal.

**Bearing.**
- There is no description of a low-current round-off floor; on B2-1 the answer is no.
- WP8 qualifies Laux's accuracy statement for low currents in double precision. The linearization is exact, but its double-precision Jacobian does not determine G at low current.
- The manuscript now cites Laux where the coupled solve is introduced, and again in the novelty paragraphs.

### S8. Bank, Rose, Fichtner, "Numerical methods for semiconductor device simulation"
*IEEE T-ED* 30(9), 1031–1041, 1983, doi:10.1109/T-ED.1983.21257. Searched in full; §II-C, §III-C and §III-D read.

- **§II-C, choice of variables and scaling:**
  - with the variables (u, v, w), v = e^{−v}, w = e^{w}, the diagonal operators are "somewhat ill-conditioned" and "poorly scaled";
  - they propose symmetric diagonal scaling;
  - changes of variables "may lead to undesirable numerical difficulties (underflow, overflow, loss of precision)".
- **§III-C:** whether one interpolation of the current-density coefficient "leads to a more stable computation with respect to roundoff is not well understood".
- **Small-signal outputs:** not treated.
- **Bearing.** A standard reference for the numerics (variables, scaling, conditioning). No pre-emption.

### S9. Polak, den Heijer, Schilders, Markowich, "Semiconductor device modelling from the numerical point of view"
*Int. J. Numer. Methods Eng.* 24(4), 763–838, 1987, doi:10.1002/nme.1620240408. Searched in full; the sections on choice of unknowns, the discretization, terminal currents and matrix properties were read.

- **Choice of unknowns.** There is a trade-off between the exponential character of the unknowns and the nonlinearity of the equations, and they introduce the mixed-variable approach.
- **"A cancellation problem" (p. 791).**
  - Where p_i and p_{i+1} are "large and almost equal", the straightforward discrete current density "will result in cancellation problems and might even lead to divergence in the non-linear iteration"; they rewrite it with Bernoulli-type functions.
  - In two and three dimensions "the same problem is encountered" and the remedy "is not as simple".
- **Terminal currents.** These "can be obtained in a straightforward way".
- **Condition numbers.** Given for scaled diagonal blocks; avalanche can make the matrices "very ill-conditioned for whatever choice of variables".
- **Small-signal outputs:** not treated.
- **Bearing on B2-2.** This is the closest prior art for the cancellation theme. It concerns the evaluation of the flux (the nonlinear residual), not the linearized output. WP8's mechanism is the Jacobian counterpart: the potential derivatives of each row cancel when applied to the uniform response of the neutral drain. The manuscript now cites Polak et al. in the relation-to-prior-work paragraph and states the difference.

### S10. Bonani, Donati Guerrieri, Ghione, "Physics-based simulation techniques for small- and large-signal device noise analysis in RF applications"
*IEEE T-ED* 50(3), 633–644, 2003, doi:10.1109/TED.2003.810477. Read from the Politecnico di Torino repository copy; searched in full.

- **Content:**
  - a review of physics-based small-signal and large-signal noise simulation;
  - linearization about the working point, Green's functions, the generalized adjoint technique, and harmonic balance for large signal.
- **Not discussed:** round-off or the accuracy of the small-signal solve at low current.
- **Bearing.** A modern review of the small-signal and noise literature, cited beside Ghione & Filicori. No pre-emption.

## 3. Outcome

| Question | Outcome |
|---|---|
| B2-1 | **Not found.** Laux discusses word-length limits only for the transient method; the others do not treat small-signal outputs. |
| B2-2 | Closest: Polak et al., cancellation in the discrete current density (nonlinear residual). WP8 is its Jacobian counterpart. Palm & Van de Wiele (WP5) is the contact-flux counterpart. |
| B2-3 | Laux: accurate and rigorously correct, with no low-current caveat. WP8 adds that caveat for double-precision Jacobians at low current. |

**Novelty status of the WP8 mechanism.** Not found in the reviewed literature: the canonical small-signal paper, three reviews of device-simulation numerics, two small-signal and noise papers, and the refinement literature. It is a new observation and explanation within that scope, still stated as "not found in the literature we reviewed". No claim is made beyond that scope.

## 4. Manuscript changes

1. **§5.3, first paragraph.** The coupled solve is identified as the zero-frequency limit of sinusoidal steady-state analysis, with Laux's accuracy result cited.
2. **§5.3, novelty paragraph.** It now says what Laux does and does not discuss, and adds Bonani et al.
3. **§7, relation to prior work.**
   - Polak et al.'s current-density cancellation is named as the closest prior art and distinguished from the WP8 mechanism.
   - The reviewed sources are listed.
4. **`references.bib`.** `Laux1985`, `BankRoseFichtner1983`, `Polak1987` and `Bonani2003` were added, each with its verification record.

The manuscript builds to 15 pages with 33 references, no undefined citations and no LaTeX warnings.

## 5. Still open
- The two Ramo–Shockley extensions to semiconductors, De Visschere (1990) and Kim, Min, Tang, Park (1991), identified in B1. These are optional and would refine one citation.
- The AC-analysis and precision sections of a commercial simulator manual (Atlas; Sentaurus if available). These are optional; they would show whether commercial tools document this limit.
