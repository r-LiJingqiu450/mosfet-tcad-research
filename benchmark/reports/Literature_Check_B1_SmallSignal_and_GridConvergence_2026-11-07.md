# Literature check B1: the WP8 claims against the papers supplied so far

**Date:** 2026-11-07 · **Author:** Claude · **Step:** B of the pre-submission plan (literature), first part.

**Sources.** Six files supplied by the author on 2026-11-07. Three of them had already been read for the WP5 literature gate (2026-11-04, `Literature_Gate_LowCurrent_Benchmark`). This time they were read again against the WP8 claims. The PDFs are not stored in the repository.

**Result.**
- No source describes the WP8 mechanism, or the failure of a coupled small-signal solve with a clean residual. The novelty statement stays as it was, still limited to the literature reviewed.
- Two references were added to the manuscript, and four passages were sharpened.
- The most important source for the novelty question, Laux's paper on small-signal analysis (1985), is still unread.

## 1. Questions fixed before reading

- **B1-1.** Does any source describe one of the following?
  - (a) that the small-signal output of a coupled (ψ, n, p) Jacobian is limited by round-off at low current;
  - (b) that the potential-column entries of the continuity rows cancel where the potential response is uniform;
  - (c) that refining the solution cannot repair this;
  - (d) a remedy?
- **B1-2.** What do the sources say about the choice of unknowns, scaling, residual checks and linear solvers? The question is whether they bear on the WP8 statements about variables and remedies.
- **B1-3.** Is the observed order 1.98–2.22 consistent with what is known about the convergence of Scharfetter–Gummel currents?
- **B1-4.** Does the refinement method used for M1 have a citable source, and do its stated assumptions match the WP8 interpretation?
- **B1-5.** Are there sources for the weighting-function (Ramo–Shockley) terminal current in semiconductors beyond the two originals already cited?

## 2. What each source establishes

### S1. Selberherr, *Analysis and Simulation of Semiconductor Devices* (Springer, 1984)

Read from the book's front matter (copyright page with ISBN and DOI, contents) and from the 1983 Habilitationsschrift of the same title (TU Wien, 463 scanned pages with OCR text). Its section numbering matches the book's contents in every section checked (5.2, 5.5, 6.1, 7.1, 8.1, 8.2, 8.5, 8.6). The thesis was searched in full; sections 5.2, 5.5, 6.1, 7.1 and 8.1 were read in full.

- **§5.2, choice of unknowns.**
  - Slotboom variables need a dynamic range of "more than 32 orders of magnitude" for ψ ∈ [−1, 1] V, so they are limited to low-voltage applications.
  - Quasi-Fermi variables are all of the same order of magnitude and exclude the negative densities that round-off can produce with (ψ, n, p), but they make the continuity equations exponentially nonlinear.
  - Selberherr favours (ψ, n, p) for computations and recommends that set.
- **§5.5, scaling.** DeMari scaling and a "better" scaling. With the latter, the residuals of the scaled equations measure the computational accuracy.
- **§6.1.**
  - The Bernoulli function must be implemented with care to avoid underflow and overflow.
  - The contact densities should not be computed by the direct formula, "because of inherent problems with cancellation".
  - For a current-controlled contact, the contact current is the integral of the current density, by the midpoint rule.
- **§7.1.** Newton's iteration is terminated on the residual norm, O(10⁻¹⁰) in scaled units.
- **§8.1.**
  - Linear systems are factored without pivoting for stability.
  - Diagonal (row, or row and column) scaling "is sufficient for the semiconductor problem".
- **Small-signal analysis.** Not treated, apart from history references.
- **Bearing.**
  - DEVSIM's unknowns, ψ and n, are the computational set Selberherr recommends, so the WP8 finding concerns standard practice. The manuscript now says so in one sentence.
  - The book discusses round-off and cancellation in other places: negative densities, the Bernoulli function, contact densities. It does not discuss them for linear-response outputs.
  - It does not pre-empt B1-1.

### S2. Demmel, Hida, Kahan, Li, Mukherjee, Riedy, "Error bounds from extra-precise iterative refinement"

Read: the technical-report version (UC Berkeley, 12 Nov. 2004, 45 pp.). The published version is *ACM TOMS* 32(2), 325–351, 2006, doi:10.1145/1141885.1141894, checked against the Scilit record.

- **Method.** Iterative refinement with the residual computed in extra precision yields accurate solutions "to all but the most ill-conditioned" systems. For componentwise accuracy the solution is held in doubled precision. WP8's M1 is this design.
- **Assumption.** "The input data A and b are assumed to be stored exactly in the working precision"; *x* denotes the exact solution of the stored system.
- **Bearing.**
  - B1-4 is answered. The manuscript now cites the paper for the method.
  - It also adds the sentence that refinement returns the solution of the stored system, takes its entries as exact, and therefore cannot recover what their rounding has lost. That is exactly the WP8 registered interpretation.
  - No pre-emption: the paper is generic.

### S3. Ghione & Filicori, *IEEE TCAD* 12(3), 425–438 (1993)

Re-read in full.

- **Content.**
  - Adjoint (interreciprocity) and Green's-function approach to the small-signal sensitivity and noise of majority-carrier devices.
  - The small-signal system is the discretized Poisson–continuity system in an admittance-like nodal form.
  - The aim is computational efficiency.
  - Nothing is said about round-off, conditioning or the accuracy of small-signal outputs at low current.
- **Bearing.**
  - The adjoint λ used in WP8's first-order sensitivity is the device form of the adjoint method. The manuscript now cites the paper where λ is defined.
  - No pre-emption of B1-1.

### S4. Farrell & Peschka (preprint 2018, 24 pp.; journal version *CAMWA* 78 (2019) 3731–3747)

The preprint was supplied now. The journal version was read again from the WP5 copy.

- **Content.**
  - With smooth doping (κ = 500), the Voronoi finite-volume (Scharfetter–Gummel) total current converges quadratically, |I_h − I_2h| ∼ h². This is in both versions.
  - With near-discontinuous doping (κ = 5·10⁵), the finite-volume potential converges only linearly, and the finite-volume error is dominated by the lack of regularity of the doping. The preprint also says the finite-volume current then converges linearly.
  - Logarithmic boundary layers form in the minority-carrier quasi-Fermi potential at ohmic contacts once the current exceeds an intrinsic limiting current.
- **Bearing on B1-3.**
  - WP8's observed order of 1.98–2.22 matches their smooth-doping case, and the deck-R wells have smooth erfc profiles. The manuscript now says this.
  - Their boundary-layer mechanism does not apply here: holes are held at equilibrium, and the currents are linear-response currents at zero bias.

### S5. Jin, Park, Min, "An application of the Ramo–Shockley theorem for calculation of the terminal currents of a MOSFET by the particle simulation"

Two-page abstract, pp. 52–53. Re-read.

- **Content.**
  - The weighting potential f_j is defined with electrode j at unit potential, the others grounded, and the charge removed.
  - It is used to compute steady-state, transient and noise terminal currents of a 0.25 µm MOSFET in particle simulation.
- **Venue.** Venue and year are not printed. The references point to 1991, so it cannot be cited as it stands.
- **Bearing on B1-5.** It names two sources for the extension of Ramo–Shockley to semiconductor devices with fixed space charge. Neither is read yet:
  - P. De Visschere, *Solid-State Electron.* 33, 455 (1990);
  - H. S. Kim, H. S. Min, T. W. Tang, Y. J. Park, *Solid-State Electron.* (1991).

  If confirmed, they are the right citations for the weighting-function current in semiconductors, beside the two originals.

### S6. Springer front matter of S1 (12 pp.)

Used only for the bibliographic data:
- © 1984 Springer-Verlag Wien;
- ISBN 978-3-7091-8754-8, e-ISBN 978-3-7091-8752-4;
- DOI 10.1007/978-3-7091-8752-4;
- the contents, which match the thesis's sections.

## 3. Outcome per question

| Question | Outcome |
|---|---|
| B1-1 mechanism | **Not found** in S1–S5. The novelty statement stays "not found in the literature we reviewed". The manuscript now names what was reviewed. Laux (1985) is the decisive unread source. |
| B1-2 variables and practice | The (ψ, n) unknowns are the recommended set (S1). The claim stays specific to DEVSIM's Jacobian on these meshes; no claim is made about other variable sets. |
| B1-3 grid order | Consistent with S4's smooth-doping result. |
| B1-4 refinement | Cited (S2); its exact-data assumption is stated in the manuscript. |
| B1-5 weighting function | Two candidate sources identified through S5; not yet read. |

## 4. Manuscript changes (`paper/main.tex`, `paper/references.bib`)

1. **§3.3, grid convergence.** One sentence comparing the observed order with S4 (smooth versus near-discontinuous doping).
2. **§5.3, "The data, not the solver".**
   - Refinement is cited to S2 as well as Higham.
   - New sentence: refinement takes the stored entries as exact and cannot recover their rounding.
3. **§5.3, "Where the sensitivity lies".**
   - λ is identified as the device form of the adjoint method (S3).
   - DEVSIM's unknowns are identified as the set S1 recommends.
   - The remedy sentence now separates what is shown from what is inferred:
     - the reference is free of the floor, because its energy form is stationary in the potentials and, by Rayleigh's monotonicity law, moves by at most the relative rounding of the conductances;
     - a higher-precision assembly and solve would lower the floor in proportion (first-order theory, untested);
     - an edge-difference assembly may avoid it (untested).
4. **§7, relation to prior work.** The sources reviewed for the mechanism are named.
5. **`references.bib`.** `Selberherr1984` and `Demmel2006` were added, each with its verification record. The header now allows per-entry check dates.

The manuscript builds to 15 pages with 29 references, no undefined citations and no LaTeX warnings.

## 5. Still needed (from the author)

**Papers, in order of importance:**
1. S. E. Laux, "Techniques for small-signal analysis of semiconductor devices," *IEEE T-ED* 32(10), 1985 (also *IEEE TCAD* 4(4), 1985). It decides B1-1.
2. R. E. Bank, D. J. Rose, W. Fichtner, "Numerical methods for semiconductor device simulation," *IEEE T-ED* 30(9), 1983.
3. S. J. Polak, C. den Heijer, W. H. A. Schilders, P. Markowich, "Semiconductor device modelling from the numerical point of view," *Int. J. Numer. Methods Eng.* 24, 1987.
4. Small-signal and noise numerics after Ghione & Filicori (the Bonani–Ghione line of work on small-signal and noise analysis with physics-based models), if available.
5. The two Ramo–Shockley extensions named in S5: De Visschere, *Solid-State Electron.* 33 (1990); Kim, Min, Tang, Park, *Solid-State Electron.* (1991).

**Manuals** (the author has Silvaco access through nanoHUB): the small-signal AC analysis and numerical-precision sections of the Silvaco Atlas and, if available, Sentaurus Device user manuals. They would show whether commercial tools document this limit.
