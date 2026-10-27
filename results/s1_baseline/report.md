# S1 baseline experiment

This run solves the intrinsic one-dimensional Poisson–Boltzmann boundary-value problem with a damped Newton method. The surface potential is prescribed at the semiconductor interface and the far boundary is neutral. The reported surface charge uses a second-order one-sided derivative.

- Commit: `c64c463ba35fb3c2f177c9558217a20387046795`
- Grids: 101, 201, and 401 points over 20 Debye lengths
- Surface-potential sweep: -0.20 V to +0.20 V in 0.05 V increments
- Converged rows: 27/27
- Maximum dimensionless residual: `1.670e-12`
- Odd-symmetry absolute charge mismatch: `1.006e-21 C m^-2`

The grid sweep is a numerical convergence check, not an independent physical validation. The model is intrinsic and omits oxide capacitance, fixed charge, doping, mobility, channel current, CLM, and MOSFET/GCA behavior. S3 cannot be reproduced from the current repository because no S3 implementation or reference data is present.
