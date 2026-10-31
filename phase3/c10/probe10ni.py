"""C-10 post-hoc exact-similarity check (NOT pre-registered; disclosed in the C-10 report).

Why: the pre-registered similarity image P10-1 matched t_hat, h_hat and psi_hat_bi. It differed from device 1 by
6-7.5 % (P10-1 REJECTED).

Post-hoc diagnosis: in deck W the junction face carries psi = psi_bi = V_T ln(N_D/n_i), so the electron density next
to the face is n_face = N_D. Relative to the channel doping this is N_D/N_A, equivalently exp(psi_hat_bi)(n_i/N_A)^2.
The electrons screen part of the junction step. n_face/N_A is therefore a further dimensionless parameter, which the
electron-free similarity statement omitted:
  - device 1: n_face/N_A = 1e4;
  - device 2' (psi_bi lowered by 0.119 V): n_face/N_A = 1.

Exact image of device 1 (1e16 cm^-3, n_i = 1e10, N_D = 1e20, 20 nm, h_w = 469 nm), with every dimensionless
parameter equal and only the mesh differing:
  - N_A = 1e18, t_ox = 2 nm, h_w = 47.3 nm;
  - n_i raised to 1e12 cm^-3. This keeps n_i/N_A = 1e-6 and psi_hat_b0 = ln(N_A/n_i).
  - The face value psi_bi = V_T ln(1e10) is unchanged, which corresponds to N_D = 1e22 = 1e4 N_A.
This is a mathematical similarity test of the simulator and of the corrected theorem, not a physical device.
Usage: as phase3/c7/endprobe.py (deck W), e.g. --deck W --L 1 --level 0 --vg 0.7835 --tox 2 --na 1e18 --hw 0.04727 --out <dir>
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from phase3.c1 import deck as D  # noqa: E402
from phase3.c10 import probe10 as P10  # noqa: E402  (installs nothing at import; robust solve_net installed below)
from phase3.c6 import wprobe as WP  # noqa: E402
from phase3.c7 import endprobe as E7  # noqa: E402

NI_IMAGE = 1.0e12


def main():
    assert '--deck' in sys.argv and sys.argv[sys.argv.index('--deck') + 1] == 'W', 'deck W only'
    D.NI = NI_IMAGE                      # used at call time by build_physics, set_bulk, onecol and the contacts
    WP.solve_net = P10.solve_net_robust  # E10-1 cancellation-free conductance
    E7.main()
    print('n_i image', D.NI, 'psi_bi (rel. intrinsic) unchanged', WP.PSI_BI)


if __name__ == '__main__':
    main()
