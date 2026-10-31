"""C-9 erratum E9-1 regression test.

The Lombardi impurity input of probe9 must be the total ionised impurity N_A + N_D + N_halo, so the halo never
drives it to the 1 cm^-3 clip of endprobe.lombardi.  The net doping used for the equilibrium solve must stay
N_D - N_A - N_halo.  probe9 is imported in a subprocess because importing it installs its short-L lateral-map
wrapper on phase3.c6.wprobe, which must not leak into other tests.
"""
from __future__ import annotations

import os
import subprocess
import sys
import textwrap

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

CODE = textwrap.dedent('''
    import sys, types
    import numpy as np
    sys.path.insert(0, __ROOT__)
    from phase3.c9 import probe9 as P
    from phase3.c7 import endprobe as E7
    P.HALO.update(Nh=2e18, xh=20e-7, sx=15e-7, yh=15e-7, sy=15e-7)
    geo = types.SimpleNamespace(L_cm=1e-4, h_w_cm=0.5e-4, sigma_cm=0.0, na=1e18)
    xs = np.linspace(-0.3e-4, 1.3e-4, 3201)
    ys = np.full_like(xs, 15e-7)
    nd0 = P._DONOR0(geo, xs, ys)
    h = P.halo(geo, xs, ys)
    assert h.max() > 1.9e18, 'pocket missing'
    np.testing.assert_allclose(P.donor9(geo, xs, ys), nd0 - h)          # net doping (equilibrium): unchanged
    ntot = geo.na + P.donor_total9(geo, xs, ys)
    np.testing.assert_allclose(ntot, geo.na + nd0 + h)                  # Lombardi input: total ionised impurity
    assert ntot.min() >= geo.na
    nv1 = geo.na + P.donor9(geo, xs, ys)                                # the v1 input (erratum E9-1)
    assert nv1.min() < 1.0, 'v1 input should fall below the 1 cm^-3 clip inside the pocket core'
    E = np.full_like(xs, 5e5)
    k = int(np.argmax(h))
    assert E7.lombardi(E, ntot)[k] > 2.0 * E7.lombardi(E, nv1)[k]
    far = int(np.argmin(np.abs(xs - 0.5e-4)))                           # channel centre: no halo, v1 == v2
    assert E7.lombardi(E, ntot)[far] == E7.lombardi(E, nv1)[far]
    print('OK')
''')


def test_lombardi_impurity_input_is_total_ionised():
    pytest.importorskip('devsim')
    r = subprocess.run([sys.executable, '-c', CODE.replace('__ROOT__', repr(ROOT))], capture_output=True, text=True,
                       env=dict(os.environ), timeout=600)
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
    assert r.stdout.strip().endswith('OK')
