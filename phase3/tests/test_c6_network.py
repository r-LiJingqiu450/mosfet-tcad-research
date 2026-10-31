"""Unit tests for the C-6 zero-bias conductance network solver (phase3/c6/wprobe.solve_net) and the
Rayleigh-monotonicity bounds used in the end-correction decomposition (pure SciPy; no DEVSIM state)."""
import numpy as np
import pytest


def grid(nx, ny, cx, cy):
    """Tensor grid; cx[j, i] = conductance of lateral edge (i,j)-(i+1,j), cy[j, i] = vertical edge (i,j)-(i,j+1)."""
    nid = lambda i, j: j * nx + i
    a, b, c = [], [], []
    for j in range(ny):
        for i in range(nx - 1):
            a.append(nid(i, j)); b.append(nid(i + 1, j)); c.append(cx[j, i])
    for j in range(ny - 1):
        for i in range(nx):
            a.append(nid(i, j)); b.append(nid(i, j + 1)); c.append(cy[j, i])
    src = np.array([nid(0, j) for j in range(ny)]); drn = np.array([nid(nx - 1, j) for j in range(ny)])
    return np.array(a), np.array(b), np.array(c, float), src, drn


def test_uniform_sheet_and_series():
    pytest.importorskip('devsim')
    from phase3.c6.wprobe import solve_net
    nx, ny = 11, 6
    cx = np.full((ny, nx - 1), 2.0); cy = np.full((ny - 1, nx), 3.0)
    a, b, c, s, d = grid(nx, ny, cx, cy)
    G, Is, _, res = solve_net(nx * ny, a, b, c, s, d)
    assert abs(G - ny * 2.0 / (nx - 1)) < 1e-12 and abs(G + Is) < 1e-12 and res < 1e-12


def test_rayleigh_bounds_random():
    pytest.importorskip('devsim')
    from phase3.c6.wprobe import solve_net
    rng = np.random.default_rng(3)
    nx, ny = 15, 9
    cx = 10 ** rng.uniform(-6, 6, (ny, nx - 1)); cy = 10 ** rng.uniform(-6, 6, (ny - 1, nx))
    a, b, c, s, d = grid(nx, ny, cx, cy)
    G, *_ = solve_net(nx * ny, a, b, c, s, d)
    G_slice = 1.0 / np.sum(1.0 / cx.sum(axis=0))          # columns shorted
    G_rows = np.sum(1.0 / np.sum(1.0 / cx, axis=1))        # vertical edges cut
    assert G_rows <= G * (1 + 1e-10) and G <= G_slice * (1 + 1e-10)
