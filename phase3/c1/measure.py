"""Measurement helpers for C-1(a)/(c): GCA conductance on 1COL widths, cross-section
SG edge-flux sums (DC and linearized AC), contact currents, lateral uniformity.

Sign and unit conventions (verified empirically in C-1(c), reported in c1c.json):
  * DEVSIM 2-D currents are per unit out-of-plane depth: A/cm (cm units).
    I' [A/um] = I [A/cm] * 1e-4;  G' [S/um] = G [S/cm] * 1e-4.
  * An edge model F_e (ElectronCurrent, A/cm^2) enters the FV balance of an edge
    through F_e * EdgeCouple_e (cm), i.e. A/cm per unit depth.
  * Cross-section current at x = x_c: I(x_c) = sum_e s_e F_e EdgeCouple_e over edges
    whose end-points straddle x_c, with s_e = +1 if x(n0) < x_c < x(n1), -1 otherwise.
    The positive direction (source -> drain or the reverse) is established in C-1(c)
    by comparison with the drain contact current.
"""
from __future__ import annotations

import numpy as np

import devsim as ds

from phase3.c1 import deck as D


def si_arrays():
    x, y = map(np.array, D.coords())
    return x, y


def column(xc, tol=1e-12):
    """Indices of Si nodes on the mesh column x = xc, sorted by depth."""
    x, y = si_arrays()
    idx = np.where(np.abs(x - xc) < tol)[0]
    return idx[np.argsort(y[idx])]


def dual_widths(yc):
    """1-D FV dual widths (= trapezoid weights incl. half cells at the ends)."""
    h = np.diff(yc)
    w = np.zeros_like(yc)
    w[:-1] += h / 2
    w[1:] += h / 2
    return w


def gca_conductance(L_cm, xc=None):
    """G_GCA,h = (q mu_n / L) * sum_i n0_i w_i on one Si column (S/cm).
    Uses the equilibrium electron density of that column (R04 5.2(D), D4)."""
    x, y = si_arrays()
    if xc is None:
        xc = x[np.argmin(np.abs(x - L_cm / 2))]
    idx = column(xc)
    n = np.array(ds.get_node_model_values(device=D.DEV, region=D.SI, name='Electrons'))[idx]
    w = dual_widths(y[idx])
    ns = float(np.sum(n * w))                       # cm^-2
    return dict(G_S_per_cm=D.Q * D.MU_N * ns / L_cm, Ns_FV_cm2=ns, column_x_cm=float(xc),
                n_nodes=int(len(idx)))


def lateral_uniformity(name='Potential'):
    """max over depth rows of (max - min across columns) of a node quantity."""
    x, y = si_arrays()
    v = np.array(ds.get_node_model_values(device=D.DEV, region=D.SI, name=name))
    worst = 0.0; worst_y = None
    for yy in np.unique(np.round(y, 15)):
        row = v[np.abs(y - yy) < 1e-15]
        d = float(row.max() - row.min())
        if name == 'Electrons':
            d = d / float(np.max(np.abs(row)))
        if d > worst:
            worst, worst_y = d, float(yy)
    return dict(quantity=name, max_spread=worst, at_depth_cm=worst_y,
                kind='relative' if name == 'Electrons' else 'absolute_V')


_EDGE = {}


def edges():
    if not _EDGE:
        ds.edge_from_node_model(device=D.DEV, region=D.SI, node_model='node_index')
        n0 = np.array(ds.get_edge_model_values(device=D.DEV, region=D.SI, name='node_index@n0')).astype(int)
        n1 = np.array(ds.get_edge_model_values(device=D.DEV, region=D.SI, name='node_index@n1')).astype(int)
        ec = np.array(ds.get_edge_model_values(device=D.DEV, region=D.SI, name='EdgeCouple'))
        _EDGE.update(n0=n0, n1=n1, ec=ec)
    return _EDGE


def crossing(xc):
    x, _ = si_arrays(); e = edges()
    a, b = x[e['n0']], x[e['n1']]
    m = ((a < xc) & (b > xc)) | ((b < xc) & (a > xc))
    s = np.where(a < xc, 1.0, -1.0)
    return np.where(m)[0], s[m]


def cross_section_dc(xc):
    idx, s = crossing(xc); e = edges()
    F = np.array(ds.get_edge_model_values(device=D.DEV, region=D.SI, name='ElectronCurrent'))
    contrib = s * F[idx] * e['ec'][idx]
    x, y = si_arrays()
    diag = np.abs(y[e['n0'][idx]] - y[e['n1'][idx]]) > 0
    return dict(I_A_per_cm=float(np.sum(contrib)), n_edges=int(len(idx)),
                n_diagonal_edges=int(np.sum(diag)),
                max_abs_EdgeCouple_diagonal=float(np.max(np.abs(e['ec'][idx][diag]))) if np.any(diag) else 0.0,
                diagonal_contribution=float(np.sum(contrib[diag])))


def cross_section_ac(xc):
    """Linearized AC electron current through x = xc:
    dF_e = sum_v [dF/dv@n0 * dv(n0) + dF/dv@n1 * dv(n1)], v in {Potential, Electrons},
    with dv the AC node solution (Real + j Imag). Conduction current only."""
    idx, s = crossing(xc); e = edges()
    out = 0j
    for v in ('Potential', 'Electrons'):
        dv = (np.array(ds.get_node_model_values(device=D.DEV, region=D.SI, name=f'{v}_Real'))
              + 1j * np.array(ds.get_node_model_values(device=D.DEV, region=D.SI, name=f'{v}_Imag')))
        for k, nn in (('n0', e['n0']), ('n1', e['n1'])):
            d = np.array(ds.get_edge_model_values(device=D.DEV, region=D.SI, name=f'ElectronCurrent:{v}@{k}'))
            out += np.sum(s * d[idx] * dv[nn[idx]] * e['ec'][idx])
    return out


def contact_currents():
    out = {}
    for c in ('source', 'drain'):
        out[f'{c}:electron'] = ds.get_contact_current(device=D.DEV, contact=c,
                                                      equation='ElectronContinuityEquation')
    for c in ('source', 'drain', 'body'):
        out[f'{c}:potential(displacement)'] = ds.get_contact_current(device=D.DEV, contact=c,
                                                                    equation='PotentialEquation') \
            if c != 'source' and c != 'drain' or _has_psi(c) else None
    out['gate:potential(displacement)'] = ds.get_contact_current(device=D.DEV, contact='gate',
                                                                equation='PotentialEquation')
    out['V1.I'] = ds.get_circuit_node_value(solution='dcop', node='V1.I')
    return out


def _has_psi(c):
    return 'PotentialEquation' in ds.get_contact_equation_list(device=D.DEV, contact=c)
