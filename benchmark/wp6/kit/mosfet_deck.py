"""DEVSIM model of a 2-D n-channel MOSFET with resolved n+ source/drain wells.

Files read from this directory:
  device.json          device description (geometry, doping, physics, contacts, bias)
  mesh.json            node lines of the fixed tensor-product mesh
  zero_bias_state.npz  a converged zero-bias solution, used as the starting point of build()

Physics: Boltzmann statistics, constant electron mobility, holes held at equilibrium p = n_i exp(-psi/V_T) (not a
solution variable), no generation/recombination, T = 300 K. Unknowns: electrostatic potential psi (silicon and
oxide, referenced to the intrinsic level) and electron density n (silicon). Electron current: Scharfetter-Gummel edge
current. Units: cm, V, cm^-3; contact currents in A/cm (per unit width out of plane).

Device 'mos', regions 'si' and 'ox', interface 'ox_si', contacts:
  gate   oxide top over the channel: psi = V_G + psi_b0 (parameter 'gate_bias' = V_G)
  source silicon edge x = x_min, 0 <= y <= contact depth: psi = built-in value, n = n_i exp(psi/V_T)
  drain  silicon edge x = x_max, 0 <= y <= contact depth: psi = built-in value + V_D, n = n_i exp((psi - V_D)/V_T);
         V_D is the voltage of circuit node 'drain_bias', driven by the voltage source 'V1' (initially 0 V)
  body   silicon bottom y = 1.2 um: psi = psi_b0, n = n_i exp(psi/V_T)
  oxc_s, oxc_d  the two Si/SiO2 corner nodes on the oxide side, tied to the source/drain contact potential

Use:
  import mosfet_deck as M
  info = M.build()                                   # builds the device and converges V_D = 0
  devsim.circuit_alter(name='V1', value=1e-3)        # change the drain voltage
  M.solve_dc()                                       # Newton until the updates are below 1e-12 V
"""
from __future__ import annotations

import json
import math
import os

os.environ.setdefault('DEVSIM_MATH_LIBS', '/usr/lib/x86_64-linux-gnu/libblas.so.3:/usr/lib/x86_64-linux-gnu/liblapack.so.3')

import numpy as np  # noqa: E402
from scipy.special import erfc, erfcinv  # noqa: E402

import devsim as ds  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
Q = 1.602176634e-19                       # C
EPS0_CM = 8.8541878128e-12 * 1e-2         # F/cm
VT = 1.380649e-23 * 300.0 / Q             # V
NI = 1.0e10                               # cm^-3
MU_N = 1350.0                             # cm^2/(V s)
ND_PEAK = 1e20                            # cm^-3
DEV, SI, OX = 'mos', 'si', 'ox'

DEVICE = json.load(open(os.path.join(HERE, 'device.json')))
MESH = json.load(open(os.path.join(HERE, 'mesh.json')))
NA = float(DEVICE['doping']['N_A_cm3'])
P0 = .5 * (NA + math.sqrt(NA * NA + 4 * NI * NI))
PSI_B0 = -VT * math.log(P0 / NI)          # bulk potential (V)
VG = float(DEVICE['bias']['V_G_V'])


def _model(region, name, eq, wrt=()):
    ds.node_model(device=DEV, region=region, name=name, equation=eq)
    for v in wrt:
        ds.node_model(device=DEV, region=region, name=f'{name}:{v}', equation=f'simplify(diff({eq},{v}))')


def _edge(region, name, eq, wrt=()):
    ds.edge_model(device=DEV, region=region, name=name, equation=eq)
    for v in wrt:
        for n in ('n0', 'n1'):
            ds.edge_model(device=DEV, region=region, name=f'{name}:{v}@{n}', equation=f'simplify(diff({eq},{v}@{n}))')


def build_mesh():
    """Tensor-product mesh from mesh.json; each rectangle is split into two triangles."""
    xs = np.asarray(MESH['x_cm'], float); y_ox = np.asarray(MESH['y_oxide_cm'], float); y_si = np.asarray(MESH['y_si_cm'], float)
    ys = np.concatenate([y_ox, y_si[1:]]); kif = len(y_ox) - 1
    nx, ny = len(xs), len(ys)
    nid = lambda i, j: j * nx + i  # noqa: E731
    coords = []
    for j in range(ny):
        for i in range(nx):
            coords += [float(xs[i]), float(ys[j]), 0.0]
    names = [OX, SI, 'ox_si', 'gate', 'source', 'drain', 'body', 'oxc_s', 'oxc_d']
    P = {n: k for k, n in enumerate(names)}
    el = []
    for j in range(ny - 1):
        reg = P[OX] if j < kif else P[SI]
        for i in range(nx - 1):
            a, b, c, d = nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)
            el += [2, reg, a, b, c, 2, reg, a, c, d]
    g0, g1 = MESH['gate_columns']
    for i in range(nx - 1):
        el += [1, P['ox_si'], nid(i, kif), nid(i + 1, kif)]
    for i in range(g0, g1):
        el += [1, P['gate'], nid(i, 0), nid(i + 1, 0)]
    jw = kif + MESH['jw_si']
    for c, i in (('source', 0), ('drain', nx - 1)):
        for j in range(kif, jw):
            el += [1, P[c], nid(i, j), nid(i, j + 1)]
    for i in range(nx - 1):
        el += [1, P['body'], nid(i, ny - 1), nid(i + 1, ny - 1)]
    el += [1, P['oxc_s'], nid(0, kif), nid(0, kif)]
    el += [1, P['oxc_d'], nid(nx - 1, kif), nid(nx - 1, kif)]
    m = 'm7'
    ds.create_gmsh_mesh(mesh=m, coordinates=coords, elements=el, physical_names=names)
    ds.add_gmsh_region(gmsh_name=OX, mesh=m, region=OX, material='Ox')
    ds.add_gmsh_region(gmsh_name=SI, mesh=m, region=SI, material='Si')
    ds.add_gmsh_interface(gmsh_name='ox_si', mesh=m, name='ox_si', region0=OX, region1=SI)
    for c, r in (('gate', OX), ('source', SI), ('drain', SI), ('body', SI), ('oxc_s', OX), ('oxc_d', OX)):
        ds.add_gmsh_contact(gmsh_name=c, mesh=m, name=c, region=r, material='metal')
    ds.finalize_mesh(mesh=m)
    ds.create_device(mesh=m, device=DEV)
    return dict(nx=nx, ny=ny)


def build_physics(update_psi='log_damp', update_n='positive'):
    for r, eps in ((SI, 11.7), (OX, 3.9)):
        ds.set_parameter(device=DEV, region=r, name='Permittivity', value=eps * EPS0_CM)
        ds.set_parameter(device=DEV, region=r, name='ElectronCharge', value=Q)
    for name, val in (('n_i', NI), ('V_t', VT), ('mu_n', MU_N), ('NetDoping', -NA)):
        ds.set_parameter(device=DEV, region=SI, name=name, value=val)
    for name, val in (('gate_bias', 0.), ('source_bias', 0.), ('psi_b0', PSI_B0), ('s_couple', 0.)):
        ds.set_parameter(device=DEV, name=name, value=val)

    # solution variables; PsiLat holds the built-in contact potential (data, not solved)
    for r in (SI, OX):
        ds.node_solution(device=DEV, region=r, name='Potential')
        ds.edge_from_node_model(device=DEV, region=r, node_model='Potential')
        ds.set_node_value(device=DEV, region=r, name='Potential', value=PSI_B0)
    ds.node_solution(device=DEV, region=SI, name='Electrons')
    ds.edge_from_node_model(device=DEV, region=SI, node_model='Electrons')
    for r in (SI, OX):
        ds.node_solution(device=DEV, region=r, name='PsiLat')
        ds.set_node_value(device=DEV, region=r, name='PsiLat', value=PSI_B0)
    ds.set_node_values(device=DEV, region=SI, name='Electrons',
                       values=[NI * math.exp(PSI_B0 / VT)] * len(ds.get_node_model_values(device=DEV, region=SI, name='Potential')))

    # oxide Poisson
    _edge(OX, 'ElectricField', '(Potential@n0-Potential@n1)*EdgeInverseLength', ('Potential',))
    _edge(OX, 'DField', 'Permittivity*ElectricField', ('Potential',))
    ds.equation(device=DEV, region=OX, name='PotentialEquation', variable_name='Potential', edge_model='DField',
                variable_update=update_psi)

    # silicon Poisson, holes at equilibrium
    _model(SI, 'Holes', 'n_i*exp(-Potential/V_t)', ('Potential',))
    _model(SI, 'PotentialNodeCharge', '-ElectronCharge*kahan3(Holes, -Electrons, NetDoping)', ('Potential', 'Electrons'))
    _edge(SI, 'ElectricField', '(Potential@n0-Potential@n1)*EdgeInverseLength', ('Potential',))
    _edge(SI, 'DField', 'Permittivity*ElectricField', ('Potential',))
    ds.equation(device=DEV, region=SI, name='PotentialEquation', variable_name='Potential',
                node_model='PotentialNodeCharge', edge_model='DField', variable_update=update_psi)

    # electron continuity with the Scharfetter-Gummel edge current
    _edge(SI, 'vdiff', '(Potential@n0 - Potential@n1)/V_t', ())
    ds.edge_model(device=DEV, region=SI, name='vdiff:Potential@n0', equation='V_t^(-1)')
    ds.edge_model(device=DEV, region=SI, name='vdiff:Potential@n1', equation='-vdiff:Potential@n0')
    ds.edge_model(device=DEV, region=SI, name='Bern01', equation='B(vdiff)')
    ds.edge_model(device=DEV, region=SI, name='Bern01:Potential@n0', equation='dBdx(vdiff) * vdiff:Potential@n0')
    ds.edge_model(device=DEV, region=SI, name='Bern01:Potential@n1', equation='-Bern01:Potential@n0')
    jn = ('ElectronCharge*mu_n*EdgeInverseLength*V_t*'
          'kahan3(Electrons@n1*Bern01,  Electrons@n1*vdiff,  -Electrons@n0*Bern01)')
    _edge(SI, 'ElectronCurrent', jn, ('Electrons', 'Potential'))
    _model(SI, 'NCharge', '-ElectronCharge*Electrons', ('Electrons',))
    ds.node_model(device=DEV, region=SI, name='ElectronGeneration', equation='0')
    ds.equation(device=DEV, region=SI, name='ElectronContinuityEquation', variable_name='Electrons',
                time_node_model='NCharge', edge_model='ElectronCurrent', node_model='ElectronGeneration',
                variable_update=update_n)

    # Si/SiO2 interface: continuous potential
    ds.interface_model(device=DEV, interface='ox_si', name='continuousPotential', equation='Potential@r0-Potential@r1')
    ds.interface_model(device=DEV, interface='ox_si', name='continuousPotential:Potential@r0', equation='1')
    ds.interface_model(device=DEV, interface='ox_si', name='continuousPotential:Potential@r1', equation='-1')
    ds.interface_equation(device=DEV, interface='ox_si', name='PotentialEquation', interface_model='continuousPotential',
                          type='continuous')

    # contact charge (displacement) models
    for r in (SI, OX):
        ds.edge_model(device=DEV, region=r, name='contactcharge_edge', equation='DField')
        for n in ('n0', 'n1'):
            ds.edge_model(device=DEV, region=r, name=f'contactcharge_edge:Potential@{n}', equation=f'DField:Potential@{n}')

    # gate
    ds.contact_node_model(device=DEV, contact='gate', name='gate_psi', equation='Potential - (gate_bias + psi_b0)')
    ds.contact_node_model(device=DEV, contact='gate', name='gate_psi:Potential', equation='1')
    ds.contact_equation(device=DEV, contact='gate', name='PotentialEquation', node_model='gate_psi',
                        edge_charge_model='contactcharge_edge')
    # body: potential
    ds.contact_node_model(device=DEV, contact='body', name='body_psi', equation='Potential - psi_b0')
    ds.contact_node_model(device=DEV, contact='body', name='body_psi:Potential', equation='1')
    ds.contact_equation(device=DEV, contact='body', name='PotentialEquation', node_model='body_psi',
                        edge_charge_model='contactcharge_edge')

    # drain voltage source
    ds.circuit_element(name='V1', n1='drain_bias', n2='0', value=0.0, acreal=1.0, acimag=0.0)

    # source and drain: electron density
    for c, b in (('source', 'source_bias'), ('drain', 'drain_bias')):
        eq = f'Electrons - n_i*exp((Potential - {b})/V_t)'
        ds.contact_node_model(device=DEV, contact=c, name=f'{c}_n', equation=eq)
        ds.contact_node_model(device=DEV, contact=c, name=f'{c}_n:Electrons', equation='1')
        ds.contact_node_model(device=DEV, contact=c, name=f'{c}_n:Potential', equation=f'-n_i*exp((Potential - {b})/V_t)/V_t')
        kw = dict(device=DEV, contact=c, name='ElectronContinuityEquation', node_model=f'{c}_n',
                  edge_current_model='ElectronCurrent')
        if c == 'drain':
            ds.contact_node_model(device=DEV, contact=c, name=f'{c}_n:{b}', equation=f'n_i*exp((Potential - {b})/V_t)/V_t')
            kw['circuit_node'] = b
        ds.contact_equation(**kw)


def add_body_electron_contact():
    """Body: electron density n = n_i exp(psi/V_T)."""
    eq = 'Electrons - n_i*exp(Potential/V_t)'
    ds.contact_node_model(device=DEV, contact='body', name='body_n', equation=eq)
    ds.contact_node_model(device=DEV, contact='body', name='body_n:Electrons', equation='1')
    ds.contact_node_model(device=DEV, contact='body', name='body_n:Potential', equation='-n_i*exp(Potential/V_t)/V_t')
    ds.contact_equation(device=DEV, contact='body', name='ElectronContinuityEquation', node_model='body_n',
                        edge_current_model='ElectronCurrent')


def set_contact_potentials():
    """Source/drain (and the oxide corner nodes): psi = PsiLat + s_couple * bias, with s_couple = 1."""
    for c, b in (('source', 'source_bias'), ('drain', 'drain_bias')):
        ds.contact_node_model(device=DEV, contact=c, name=f'{c}_psi', equation=f'Potential - (PsiLat + s_couple*{b})')
        ds.contact_node_model(device=DEV, contact=c, name=f'{c}_psi:Potential', equation='1')
        kw = dict(device=DEV, contact=c, name='PotentialEquation', node_model=f'{c}_psi', edge_charge_model='contactcharge_edge')
        if c == 'drain':
            ds.contact_node_model(device=DEV, contact=c, name=f'{c}_psi:{b}', equation='-s_couple')
            kw['circuit_node'] = b
        ds.contact_equation(**kw)
    for c in ('oxc_s', 'oxc_d'):
        b = 'drain_bias' if c == 'oxc_d' else 'source_bias'
        ds.contact_node_model(device=DEV, contact=c, name=f'{c}_psi', equation=f'Potential - (PsiLat + s_couple*{b})')
        ds.contact_node_model(device=DEV, contact=c, name=f'{c}_psi:Potential', equation='1')
        kw = dict(device=DEV, contact=c, name='PotentialEquation', node_model=f'{c}_psi', edge_charge_model='contactcharge_edge')
        if c == 'oxc_d':
            ds.contact_node_model(device=DEV, contact=c, name=f'{c}_psi:{b}', equation='-s_couple')
            kw['circuit_node'] = b
        ds.contact_equation(**kw)


def donor(xs, ys):
    """Donor density N_D(x, y) = 1e20 f_x(x) f_y(y) (cm^-3) of the two wells (see device.json)."""
    dop = DEVICE['doping']
    sigma, sigma_y, h_w = dop['sigma_nm'] * 1e-7, dop['sigma_y_nm'] * 1e-7, dop['h_w_um'] * 1e-4
    L = MESH['x_junctions_cm'][1]
    s2 = math.sqrt(2.0) * sigma
    xc = -s2 * float(erfcinv(2.0 * NA / ND_PEAK))
    fx = 0.5 * erfc((xs - xc) / s2) + 0.5 * erfc(((L - xs) - xc) / s2)
    fy = 0.5 * erfc((ys - h_w) / (math.sqrt(2.0) * sigma_y))
    return ND_PEAK * fx * fy


def coords(region=SI):
    return (np.array(ds.get_node_model_values(device=DEV, region=region, name='x')),
            np.array(ds.get_node_model_values(device=DEV, region=region, name='y')))


def _solve_one():
    return ds.solve(type='dc', absolute_error=1e300, relative_error=1e300, maximum_iterations=1, solver_type='direct', info=True)


def _errors(info):
    it = info['iterations'][-1]
    out = dict(circuit_abs=it['circuit']['absolute_error'])
    for reg in it['devices'][0]['regions']:
        for e in reg['equations']:
            out[f"{reg['name']}:{e['name']}"] = dict(abs=e['absolute_error'], rel=e['relative_error'])
    return out


def solve_dc(tol=1e-12, max_it=100, extra=2):
    """Newton iterations until max |d psi| <= tol (V) and V_T * max |d n|/n <= tol, then `extra` more iterations.
    Returns the iteration history; raises RuntimeError if it does not converge."""
    hist = []
    for k in range(max_it):
        e = _errors(_solve_one())
        dpsi = max(e['si:PotentialEquation']['abs'], e['ox:PotentialEquation']['abs'])
        dphin = VT * e['si:ElectronContinuityEquation']['rel']
        hist.append(dict(dpsi=dpsi, dphin=dphin))
        if not (math.isfinite(dpsi) and math.isfinite(dphin)) or dpsi > 50:
            break
        if dpsi <= tol and dphin <= tol:
            for _ in range(extra):
                e = _errors(_solve_one())
                hist.append(dict(dpsi=max(e['si:PotentialEquation']['abs'], e['ox:PotentialEquation']['abs']),
                                 dphin=VT * e['si:ElectronContinuityEquation']['rel'], extra=True))
            return dict(converged=True, iterations=k + 1, history=hist)
    raise RuntimeError(f'Newton did not converge: {hist[-3:]}')


def build():
    """Build the device and converge the zero-bias state (V_G from device.json, V_D = V_S = V_B = 0)."""
    m = build_mesh()
    build_physics()
    ds.set_parameter(device=DEV, name='s_couple', value=1.0)
    add_body_electron_contact()
    xs, ys = coords(SI)
    xo, yo = coords(OX)
    net = donor(xs, ys) - NA
    ds.node_solution(device=DEV, region=SI, name='NetDop')
    ds.set_node_values(device=DEV, region=SI, name='NetDop', values=list(net))
    _model(SI, 'PotentialNodeCharge', '-ElectronCharge*kahan3(Holes, -Electrons, NetDop)', ('Potential', 'Electrons'))
    # built-in contact potential: V_T asinh(N_net / 2 n_i) on the source/drain contacts; the oxide corner nodes take the
    # largest of these values
    neutral = VT * np.arcsinh(net / (2 * NI))
    in_contact = ys <= MESH['contact_extent_cm'] * (1 + 1e-12)
    cs = ((np.abs(xs - xs.min()) < 1e-15) | (np.abs(xs - xs.max()) < 1e-15)) & in_contact
    psi_w = float(np.max(neutral[cs]))
    ox_corner = ((np.abs(xo - xs.min()) < 1e-15) | (np.abs(xo - xs.max()) < 1e-15)) & (np.abs(yo) < 1e-15)
    st = np.load(os.path.join(HERE, 'zero_bias_state.npz'))
    psi_si, psi_ox, n_si = st['psi_si'], st['psi_ox'], st['n_si']
    assert len(psi_si) == len(xs) and len(psi_ox) == len(xo)
    for r, v, mask, lat in ((SI, psi_si, cs, neutral), (OX, psi_ox, ox_corner, np.full(len(psi_ox), psi_w))):
        ds.set_node_values(device=DEV, region=r, name='Potential', values=list(v))
        ds.set_node_values(device=DEV, region=r, name='PsiLat', values=list(np.where(mask, lat, v)))
    ds.set_node_values(device=DEV, region=SI, name='Electrons', values=list(n_si))
    ds.set_parameter(device=DEV, name='gate_bias', value=VG)
    ds.circuit_alter(name='V1', value=0.)
    set_contact_potentials()
    r = solve_dc()
    return dict(mesh=m, n_si_nodes=int(len(xs)), V_G=VG, psi_b0=PSI_B0, newton_iterations=r['iterations'])
