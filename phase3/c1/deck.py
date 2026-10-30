"""Minimal 2-D DEVSIM deck for C-1 feasibility checks (NOT a T0/T1/T2 deck).

Physics: Research_05 v2.1 section 0.2 (PM checklist), family F0.
  T = 300 K, V_T = kT/q from common/constants.py, n_i = 1e10 cm^-3,
  eps_Si = 11.7 eps0, eps_ox = 3.9 eps0, Boltzmann, mu_n = 1350 cm^2/Vs const,
  no generation/recombination, holes frozen p = n_i exp(-psi/V_T) (node model,
  not a solution variable), psi referenced to the intrinsic level,
  gate contact psi_gate = V_G + psi_b0 (V_FB = 0).

Units: DEVSIM in cm. Densities cm^-3, potentials V, contact currents A/cm
(per unit out-of-plane depth). I' [A/um] = I_DEVSIM [A/cm] * 1e-4 (R04 3.2).

Coordinates (DEVSIM): x along the channel (0 = source, L = drain);
y = depth (y < 0 oxide, y = 0 Si/SiO2 interface, y = x_max bottom).
This is the R04 (y, x) pair with axes swapped.

Contacts:
  gate   (oxide top)       : Potential = gate_bias + psi_b0
  body   (Si bottom,       : Potential = psi_b0; NO electron contact equation
          corners excluded)  -> zero normal electron flux (natural BC)
  source (Si x = 0)        : Electrons = n_i exp((Potential - source_bias)/V_t)
  drain  (Si x = L)        : Electrons = n_i exp((Potential - drain_bias)/V_t),
                             integrated into circuit node 'drain_bias'
  R-I(a) mode adds, on source and drain:
                             Potential = PsiLat + s_couple * <bias>
  where PsiLat holds lateral psi data (set from the V_D = 0 Neumann solution,
  i.e. the laterally uniform 1-column solution on identical y-lines), and
  s_couple is an artificial psi-data sensitivity used for the decoupling check.
  R-I(b) mode: no Potential contact equation on source/drain (Neumann psi).
"""
from __future__ import annotations

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import devsim as ds

import phase3  # noqa: F401
from common import constants as C

EPS0_CM = C.EPS0 * 1e-2          # F/cm
Q = C.Q_E
VT = C.VT
NI = 1.0e10                      # cm^-3
NA = 1.0e16                      # cm^-3
MU_N = 1350.0                    # cm^2/Vs
P0 = .5 * (NA + math.sqrt(NA * NA + 4 * NI * NI))
PSI_B0 = -VT * math.log(P0 / NI)

DEV = 'mos'
CORNER = {'on': False, 'mode': 'legacy_pin'}
SI = 'si'
OX = 'ox'


def pm_checklist():
    return dict(T_K=C.TEMP, V_T_V=VT, n_i_cm3=NI, eps_Si_rel=11.7, eps_ox_rel=3.9,
                statistics='Boltzmann', mu_n_cm2_Vs=MU_N, mobility='constant',
                generation_recombination='none', holes='frozen p=n_i*exp(-psi/V_T)',
                interface_charge='none', NA_cm3=NA, psi_b0_V=PSI_B0,
                gate='psi_gate = V_G + psi_b0 (V_FB = 0)', units='cm; I [A/cm]; I\'[A/um]=I*1e-4')


def _model(region, name, eq, wrt=()):
    ds.node_model(device=DEV, region=region, name=name, equation=eq)
    for v in wrt:
        ds.node_model(device=DEV, region=region, name=f'{name}:{v}',
                      equation=f'simplify(diff({eq},{v}))')


def _edge(region, name, eq, wrt=()):
    ds.edge_model(device=DEV, region=region, name=name, equation=eq)
    for v in wrt:
        for n in ('n0', 'n1'):
            ds.edge_model(device=DEV, region=region, name=f'{name}:{v}@{n}',
                          equation=f'simplify(diff({eq},{v}@{n}))')


def build_mesh(L_um=1.0, xmax_um=1.2, tox_nm=20.0, nx=21, y_ps0_nm=0.5, y_psmax_nm=20.0,
               ox_cells=5, corner_contacts=False, ox_corner_depth_cm=3e-8):
    """Tensor-product mesh from DEVSIM's built-in 2-D mesher (right triangles)."""
    L, xm, tox = L_um * 1e-4, xmax_um * 1e-4, tox_nm * 1e-7
    m = 'm'
    ds.create_2d_mesh(mesh=m)
    dx = L / (nx - 1)
    ds.add_2d_mesh_line(mesh=m, dir='x', pos=0., ps=dx)
    ds.add_2d_mesh_line(mesh=m, dir='x', pos=L, ps=dx)
    ds.add_2d_mesh_line(mesh=m, dir='y', pos=-tox, ps=tox / ox_cells)
    ds.add_2d_mesh_line(mesh=m, dir='y', pos=0., ps=y_ps0_nm * 1e-7)
    ds.add_2d_mesh_line(mesh=m, dir='y', pos=xm, ps=y_psmax_nm * 1e-7)
    # Thin inert padding ('gas', no equations) outside every contact boundary.
    # DEVSIM's 2-D mesher detects contacts/interfaces as region boundaries
    # (pattern from devsim testing/mesh2d.py and examples/diode/diode_common.py).
    pad = 1e-8
    ds.add_2d_mesh_line(mesh=m, dir='x', pos=-pad, ps=pad)
    ds.add_2d_mesh_line(mesh=m, dir='x', pos=L + pad, ps=pad)
    ds.add_2d_mesh_line(mesh=m, dir='y', pos=-tox - pad, ps=pad)
    ds.add_2d_mesh_line(mesh=m, dir='y', pos=xm + pad, ps=pad)
    ds.add_2d_region(mesh=m, material='gas', region='pad_l', xl=-pad, xh=0, yl=-tox - pad, yh=xm + pad)
    ds.add_2d_region(mesh=m, material='gas', region='pad_r', xl=L, xh=L + pad, yl=-tox - pad, yh=xm + pad)
    ds.add_2d_region(mesh=m, material='gas', region='pad_t', xl=0, xh=L, yl=-tox - pad, yh=-tox)
    ds.add_2d_region(mesh=m, material='gas', region='pad_b', xl=0, xh=L, yl=xm, yh=xm + pad)
    ds.add_2d_region(mesh=m, material='Ox', region=OX, xl=0, xh=L, yl=-tox, yh=0.)
    ds.add_2d_region(mesh=m, material='Si', region=SI, xl=0, xh=L, yl=0., yh=xm)
    ds.add_2d_interface(mesh=m, name='ox_si', region0=OX, region1=SI, xl=0, xh=L, yl=0, yh=0, bloat=1e-10)
    ds.add_2d_contact(mesh=m, name='gate', material='metal', region=OX, xl=0, xh=L, yl=-tox, yh=-tox, bloat=1e-10)
    ds.add_2d_contact(mesh=m, name='source', material='metal', region=SI, xl=0, xh=0, yl=0, yh=xm, bloat=1e-10)
    ds.add_2d_contact(mesh=m, name='drain', material='metal', region=SI, xl=L, xh=L, yl=0, yh=xm, bloat=1e-10)
    eps = .5 * dx
    ds.add_2d_contact(mesh=m, name='body', material='metal', region=SI, xl=eps, xh=L - eps, yl=xm, yh=xm, bloat=1e-10)
    if corner_contacts:
        # DEVSIM drops Si/SiO2 interface pairs at nodes that are also contact nodes
        # (mesh log: 'nx-2 node pairs'). The two Si corner nodes then lose the gate
        # coupling. Repair for R-I(a): pin the oxide corner nodes to the same lateral
        # data with one-node oxide contacts (oxide lateral edges otherwise Neumann).
        # A one-node box yields no contact edge, so the shortest oxide segment at each
        # corner is pinned (2 oxide nodes; 0.3 nm box covers one 0.24 nm segment).
        d = ox_corner_depth_cm
        ds.add_2d_contact(mesh=m, name='oxc_s', material='metal', region=OX, xl=0, xh=0, yl=-d, yh=0, bloat=1e-10)
        ds.add_2d_contact(mesh=m, name='oxc_d', material='metal', region=OX, xl=L, xh=L, yl=-d, yh=0, bloat=1e-10)
    CORNER['on'] = corner_contacts
    CORNER['mode'] = 'legacy_pin'
    ds.finalize_mesh(mesh=m)
    ds.create_device(mesh=m, device=DEV)
    return dict(L_cm=L, xmax_cm=xm, tox_cm=tox, nx=nx, y_ps0_nm=y_ps0_nm,
                y_psmax_nm=y_psmax_nm, ox_cells=ox_cells,
                corner_contacts=corner_contacts)


def build_physics(update_psi='log_damp', update_n='positive'):
    for r, eps in ((SI, 11.7), (OX, 3.9)):
        ds.set_parameter(device=DEV, region=r, name='Permittivity', value=eps * EPS0_CM)
        ds.set_parameter(device=DEV, region=r, name='ElectronCharge', value=Q)
    for name, val in (('n_i', NI), ('V_t', VT), ('mu_n', MU_N), ('NetDoping', -NA)):
        ds.set_parameter(device=DEV, region=SI, name=name, value=val)
    for name, val in (('gate_bias', 0.), ('source_bias', 0.), ('psi_b0', PSI_B0),
                      ('s_couple', 0.)):
        ds.set_parameter(device=DEV, name=name, value=val)

    # ---- solution variables
    for r in (SI, OX):
        ds.node_solution(device=DEV, region=r, name='Potential')
        ds.edge_from_node_model(device=DEV, region=r, node_model='Potential')
        ds.set_node_value(device=DEV, region=r, name='Potential', value=PSI_B0)
    ds.node_solution(device=DEV, region=SI, name='Electrons')
    ds.edge_from_node_model(device=DEV, region=SI, node_model='Electrons')
    for r in (SI, OX):                                         # lateral psi data (not solved)
        ds.node_solution(device=DEV, region=r, name='PsiLat')
        ds.set_node_value(device=DEV, region=r, name='PsiLat', value=PSI_B0)
    ds.set_node_values(device=DEV, region=SI, name='Electrons',
                       values=[NI * math.exp(PSI_B0 / VT)] * len(
                           ds.get_node_model_values(device=DEV, region=SI, name='Potential')))

    # ---- oxide Poisson
    _edge(OX, 'ElectricField', '(Potential@n0-Potential@n1)*EdgeInverseLength', ('Potential',))
    _edge(OX, 'DField', 'Permittivity*ElectricField', ('Potential',))
    ds.equation(device=DEV, region=OX, name='PotentialEquation', variable_name='Potential',
                edge_model='DField', variable_update=update_psi)

    # ---- silicon Poisson with frozen holes
    _model(SI, 'Holes', 'n_i*exp(-Potential/V_t)', ('Potential',))
    _model(SI, 'PotentialNodeCharge', '-ElectronCharge*kahan3(Holes, -Electrons, NetDoping)',
           ('Potential', 'Electrons'))
    _edge(SI, 'ElectricField', '(Potential@n0-Potential@n1)*EdgeInverseLength', ('Potential',))
    _edge(SI, 'DField', 'Permittivity*ElectricField', ('Potential',))
    ds.equation(device=DEV, region=SI, name='PotentialEquation', variable_name='Potential',
                node_model='PotentialNodeCharge', edge_model='DField', variable_update=update_psi)

    # ---- electron continuity, SG flux exactly as devsim/python_packages/simple_dd.py
    _edge(SI, 'vdiff', '(Potential@n0 - Potential@n1)/V_t', ())
    ds.edge_model(device=DEV, region=SI, name='vdiff:Potential@n0', equation='V_t^(-1)')
    ds.edge_model(device=DEV, region=SI, name='vdiff:Potential@n1', equation='-vdiff:Potential@n0')
    ds.edge_model(device=DEV, region=SI, name='Bern01', equation='B(vdiff)')
    ds.edge_model(device=DEV, region=SI, name='Bern01:Potential@n0',
                  equation='dBdx(vdiff) * vdiff:Potential@n0')
    ds.edge_model(device=DEV, region=SI, name='Bern01:Potential@n1', equation='-Bern01:Potential@n0')
    jn = ('ElectronCharge*mu_n*EdgeInverseLength*V_t*'
          'kahan3(Electrons@n1*Bern01,  Electrons@n1*vdiff,  -Electrons@n0*Bern01)')
    _edge(SI, 'ElectronCurrent', jn, ('Electrons', 'Potential'))
    _model(SI, 'NCharge', '-ElectronCharge*Electrons', ('Electrons',))
    ds.node_model(device=DEV, region=SI, name='ElectronGeneration', equation='0')
    ds.equation(device=DEV, region=SI, name='ElectronContinuityEquation', variable_name='Electrons',
                time_node_model='NCharge', edge_model='ElectronCurrent',
                node_model='ElectronGeneration', variable_update=update_n)

    # ---- interface: continuous potential
    ds.interface_model(device=DEV, interface='ox_si', name='continuousPotential',
                       equation='Potential@r0-Potential@r1')
    ds.interface_model(device=DEV, interface='ox_si', name='continuousPotential:Potential@r0',
                       equation='1')
    ds.interface_model(device=DEV, interface='ox_si', name='continuousPotential:Potential@r1',
                       equation='-1')
    ds.interface_equation(device=DEV, interface='ox_si', name='PotentialEquation',
                          interface_model='continuousPotential', type='continuous')

    # ---- contact charge models (displacement current at circuit contacts)
    for r in (SI, OX):
        ds.edge_model(device=DEV, region=r, name='contactcharge_edge', equation='DField')
        for n in ('n0', 'n1'):
            ds.edge_model(device=DEV, region=r, name=f'contactcharge_edge:Potential@{n}',
                          equation=f'DField:Potential@{n}')

    # gate (oxide)
    ds.contact_node_model(device=DEV, contact='gate', name='gate_psi',
                          equation='Potential - (gate_bias + psi_b0)')
    ds.contact_node_model(device=DEV, contact='gate', name='gate_psi:Potential', equation='1')
    ds.contact_equation(device=DEV, contact='gate', name='PotentialEquation', node_model='gate_psi',
                        edge_charge_model='contactcharge_edge')
    # body: psi Dirichlet only; electron flux natural (zero)
    ds.contact_node_model(device=DEV, contact='body', name='body_psi', equation='Potential - psi_b0')
    ds.contact_node_model(device=DEV, contact='body', name='body_psi:Potential', equation='1')
    ds.contact_equation(device=DEV, contact='body', name='PotentialEquation', node_model='body_psi',
                        edge_charge_model='contactcharge_edge')

    # drain circuit
    ds.circuit_element(name='V1', n1='drain_bias', n2='0', value=0.0, acreal=1.0, acimag=0.0)

    # electron reservoirs (phi_n Dirichlet on whole lateral edge)
    for c, b in (('source', 'source_bias'), ('drain', 'drain_bias')):
        eq = f'Electrons - n_i*exp((Potential - {b})/V_t)'
        ds.contact_node_model(device=DEV, contact=c, name=f'{c}_n', equation=eq)
        ds.contact_node_model(device=DEV, contact=c, name=f'{c}_n:Electrons', equation='1')
        ds.contact_node_model(device=DEV, contact=c, name=f'{c}_n:Potential',
                              equation=f'-n_i*exp((Potential - {b})/V_t)/V_t')
        kw = dict(device=DEV, contact=c, name='ElectronContinuityEquation', node_model=f'{c}_n',
                  edge_current_model='ElectronCurrent')
        if c == 'drain':
            ds.contact_node_model(device=DEV, contact=c, name=f'{c}_n:{b}',
                                  equation=f'n_i*exp((Potential - {b})/V_t)/V_t')
            kw['circuit_node'] = b
        ds.contact_equation(**kw)


def build_mesh_explicit(xs, y_ox, y_si, corner='tie'):
    """Explicit tensor-product mesh via create_gmsh_mesh (node coordinates exact).

    xs: lateral node coordinates (cm, ascending, xs[0]=0, xs[-1]=L); y_ox: oxide lines
    ascending ending at 0; y_si: Si lines from 0 ascending. Each rectangle is split
    into two right triangles (diagonal EdgeCouple = 0). Contacts / interface are
    defined by explicit boundary edges, so a node may belong to a contact AND the
    interface (DEVSIM Device.cc allows this with an INFO warning).

    corner = 'tie'     : one-node oxide contacts oxc_s/oxc_d ON the triple points
                         (revised construction; discrete continuity, see report)
             'exclude' : Si lateral contacts start below the interface node; the Si
                         corner stays an interface node (endpoint separation)
             'none'    : no corner treatment (control: reproduces the original defect)
    """
    import numpy as np
    xs = np.asarray(xs, float); y_ox = np.asarray(y_ox, float); y_si = np.asarray(y_si, float)
    assert abs(y_ox[-1]) == 0.0 and y_si[0] == 0.0
    ys = np.concatenate([y_ox, y_si[1:]]); kif = len(y_ox) - 1          # interface row
    nx, ny = len(xs), len(ys)
    nid = lambda i, j: j * nx + i
    coords = []
    for j in range(ny):
        for i in range(nx):
            coords += [float(xs[i]), float(ys[j]), 0.0]
    names = [OX, SI, 'ox_si', 'gate', 'source', 'drain', 'body'] + (['oxc_s', 'oxc_d'] if corner == 'tie' else [])
    P = {n: k for k, n in enumerate(names)}
    el = []
    for j in range(ny - 1):
        reg = P[OX] if j < kif else P[SI]
        for i in range(nx - 1):
            a, b, c, d = nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)
            el += [2, reg, a, b, c, 2, reg, a, c, d]
    for i in range(nx - 1):
        el += [1, P['ox_si'], nid(i, kif), nid(i + 1, kif)]
        el += [1, P['gate'], nid(i, 0), nid(i + 1, 0)]
    j0 = kif + 1 if corner == 'exclude' else kif
    for j in range(j0, ny - 1):
        el += [1, P['source'], nid(0, j), nid(0, j + 1)]
        el += [1, P['drain'], nid(nx - 1, j), nid(nx - 1, j + 1)]
    for i in range(1, nx - 2):                       # body: corners excluded
        el += [1, P['body'], nid(i, ny - 1), nid(i + 1, ny - 1)]
    if corner == 'tie':
        el += [1, P['oxc_s'], nid(0, kif), nid(0, kif)]          # one-node contacts
        el += [1, P['oxc_d'], nid(nx - 1, kif), nid(nx - 1, kif)]
    m = 'mx'
    ds.create_gmsh_mesh(mesh=m, coordinates=coords, elements=el, physical_names=names)
    ds.add_gmsh_region(gmsh_name=OX, mesh=m, region=OX, material='Ox')
    ds.add_gmsh_region(gmsh_name=SI, mesh=m, region=SI, material='Si')
    ds.add_gmsh_interface(gmsh_name='ox_si', mesh=m, name='ox_si', region0=OX, region1=SI)
    contacts = [('gate', OX), ('source', SI), ('drain', SI), ('body', SI)]
    if corner == 'tie':
        contacts += [('oxc_s', OX), ('oxc_d', OX)]
    for c, r in contacts:
        ds.add_gmsh_contact(gmsh_name=c, mesh=m, name=c, region=r, material='metal')
    ds.finalize_mesh(mesh=m)
    ds.create_device(mesh=m, device=DEV)
    CORNER['on'] = corner == 'tie'
    CORNER['mode'] = 'tie' if corner == 'tie' else 'none'
    return dict(kind='explicit', corner=corner, L_cm=float(xs[-1]), nx=nx, n_ox_lines=len(y_ox),
                n_si_lines=len(y_si), xmax_cm=float(y_si[-1]), tox_cm=float(-y_ox[0]),
                h_si_first_cm=float(y_si[1] - y_si[0]))


def redeclare_updates(update_psi='default', update_n='default'):
    """Re-declare the bulk equations with a different variable_update (same models).
    'default' = undamped Newton update x <- x + dx (DEVSIM Equation.cc
    DefaultSolutionUpdate); 'log_damp'/'positive' damp the update."""
    ds.equation(device=DEV, region=OX, name='PotentialEquation', variable_name='Potential',
                edge_model='DField', variable_update=update_psi)
    ds.equation(device=DEV, region=SI, name='PotentialEquation', variable_name='Potential',
                node_model='PotentialNodeCharge', edge_model='DField', variable_update=update_psi)
    ds.equation(device=DEV, region=SI, name='ElectronContinuityEquation', variable_name='Electrons',
                time_node_model='NCharge', edge_model='ElectronCurrent',
                node_model='ElectronGeneration', variable_update=update_n)


def set_lateral_dirichlet(enable: bool):
    """Switch source/drain psi between R-I(a) (Dirichlet = PsiLat + s*bias) and
    R-I(b) (no Potential contact equation = homogeneous Neumann)."""
    for c, b in (('source', 'source_bias'), ('drain', 'drain_bias')):
        if enable:
            ds.contact_node_model(device=DEV, contact=c, name=f'{c}_psi',
                                  equation=f'Potential - (PsiLat + s_couple*{b})')
            ds.contact_node_model(device=DEV, contact=c, name=f'{c}_psi:Potential', equation='1')
            kw = dict(device=DEV, contact=c, name='PotentialEquation', node_model=f'{c}_psi',
                      edge_charge_model='contactcharge_edge')
            if c == 'drain':
                ds.contact_node_model(device=DEV, contact=c, name=f'{c}_psi:{b}',
                                      equation='-s_couple')
                kw['circuit_node'] = b
            ds.contact_equation(**kw)
        else:
            try:
                ds.delete_contact_equation(device=DEV, contact=c, name='PotentialEquation')
            except ds.error:
                pass
    if CORNER['on']:
        for c in ('oxc_s', 'oxc_d'):
            if enable and CORNER['mode'] == 'legacy_pin':
                # C-1 attempt-3 construction (built-in mesher): pins the oxide node ADJACENT
                # to the triple point; kept unchanged as regression evidence.
                ds.contact_node_model(device=DEV, contact=c, name=f'{c}_psi', equation='Potential - PsiLat')
                ds.contact_node_model(device=DEV, contact=c, name=f'{c}_psi:Potential', equation='1')
                ds.contact_equation(device=DEV, contact=c, name='PotentialEquation', node_model=f'{c}_psi')
            elif enable and CORNER['mode'] == 'tie':
                # Revised construction (explicit mesh): the one-node oxide contact sits ON the
                # triple point. Its value equals the Si corner Dirichlet value, i.e. it is the
                # discrete potential-continuity condition at a node whose Si twin is Dirichlet.
                # Its displacement flux is assigned to the same terminal (drain via circuit).
                b = 'drain_bias' if c == 'oxc_d' else 'source_bias'
                ds.contact_node_model(device=DEV, contact=c, name=f'{c}_psi',
                                      equation=f'Potential - (PsiLat + s_couple*{b})')
                ds.contact_node_model(device=DEV, contact=c, name=f'{c}_psi:Potential', equation='1')
                kw = dict(device=DEV, contact=c, name='PotentialEquation', node_model=f'{c}_psi',
                          edge_charge_model='contactcharge_edge')
                if c == 'oxc_d':
                    ds.contact_node_model(device=DEV, contact=c, name=f'{c}_psi:{b}', equation='-s_couple')
                    kw['circuit_node'] = b
                ds.contact_equation(**kw)
            else:
                try:
                    ds.delete_contact_equation(device=DEV, contact=c, name='PotentialEquation')
                except ds.error:
                    pass


def init_electrons_equilibrium():
    psi = ds.get_node_model_values(device=DEV, region=SI, name='Potential')
    ds.set_node_values(device=DEV, region=SI, name='Electrons',
                       values=[NI * math.exp(p / VT) for p in psi])


def solve(abs_err=1e10, rel_err=1e-12, it=40, **kw):
    return ds.solve(type='dc', absolute_error=abs_err, relative_error=rel_err,
                    maximum_iterations=it, solver_type='direct', info=True, **kw)


class NotConverged(RuntimeError):
    pass


def _state():
    s = {}
    for r, var in ((SI, 'Potential'), (SI, 'Electrons'), (OX, 'Potential')):
        s[(r, var)] = list(ds.get_node_model_values(device=DEV, region=r, name=var))
    for node in ds.get_circuit_node_list():
        s[('circuit', node)] = ds.get_circuit_node_value(solution='dcop', node=node)
    return s


def _restore(s):
    for (r, var), v in s.items():
        if r == 'circuit':
            ds.set_circuit_node_value(solution='dcop', node=var, value=v)
        else:
            ds.set_node_values(device=DEV, region=r, name=var, values=v)


TAU = 1e-10      # R05 D-3: ||dpsi||_inf <= 1e-10 V and ||dphi_n||_inf <= 1e-10 V


def newton(tau=TAU, max_it=40, extra=0):
    """Newton iterations one at a time with the R05 D-3 stopping rule.

    DEVSIM's own test uses rel = |dx|/(|x| + 1e-10) (Equation.cc, minError not
    exposed to Python), which has no floor where the intrinsic-referenced psi
    passes through 0, and its device-level abs error mixes V and cm^-3. Here:
      dpsi  = max over si/ox PotentialEquation of max|d psi|            [V]
      dphin = V_T * (si ElectronContinuityEquation rel) ~ max|d n|/n * V_T  [V]
    Converged when the last update has dpsi <= tau and dphin <= tau.
    `extra` further iterations are then taken and their updates recorded.
    """
    hist = []
    for k in range(max_it):
        info = solve(abs_err=1e300, rel_err=1e300, it=1)
        e = equation_errors(info)
        dpsi = max(e['si:PotentialEquation']['abs'], e['ox:PotentialEquation']['abs'])
        dphin = VT * e['si:ElectronContinuityEquation']['rel']
        hist.append(dict(dpsi=dpsi, dphin=dphin, dcirc=e['circuit_abs']))
        if not (math.isfinite(dpsi) and math.isfinite(dphin)) or dpsi > 50:
            break
        if dpsi <= tau and dphin <= tau:
            for _ in range(extra):
                info = solve(abs_err=1e300, rel_err=1e300, it=1)
                e = equation_errors(info)
                hist.append(dict(dpsi=max(e['si:PotentialEquation']['abs'], e['ox:PotentialEquation']['abs']),
                                 dphin=VT * e['si:ElectronContinuityEquation']['rel'],
                                 dcirc=e['circuit_abs'], extra=True))
            return dict(converged=True, iterations=k + 1, history=hist, final=hist[-1])
    return dict(converged=False, iterations=len(hist), history=hist, final=hist[-1] if hist else None)


def ramp_gate(vg_target, step=.05, min_step=1e-3, tau=TAU):
    """Gate continuation with step halving. Every attempt (incl. failures) is logged;
    a failed attempt restores the last converged state before retrying."""
    v = ds.get_parameter(device=DEV, name='gate_bias')
    log = []
    while abs(vg_target - v) > 1e-15:
        h = step
        while True:
            vn = v + max(min(vg_target - v, h), -h)
            saved = _state()
            ds.set_parameter(device=DEV, name='gate_bias', value=vn)
            r = newton(tau)
            log.append(dict(vg=vn, step=h, converged=r['converged'], iterations=r['iterations'],
                            final=r['final']))
            if r['converged']:
                v = vn
                break
            _restore(saved)
            ds.set_parameter(device=DEV, name='gate_bias', value=v)
            h /= 2
            if h < min_step:
                raise NotConverged(f'gate ramp failed at {vn} V; log={log}')
    return log


def converge(tau=TAU, extra=2, max_it=40):
    r = newton(tau, max_it=max_it, extra=extra)
    if not r['converged']:
        raise NotConverged(f'newton failed: {r["history"][-3:]}')
    return r


_STARTED = {'done': False}


def equilibrium(vg=1.25, tau=TAU):
    """Build the V_D = 0 R-I(b) state at V_G, copy its psi into PsiLat (D2),
    then enable R-I(a) lateral Dirichlet data. Returns run log.
    May be called repeatedly (V1 must be 0); lateral data are rebuilt each time."""
    log = {}
    if not _STARTED['done']:
        log['vg0'] = converge(tau, extra=0)                # V_G = 0 start
        _STARTED['done'] = True
    set_lateral_dirichlet(False)
    log['neumann_ramp'] = ramp_gate(vg, tau=tau)
    log['neumann_final'] = converge(tau)
    psi = ds.get_node_model_values(device=DEV, region=SI, name='Potential')
    ds.set_node_values(device=DEV, region=SI, name='PsiLat', values=list(psi))
    set_lateral_dirichlet(True)
    log['ria_final'] = converge(tau)
    return log


def equilibrium_1col(vg, tau=TAU):
    """R-I(a) V_D = 0 state from independent 1COL data (onecol.py) on identical y-lines.
    Sets PsiLat (Si and oxide), seeds Potential/Electrons with the laterally uniform
    1COL profile, enables lateral Dirichlet data (incl. oxide corners if built) and
    runs Newton. If D1-D4 hold the seed is the discrete solution and the first
    update is at round-off level."""
    import numpy as np
    from phase3.c1 import onecol
    xs, ys = map(np.array, coords(SI)); xo, yo = map(np.array, coords(OX))
    xc = xs[np.argmin(np.abs(xs - xs.max() / 2))]
    y_si = np.sort(ys[np.abs(xs - xc) < 1e-12]); y_ox = np.sort(yo[np.abs(xo - xc) < 1e-12])
    sol = onecol.solve_1col(y_ox, y_si, vg)
    lut = {round(float(yy), 15): float(pp) for yy, pp in zip(sol['y'], sol['psi'])}
    psi_si = [lut[round(float(yy), 15)] for yy in ys]
    psi_ox = [lut[round(float(yy), 15)] for yy in yo]
    ds.set_parameter(device=DEV, name='gate_bias', value=vg)
    for r, v in ((SI, psi_si), (OX, psi_ox)):
        ds.set_node_values(device=DEV, region=r, name='PsiLat', values=v)
        ds.set_node_values(device=DEV, region=r, name='Potential', values=v)
    ds.set_node_values(device=DEV, region=SI, name='Electrons', values=[NI * math.exp(p / VT) for p in psi_si])
    set_lateral_dirichlet(True)
    r = newton(tau, extra=2)
    if not r['converged']:
        raise NotConverged(f'equilibrium_1col failed: {r["history"][-3:]}')
    return dict(onecol=dict(converged=sol['converged'], iterations=sol['iterations'],
                            last_update_V=sol['last_update_V'], n_si_lines=len(y_si),
                            n_ox_lines=len(y_ox), psi_s_V=float(sol['psi'][sol['k']])),
                newton=r)


def equation_errors(info):
    it = info['iterations'][-1]
    out = dict(circuit_abs=it['circuit']['absolute_error'])
    for reg in it['devices'][0]['regions']:
        for e in reg['equations']:
            out[f"{reg['name']}:{e['name']}"] = dict(abs=e['absolute_error'],
                                                     rel=e['relative_error'])
    return out


def coords(region=SI):
    x = ds.get_node_model_values(device=DEV, region=region, name='x')
    y = ds.get_node_model_values(device=DEV, region=region, name='y')
    return list(x), list(y)


def contact_nodes(contact, region=SI):
    nodes = set()
    for e in ds.get_element_node_list(device=DEV, region=region, contact=contact):
        nodes.update(e)
    return sorted(nodes)
