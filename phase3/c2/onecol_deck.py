"""C-2 one-column (1COL) DEVSIM deck, two formulations, D-1 node coordinates.

Physics: R05 v2.1 section 0.2 PM checklist, family F0 (identical constants to the C-1
deck: phase3/c1/deck.py). Units cm. Depth coordinate = DEVSIM y (0 = Si/SiO2
interface, y > 0 silicon, y < 0 oxide); lateral = DEVSIM x.

Formulation P (reduced, primary for T0/T2):
  unknown psi only; n = n_i exp((psi - phi_n)/V_T) with a spatially UNIFORM phi_n
  (parameter), p = n_i exp(-psi/V_T) (frozen holes). Column one cell wide (1 nm).
  Lateral edges: no contact => homogeneous Neumann for psi.
Formulation F (full, equivalence check):
  unknowns psi and n; electron continuity with the SG flux of the C-1 deck
  (identical to devsim/python_packages/simple_dd.py). Column TWO cells wide so the
  middle column is a genuine interior where continuity is solved. Lateral Si edges:
  electron contacts n = n_i exp((psi - phi_n)/V_T) (phi_n Dirichlet, same value on
  both sides), NO potential contact (Neumann psi). Because no lateral contact solves
  PotentialEquation, DEVSIM keeps every Si/SiO2 interface pair active (corner pairs
  included; InterfaceEquation.cc RemoveContactNodesFromList).
Common: gate contact on the oxide top psi = V_G + psi_b0; body contact on the Si
  bottom psi = psi_b0 (all bottom nodes); electrons have no bottom contact (zero
  normal flux, natural); interface: continuous potential.

Explicit mesh via create_gmsh_mesh (exact D-1 coordinates; no inserted nodes).
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np

import devsim as ds

from phase3 import meshmap as MM, runlog
from phase3.c1 import deck as C1          # shared constants only (PM checklist)

DEV, SI, OX = 'col', 'si', 'ox'
VT, NI, NA, Q, MU_N, PSI_B0, EPS0 = C1.VT, C1.NI, C1.NA, C1.Q, C1.MU_N, C1.PSI_B0, C1.EPS0_CM
WIDTH_CM = 1e-7                           # 1 nm (R05 C-2)


def build(level, form='P'):
    y_si = MM.f0_depth_map().nodes(level)
    y_ox = MM.oxide_map(20.0, 5).nodes(level)
    xs = np.array([0.0, WIDTH_CM]) if form == 'P' else np.array([0.0, WIDTH_CM / 2, WIDTH_CM])
    ys = np.concatenate([y_ox, y_si[1:]]); kif = len(y_ox) - 1
    nx, ny = len(xs), len(ys)
    nid = lambda i, j: j * nx + i
    coords = [v for j in range(ny) for i in range(nx) for v in (float(xs[i]), float(ys[j]), 0.0)]
    names = [OX, SI, 'ox_si', 'gate', 'body'] + (['left', 'right'] if form == 'F' else [])
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
        el += [1, P['body'], nid(i, ny - 1), nid(i + 1, ny - 1)]
    if form == 'F':
        for j in range(kif, ny - 1):
            el += [1, P['left'], nid(0, j), nid(0, j + 1)]
            el += [1, P['right'], nid(nx - 1, j), nid(nx - 1, j + 1)]
    ds.create_gmsh_mesh(mesh='m', coordinates=coords, elements=el, physical_names=names)
    ds.add_gmsh_region(gmsh_name=OX, mesh='m', region=OX, material='Ox')
    ds.add_gmsh_region(gmsh_name=SI, mesh='m', region=SI, material='Si')
    ds.add_gmsh_interface(gmsh_name='ox_si', mesh='m', name='ox_si', region0=OX, region1=SI)
    for c, r in [('gate', OX), ('body', SI)] + ([('left', SI), ('right', SI)] if form == 'F' else []):
        ds.add_gmsh_contact(gmsh_name=c, mesh='m', name=c, region=r, material='metal')
    ds.finalize_mesh(mesh='m'); ds.create_device(mesh='m', device=DEV)
    _physics(form)
    return dict(level=level, form=form, n_si_lines=len(y_si), n_ox_lines=len(y_ox), nx=nx,
                h_si_first_cm=float(y_si[1] - y_si[0]), h_si_max_cm=float(np.max(np.diff(y_si))),
                quality_depth=MM.f0_depth_map().quality(level), quality_oxide=MM.oxide_map().quality(level))


def _nm(r, name, eq, wrt=()):
    ds.node_model(device=DEV, region=r, name=name, equation=eq)
    for v in wrt:
        ds.node_model(device=DEV, region=r, name=f'{name}:{v}', equation=f'simplify(diff({eq},{v}))')


def _em(r, name, eq, wrt=()):
    ds.edge_model(device=DEV, region=r, name=name, equation=eq)
    for v in wrt:
        for n in ('n0', 'n1'):
            ds.edge_model(device=DEV, region=r, name=f'{name}:{v}@{n}', equation=f'simplify(diff({eq},{v}@{n}))')


def _physics(form):
    for r, e in ((SI, 11.7), (OX, 3.9)):
        ds.set_parameter(device=DEV, region=r, name='Permittivity', value=e * EPS0)
        ds.set_parameter(device=DEV, region=r, name='ElectronCharge', value=Q)
    for n, v in (('n_i', NI), ('V_t', VT), ('mu_n', MU_N), ('NetDoping', -NA)):
        ds.set_parameter(device=DEV, region=SI, name=n, value=v)
    for n, v in (('gate_bias', 0.0), ('psi_b0', PSI_B0), ('phin', 0.0)):
        ds.set_parameter(device=DEV, name=n, value=v)
    for r in (SI, OX):
        ds.node_solution(device=DEV, region=r, name='Potential')
        ds.edge_from_node_model(device=DEV, region=r, node_model='Potential')
        ds.set_node_value(device=DEV, region=r, name='Potential', value=PSI_B0)
        _em(r, 'ElectricField', '(Potential@n0-Potential@n1)*EdgeInverseLength', ('Potential',))
        _em(r, 'DField', 'Permittivity*ElectricField', ('Potential',))
    ds.equation(device=DEV, region=OX, name='PotentialEquation', variable_name='Potential',
                edge_model='DField', variable_update='log_damp')
    _nm(SI, 'Holes', 'n_i*exp(-Potential/V_t)', ('Potential',))
    if form == 'P':
        _nm(SI, 'Electrons', 'n_i*exp((Potential-phin)/V_t)', ('Potential',))
        _nm(SI, 'PotentialNodeCharge', '-ElectronCharge*kahan3(Holes, -Electrons, NetDoping)', ('Potential',))
    else:
        ds.node_solution(device=DEV, region=SI, name='Electrons')
        ds.set_node_value(device=DEV, region=SI, name='Electrons', value=NI * math.exp(PSI_B0 / VT))
        _nm(SI, 'PotentialNodeCharge', '-ElectronCharge*kahan3(Holes, -Electrons, NetDoping)',
            ('Potential', 'Electrons'))
    ds.equation(device=DEV, region=SI, name='PotentialEquation', variable_name='Potential',
                node_model='PotentialNodeCharge', edge_model='DField', variable_update='log_damp')
    # SG electron flux (evaluated in both formulations; solved only in F)
    ds.edge_from_node_model(device=DEV, region=SI, node_model='Electrons')
    _em(SI, 'vdiff', '(Potential@n0 - Potential@n1)/V_t')
    ds.edge_model(device=DEV, region=SI, name='Bern01', equation='B(vdiff)')
    jn = ('ElectronCharge*mu_n*EdgeInverseLength*V_t*'
          'kahan3(Electrons@n1*Bern01,  Electrons@n1*vdiff,  -Electrons@n0*Bern01)')
    # term scale for cancellation diagnostics: |q mu V_t/h| * max(n0 B, n1(B+v))
    ds.edge_model(device=DEV, region=SI, name='ElectronCurrentScale',
                  equation='ElectronCharge*mu_n*EdgeInverseLength*V_t*'
                           'max(abs(Electrons@n0*Bern01), abs(Electrons@n1*(Bern01+vdiff)))')
    if form == 'F':
        ds.edge_model(device=DEV, region=SI, name='vdiff:Potential@n0', equation='V_t^(-1)')
        ds.edge_model(device=DEV, region=SI, name='vdiff:Potential@n1', equation='-vdiff:Potential@n0')
        ds.edge_model(device=DEV, region=SI, name='Bern01:Potential@n0', equation='dBdx(vdiff) * vdiff:Potential@n0')
        ds.edge_model(device=DEV, region=SI, name='Bern01:Potential@n1', equation='-Bern01:Potential@n0')
        _em(SI, 'ElectronCurrent', jn, ('Electrons', 'Potential'))
        _nm(SI, 'NCharge', '-ElectronCharge*Electrons', ('Electrons',))
        ds.node_model(device=DEV, region=SI, name='ElectronGeneration', equation='0')
        ds.equation(device=DEV, region=SI, name='ElectronContinuityEquation', variable_name='Electrons',
                    time_node_model='NCharge', edge_model='ElectronCurrent',
                    node_model='ElectronGeneration', variable_update='positive')
        for c in ('left', 'right'):
            eq = 'Electrons - n_i*exp((Potential - phin)/V_t)'
            ds.contact_node_model(device=DEV, contact=c, name=f'{c}_n', equation=eq)
            ds.contact_node_model(device=DEV, contact=c, name=f'{c}_n:Electrons', equation='1')
            ds.contact_node_model(device=DEV, contact=c, name=f'{c}_n:Potential',
                                  equation='-n_i*exp((Potential - phin)/V_t)/V_t')
            ds.contact_equation(device=DEV, contact=c, name='ElectronContinuityEquation',
                                node_model=f'{c}_n', edge_current_model='ElectronCurrent')
    else:
        ds.edge_model(device=DEV, region=SI, name='ElectronCurrent', equation=jn)
    ds.interface_model(device=DEV, interface='ox_si', name='cont', equation='Potential@r0-Potential@r1')
    ds.interface_model(device=DEV, interface='ox_si', name='cont:Potential@r0', equation='1')
    ds.interface_model(device=DEV, interface='ox_si', name='cont:Potential@r1', equation='-1')
    ds.interface_equation(device=DEV, interface='ox_si', name='PotentialEquation', interface_model='cont',
                          type='continuous')
    for c, r, eq in (('gate', OX, 'Potential - (gate_bias + psi_b0)'), ('body', SI, 'Potential - psi_b0')):
        ds.contact_node_model(device=DEV, contact=c, name=f'{c}_psi', equation=eq)
        ds.contact_node_model(device=DEV, contact=c, name=f'{c}_psi:Potential', equation='1')
        ds.edge_model(device=DEV, region=r, name='cc_edge', equation='DField')
        for n in ('n0', 'n1'):
            ds.edge_model(device=DEV, region=r, name=f'cc_edge:Potential@{n}', equation=f'DField:Potential@{n}')
        ds.contact_equation(device=DEV, contact=c, name='PotentialEquation', node_model=f'{c}_psi',
                            edge_charge_model='cc_edge')


def errors(info):
    it = info['iterations'][-1]; e = {}
    for reg in it['devices'][0]['regions']:
        for q in reg['equations']:
            e[f"{reg['name']}:{q['name']}"] = (q['absolute_error'], q['relative_error'])
    return e


def newton(tau, form, max_it=60, extra=2):
    """R05 D-3 absolute stopping: ||dpsi||_inf <= tau and (F only) V_T*max|dn|/n <= tau."""
    hist = []
    for k in range(max_it):
        info = ds.solve(type='dc', absolute_error=1e300, relative_error=1e300, maximum_iterations=1,
                        solver_type='direct', info=True)
        e = errors(info)
        dpsi = max(e['si:PotentialEquation'][0], e['ox:PotentialEquation'][0])
        dphin = VT * e['si:ElectronContinuityEquation'][1] if form == 'F' else 0.0
        hist.append([dpsi, dphin])
        if not (math.isfinite(dpsi) and dpsi < 50):
            break
        if dpsi <= tau and dphin <= tau:
            for _ in range(extra):
                info = ds.solve(type='dc', absolute_error=1e300, relative_error=1e300, maximum_iterations=1,
                                solver_type='direct', info=True)
                e = errors(info)
                hist.append([max(e['si:PotentialEquation'][0], e['ox:PotentialEquation'][0]),
                             VT * e['si:ElectronContinuityEquation'][1] if form == 'F' else 0.0])
            return dict(converged=True, iterations=k + 1, final=hist[-1], first=hist[0])
    return dict(converged=False, iterations=len(hist), final=hist[-1] if hist else None, first=hist[0] if hist else None)


def measure(form):
    x = np.array(ds.get_node_model_values(device=DEV, region=SI, name='x'))
    y = np.array(ds.get_node_model_values(device=DEV, region=SI, name='y'))
    psi = np.array(ds.get_node_model_values(device=DEV, region=SI, name='Potential'))
    n = np.array(ds.get_node_model_values(device=DEV, region=SI, name='Electrons'))
    p = np.array(ds.get_node_model_values(device=DEV, region=SI, name='Holes'))
    phin = ds.get_parameter(device=DEV, name='phin'); vg = ds.get_parameter(device=DEV, name='gate_bias')
    col = np.where(x == 0.0)[0]; col = col[np.argsort(y[col])]
    yc = y[col]; h = np.diff(yc); w = np.zeros_like(yc); w[:-1] += h / 2; w[1:] += h / 2   # D4 dual widths
    ns = float(np.sum(n[col] * w))                                                         # cm^-2
    qs = float(Q * np.sum((p[col] - n[col] - NA) * w))                                      # C/cm^2
    psis = float(psi[col[0]])
    # oxide displacement at the interface (Gauss check): eps_ox * (psi_gate - psi_s)/t_ox on a
    # linear oxide (exact for the charge-free 1-D oxide) vs -Qs
    xo = np.array(ds.get_node_model_values(device=DEV, region=OX, name='x'))
    yo = np.array(ds.get_node_model_values(device=DEV, region=OX, name='y'))
    po = np.array(ds.get_node_model_values(device=DEV, region=OX, name='Potential'))
    oc = np.where(xo == 0.0)[0]; oc = oc[np.argsort(yo[oc])]
    d_ox = float(3.9 * EPS0 * (po[oc[0]] - po[oc[-1]]) / (yo[oc[-1]] - yo[oc[0]]))
    # quasi-Fermi relations and lateral uniformity
    rel_n = float(np.max(np.abs(n / (NI * np.exp((psi - phin) / VT)) - 1)))
    rel_p = float(np.max(np.abs(p / (NI * np.exp(-psi / VT)) - 1)))
    lat = 0.0
    for yy in np.unique(y):
        m = y == yy; lat = max(lat, float(np.ptp(psi[m])))
    # SG electron flux and continuity residual (sum of edge fluxes per node)
    F = np.array(ds.get_edge_model_values(device=DEV, region=SI, name='ElectronCurrent'))
    S = np.array(ds.get_edge_model_values(device=DEV, region=SI, name='ElectronCurrentScale'))
    ec = np.array(ds.get_edge_model_values(device=DEV, region=SI, name='EdgeCouple'))
    ds.edge_from_node_model(device=DEV, region=SI, node_model='node_index')
    n0 = np.array(ds.get_edge_model_values(device=DEV, region=SI, name='node_index@n0')).astype(int)
    n1 = np.array(ds.get_edge_model_values(device=DEV, region=SI, name='node_index@n1')).astype(int)
    div = np.zeros(len(x)); np.add.at(div, n0, F * ec); np.add.at(div, n1, -F * ec)
    scale = np.zeros(len(x)); np.add.at(scale, n0, S * ec); np.add.at(scale, n1, S * ec)
    act = scale > 0
    return dict(VG_V=vg, phin_V=phin, psi_s_V=psis, phi_s_bulkref_V=psis - PSI_B0, Ns_FV_cm2=ns, Qs_FV_C_cm2=qs,
                D_ox_C_cm2=d_ox, gauss_rel=abs(d_ox + qs) / max(abs(qs), 1e-300),
                gate_bc_err_V=float(po[oc[0]] - (vg + PSI_B0)),
                bottom_bc_err_V=float(psi[col[-1]] - PSI_B0),
                n_relation_max_rel=rel_n, p_relation_max_rel=rel_p, lateral_psi_spread_V=lat,
                sg_flux_max_abs_A_cm2=float(np.max(np.abs(F))),
                sg_flux_rel_to_term_scale=float(np.max(np.abs(F) / np.where(S > 0, S, 1))),
                continuity_residual_rel=float(np.max(np.abs(div[act]) / scale[act])),
                psi_profile=None)


def run(level, form, tau, outdir, vg_step=0.05):
    rec = runlog.header(dict(level=level, form=form, tau=tau, width_cm=WIDTH_CM, pm=C1.pm_checklist()))
    t0 = time.time()
    rec['mesh'] = build(level, form)
    rec['n_unknowns_hint'] = len(ds.get_node_model_values(device=DEV, region=SI, name='x'))
    states, fails = [], []
    r = newton(tau, form)
    if not r['converged']:
        fails.append(dict(vg=0.0, phin=0.0, newton=r))
    # V_G sweep at phi_n = 0 (R05 C-2)
    vgs = np.round(np.arange(0, 2.0 + 1e-9, vg_step), 10)
    for vg in vgs:
        ds.set_parameter(device=DEV, name='gate_bias', value=float(vg))
        r = newton(tau, form)
        m = measure(form); m['newton'] = r; m['status'] = 'ok' if r['converged'] else 'failed'
        states.append(m)
        if not r['converged']:
            fails.append(dict(vg=float(vg), phin=0.0, newton=r))
    # phi_n states at V_G in {0.8, 1.25, 2.0}: continuation in phi_n (0.05 V steps) from phi_n = 0
    snap = {}
    for vg in (0.8, 1.25, 2.0):
        ds.set_parameter(device=DEV, name='gate_bias', value=vg); ds.set_parameter(device=DEV, name='phin', value=0.0)
        newton(tau, form)
        for target in (0.05, 0.5):
            cur = ds.get_parameter(device=DEV, name='phin')
            for v in np.round(np.arange(cur + 0.05, target + 1e-9, 0.05), 10):
                ds.set_parameter(device=DEV, name='phin', value=float(v)); r = newton(tau, form)
                if not r['converged']:
                    fails.append(dict(vg=vg, phin=float(v), newton=r))
            m = measure(form); m['newton'] = r; m['status'] = 'ok' if r['converged'] else 'failed'
            states.append(m)
        ds.set_parameter(device=DEV, name='phin', value=0.0); newton(tau, form)
    rec.update(states=states, failures=fails, wall_s=time.time() - t0, end_utc=runlog.now())
    # psi profile at T0 points for the record
    prof = {}
    for vg in (0.0, 1.0, 2.0):
        ds.set_parameter(device=DEV, name='gate_bias', value=vg); newton(tau, form)
        x = np.array(ds.get_node_model_values(device=DEV, region=SI, name='x'))
        y = np.array(ds.get_node_model_values(device=DEV, region=SI, name='y'))
        psi = np.array(ds.get_node_model_values(device=DEV, region=SI, name='Potential'))
        c = np.where(x == 0.0)[0]; c = c[np.argsort(y[c])]
        prof[f'VG={vg}'] = dict(y_cm=y[c].tolist(), psi_V=psi[c].tolist())
    rec['psi_profiles'] = prof
    fn = os.path.join(outdir, f'{form}_L{level}_tau{tau:g}.json')
    runlog.write_json(fn, rec)
    return fn, len(fails)


def run_points(level, form, tau, points, outdir, tag):
    """Selected (V_G, phi_n) points only (used for U_iter re-runs at tau/100)."""
    rec = runlog.header(dict(level=level, form=form, tau=tau, points=points, tag=tag))
    rec['mesh'] = build(level, form); states = []
    newton(tau, form)
    cur = 0.0
    for vg, ph in points:
        # gate continuation at phi_n = 0, then phi_n continuation
        ds.set_parameter(device=DEV, name='phin', value=0.0); newton(tau, form)
        for v in np.round(np.arange(cur, vg + 1e-9, 0.05), 10) if vg >= cur else np.round(np.arange(cur, vg - 1e-9, -0.05), 10):
            ds.set_parameter(device=DEV, name='gate_bias', value=float(v)); newton(tau, form)
        ds.set_parameter(device=DEV, name='gate_bias', value=float(vg)); newton(tau, form); cur = vg
        for v in np.round(np.arange(0.05, ph + 1e-9, 0.05), 10):
            ds.set_parameter(device=DEV, name='phin', value=float(v)); newton(tau, form)
        ds.set_parameter(device=DEV, name='phin', value=float(ph)); r = newton(tau, form)
        m = measure(form); m['newton'] = r; states.append(m)
    rec.update(states=states, end_utc=runlog.now())
    fn = os.path.join(outdir, f'{form}_L{level}_tau{tau:g}_{tag}.json')
    runlog.write_json(fn, rec)
    return fn


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--level', type=int, required=True); ap.add_argument('--form', choices=['P', 'F'], required=True)
    ap.add_argument('--tau', type=float, default=1e-10); ap.add_argument('--out', required=True)
    ap.add_argument('--points', default=None, help='JSON list of [VG, phin] for a point run')
    ap.add_argument('--tag', default='pts')
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    if a.points:
        print('C2_POINTS_DONE', run_points(a.level, a.form, a.tau, json.loads(a.points), a.out, a.tag))
    else:
        fn, nf = run(a.level, a.form, a.tau, a.out)
        print('C2_RUN_DONE', fn, 'failures', nf)
        sys.exit(0 if nf == 0 else 3)
