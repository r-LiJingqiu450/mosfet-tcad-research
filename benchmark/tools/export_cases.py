#!/usr/bin/env python3
"""WP1 exporter: freeze existing Phase-3 results into benchmark case files.

Reads only existing JSON/CSV result files. Runs no device simulation and writes only under benchmark/.
Every exported number is copied from a source file or obtained from copied numbers by the closed-form
end-correction relation dL = (L/2)(1 - G_GCA/G) used by the producing scripts (phase3/c7/endprobe.py).

Determinism: no timestamps, sorted keys, fixed iteration order. Two consecutive runs must produce
byte-identical case files and MANIFEST.sha256 (checked by benchmark/tests).

Usage (from the repository root):  python benchmark/tools/export_cases.py
"""
from __future__ import annotations

import collections
import glob
import hashlib
import json
import math
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BENCH = os.path.join(ROOT, 'benchmark')
CASES = os.path.join(BENCH, 'cases')
BASE_COMMIT = 'c05139b4ae037ece34257ff0382e84a91b059773'
EXPORTER = 'benchmark/tools/export_cases.py'
SCHEMA_VERSION = '1.0.0'

# Physical constants exactly as used by the producing code (semiconductor_sim/common/constants.py, phase3/c1/deck.py)
Q = 1.602176634e-19
KB = 1.380649e-23
T_K = 300.0
VT = KB * T_K / Q
NI = 1.0e10
ND = 1.0e20
EPS_SI_REL, EPS_OX_REL = 11.7, 3.9
MU_N = 1350.0
BRACKET_RTOL = 1e-12          # same relative slack as phase3/c7/endprobe.py and phase3/c10/c10_analysis.py
COVERAGE_K = 2.0

A_W_VG = (0.488174, 0.621433, 0.7835)


# ----------------------------------------------------------------------------------------------- helpers
def rel(p):
    return os.path.relpath(p, ROOT).replace(os.sep, '/')


def sha256(path):
    h = hashlib.sha256()
    with open(os.path.join(ROOT, path), 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def git(*args):
    return subprocess.run(['git', '-C', ROOT, *args], check=True, capture_output=True, text=True).stdout.strip()


_BLOB_OK = {}


def check_blob(path):
    """The working-tree file must be byte-identical to the file at BASE_COMMIT (no silent edits)."""
    if path not in _BLOB_OK:
        at_base = git('rev-parse', f'{BASE_COMMIT}:{path}')
        here = git('hash-object', path)
        if at_base != here:
            raise SystemExit(f'source file differs from {BASE_COMMIT[:7]}: {path}')
        _BLOB_OK[path] = True
    return True


def load(path):
    check_blob(path)
    with open(os.path.join(ROOT, path)) as f:
        return json.load(f)


def src(path, role):
    check_blob(path)
    return dict(path=path, sha256=sha256(path), role=role)


def dl_nm(G, G_gca, L_nm):
    return 0.5 * L_nm * (1.0 - G_gca / G)


def g_tol_rel(tol_dl_nm, dl_ref_nm, L_nm):
    """Relative G tolerance equivalent to an absolute dL tolerance: d ln G = 2 d(dL) / (L - 2 dL)."""
    return 2.0 * tol_dl_nm / (L_nm - 2.0 * dl_ref_nm)


def fmt(x):
    return f'{x:g}'


def psi_b0_of(na):
    """Bulk equilibrium potential psi_b0 = -V_T ln(p0/n_i), p0 = (N_A + sqrt(N_A^2 + 4 n_i^2))/2 (phase3/c7/endprobe.set_bulk)."""
    p0 = 0.5 * (na + math.sqrt(na * na + 4 * NI * NI))
    return -VT * math.log(p0 / NI)


def physics(psi_b0):
    return dict(
        T_K=T_K, V_T_V=VT, n_i_cm3=NI, eps_si_rel=EPS_SI_REL, eps_ox_rel=EPS_OX_REL, mu_n_cm2_per_Vs=MU_N,
        statistics='Boltzmann (n = n_i exp((psi - phi_n)/V_T))',
        holes='not solved; equilibrium p = n_i exp(-psi/V_T) everywhere (hole quasi-Fermi level 0 V)',
        generation_recombination='none',
        mobility='constant mu_n (field-independent); zero-bias network edge conductance c_e = q mu_n EdgeCouple EIL n0 B(Delta)',
        psi_reference='intrinsic Fermi level (psi = 0 at n = p = n_i)',
        gate_convention='gate contact psi_gate = V_G + psi_b0 (V_FB = 0; no work-function difference, no oxide charge)',
        psi_b0_V=psi_b0,
        oxide='SiO2, no fixed or trapped charge, no carriers',
        units='2-D simulation; currents in A/cm and conductances in S/cm per unit device width')


def geometry(rec, deck):
    a, g = rec['args'], rec['geometry']
    out = dict(L_um=a['L'], t_ox_nm=a['tox'], h_w_um_requested=a['hw'], h_w_um_realized=g['h_w_cm'] * 1e4,
               si_depth_um=1.2, coordinates='x along the channel (source side first), y depth into Si (y = 0 at Si/SiO2 interface, oxide at -t_ox <= y <= 0)')
    if deck == 'W':
        out['si_x_range_um'] = [0.0, a['L']]
        out['note'] = 'Ward-truncated deck: silicon only between the metallurgical junction faces x = 0 and x = L'
    else:
        out['si_x_range_um'] = [-a['lw'], a['L'] + a['lw']]
        out['l_w_um'] = a['lw']
        out['gate_overlap_nm'] = a.get('ov', 0.0)
        out['gate_x_range_um'] = [-a.get('ov', 0.0) * 1e-3, a['L'] + a.get('ov', 0.0) * 1e-3]
    return out


def doping(rec, deck):
    a = rec['args']
    if deck == 'W':
        return dict(N_A_cm3=a['na'], profile='uniform acceptor N_A; no donors inside the domain (wells represented by boundary data)',
                    N_D_well_cm3_for_psi_bi=ND)
    s = a.get('sigma', 0.0)
    if s == 0.0:
        lat = 'abrupt: N_D(x, y) = 1e20 for (x < 0 or x > L) and y <= h_w, else 0'
    else:
        lat = ('lateral erfc: N_D(x, y) = 1e20 * [0.5 erfc((x - x_c)/(sqrt(2) sigma)) + 0.5 erfc(((L - x) - x_c)/(sqrt(2) sigma))] * [y <= h_w], '
               'x_c = -sqrt(2) sigma erfcinv(2 N_A / N_D), so that N_D = N_A (metallurgical junction) at x = 0 and x = L')
    return dict(N_A_cm3=a['na'], N_D_cm3=ND, sigma_nm=s, profile=lat, vertical='box: donors only for 0 <= y <= h_w (realized mesh line)')


def contacts(deck, s_couple=None):
    common = dict(
        gate='oxide top y = -t_ox: psi = V_G + psi_b0',
        body='Si bottom y = 1.2 um: psi = psi_b0; zero normal electron flux (no electron contact equation)',
        oxide_sides='homogeneous Neumann, except the two Si/SiO2 triple points at the source/drain faces, whose potential is tied to the adjacent silicon Dirichlet value')
    if deck == 'W':
        common.update(
            source='junction face x = 0, 0 <= y <= h_w: psi = psi_bi + psi_shift with psi_bi = V_T ln(N_D/n_i), N_D = 1e20; phi_n = 0',
            drain='junction face x = L, 0 <= y <= h_w: psi = psi_bi + psi_shift; phi_n = V_D (zero-bias network: unit perturbation, see observables)',
            faces_below_well='x = 0 and x = L for y > h_w: homogeneous Neumann for psi and zero normal electron flux')
    else:
        common.update(
            source='ohmic, far edge x = -l_w, 0 <= y <= h_w: psi = V_T asinh(N_net/(2 n_i)) (charge-neutral value); phi_n = 0',
            drain='ohmic, far edge x = L + l_w, 0 <= y <= h_w: psi = V_T asinh(N_net/(2 n_i)) + s_couple * V_D; phi_n = V_D',
            faces_below_well='far edges for y > h_w: homogeneous Neumann for psi and zero normal electron flux',
            oxide_top_outside_gate='homogeneous Neumann')
        if s_couple is not None:
            common['s_couple'] = s_couple
            common['s_couple_note'] = ('s_couple = 0 keeps the drain-contact psi at its equilibrium value when V_D is applied (non-standard); '
                                       's_couple = 1 shifts it by V_D (standard ohmic contact)')
    return common


def mesh(rec, script):
    g, m = rec['geometry'], rec['mesh']
    return dict(level=g['level'], n_lateral=g['n_lateral'], n_depth_si=g['n_depth_si'], n_ox=g['n_ox'], nx=m['nx'], ny=m['ny'],
                dx_min_nm=g['dx_min_nm'], dx_max_nm=g['dx_max_nm'], max_adjacent_ratio=g['max_adjacent_ratio'],
                generator=f'{script} (tensor-product node lines, each rectangle split into two triangles; refinement factor sqrt(2) per level)')


def cmd_e7(a):
    s = (f"python phase3/c7/endprobe.py --deck {a['deck']} --L {fmt(a['L'])} --level {a['level']} --vg {fmt(a['vg'])} "
         f"--tox {fmt(a['tox'])} --na {a['na']:.0e} --hw {fmt(a['hw'])}")
    if a['deck'] == 'R':
        s += f" --lw {fmt(a['lw'])} --ov {fmt(a['ov'])} --sigma {fmt(a['sigma'])}"
    else:
        s += f" --psibi-shift {fmt(a.get('psibi_shift', 0.0))}"
    if a.get('devsim'):
        s += ' --devsim'
    return s + ' --out <dir>   (reconstructed from the recorded args)'


def zb_state(path):
    """Zero-bias network quantities from a phase3/c7/endprobe.py record (constant mobility)."""
    r = load(path)
    c = r['result']['const_mu']
    a = r['args']
    L_nm = a['L'] * 1e3
    G_gca = r['result']['G_GCA_1COL']
    G_src = -c['src_current']
    G_drn = c['G_net']
    st = dict(rec=r, args=a, L_nm=L_nm, G_gca=G_gca, G_src=G_src, G_drn=G_drn, G_slice=c['G_slice'], G_rows=c['G_rows'],
              net_rel_residual=c['net_rel_residual'], bounds_ok_drain=bool(c['bounds_ok']))
    st['dL_src'] = dl_nm(G_src, G_gca, L_nm)
    st['dL_drn'] = dl_nm(G_drn, G_gca, L_nm)
    st['dL_slice'] = dl_nm(c['G_slice'], G_gca, L_nm)
    st['dL_rows'] = dl_nm(c['G_rows'], G_gca, L_nm)
    st['bounds_ok_src'] = bool(c['G_rows'] * (1 - BRACKET_RTOL) <= G_src <= c['G_slice'] * (1 + BRACKET_RTOL))
    return st


def zb_observables(st):
    return dict(
        G0_S_per_cm=dict(value=st['G_src'], path='source-side terminal current of the zero-bias network (phi_n = 0 source, 1 drain); '
                                                  'cancellation-free side, equal to the edge-energy / characteristic-potential form'),
        dL_nm=dict(value=st['dL_src'], definition='per-end correction dL = (L/2)(1 - G_GCA/G0)'),
        G_GCA_S_per_cm=dict(value=st['G_gca'], definition='long-channel (1-column) sheet conductance / L used as the dL reference'),
        bracket=dict(G_rows_S_per_cm=st['G_rows'], G_slice_S_per_cm=st['G_slice'], dL_rows_nm=st['dL_rows'], dL_slice_nm=st['dL_slice'],
                     rule=f'G_rows (1 - {BRACKET_RTOL:g}) <= G <= G_slice (1 + {BRACKET_RTOL:g})',
                     scope='Rayleigh bounds of THIS discrete network on THIS mesh; valid only for same-mesh submissions, not a bound on the continuum value'),
        network_relative_residual=st['net_rel_residual'])


def not_established(basis):
    return dict(status='not_established', components={}, coverage_factor=COVERAGE_K, tolerance='not_established', basis=basis)


def case_id(prefix, a, extra=''):
    s = f"{prefix}__{a['deck']}_L{fmt(a['L'])}_NA{a['na']:.0e}_tox{fmt(a['tox'])}_hw{fmt(a['hw'])}"
    if a['deck'] == 'W':
        s += f"_dpsi{fmt(a.get('psibi_shift', 0.0))}"
    else:
        s += f"_ov{fmt(a.get('ov', 0.0))}_sg{fmt(a.get('sigma', 0.0))}"
    s += f"_VG{fmt(a['vg'])}_lev{a['level']}{extra}"
    return s.replace('+', '')


def zb_case(tier, cid, st, sources, uncertainty, known_bad, cross, notes, related=()):
    r, a = st['rec'], st['args']
    assert abs(r['bulk']['psi_b0'] - psi_b0_of(a['na'])) < 1e-12, (cid, r['bulk']['psi_b0'])
    bias = dict(V_G_V=a['vg'], V_D_V=0.0, V_S_V=0.0)
    if a['deck'] == 'W':
        bias['psi_bi_shift_V'] = a.get('psibi_shift', 0.0)
    return dict(
        schema_version=SCHEMA_VERSION, id=cid, tier=tier, deck=a['deck'], family='zero_bias_network', related_cases=sorted(related),
        device=dict(geometry=geometry(r, a['deck']), doping=doping(r, a['deck']), contacts=contacts(a['deck']),
                    physics=physics(r['bulk']['psi_b0']), bias=bias),
        mesh=mesh(r, 'phase3/c7/endprobe.py'), observables=zb_observables(st), uncertainty=uncertainty,
        known_bad_outputs=known_bad, cross_paths=cross,
        provenance=dict(base_commit=BASE_COMMIT, sources=sources, producing_script='phase3/c7/endprobe.py', command=cmd_e7(a), exporter=EXPORTER),
        notes=notes)


def drain_entry(st, evidence):
    return dict(name='drain_side_flux_sum', same_mesh=True, G0_S_per_cm=st['G_drn'], dL_nm=st['dL_drn'],
                inside_bracket=st['bounds_ok_drain'], evidence=evidence)


ZB_NOTE = ('Zero-bias linear-response conductance: at equilibrium grad(phi_n) = 0, so the O(V_D) current is that of the medium '
           'sigma0 = q mu_n n0(x, y) and does not depend on the electrostatic response (phase3/reports/C6_Theory_ZeroBias_EndCorrection_2026-10-31.md, P4). '
           'Any drift-diffusion code can therefore be scored on G0 = lim I_D/V_D (V_D -> 0).')


# ----------------------------------------------------------------------------------------------- tiers
def export_bw():
    c10 = 'results/phase3/C10/analysis/c10_result.json'
    rows = load(c10)['all_states']
    assert len(rows) == 28, len(rows)
    cases = []
    for row in rows:
        st = zb_state(row['file'])
        assert abs(st['dL_src'] - row['dL']) <= 1e-9 * max(1.0, abs(row['dL'])), (row['file'], st['dL_src'], row['dL'])
        assert st['bounds_ok_drain'] == row['bounds_ok_drain'] and st['bounds_ok_src'] == row['bounds_ok_src'], row['file']
        a = st['args']
        cid = case_id('BW', a)
        if st['bounds_ok_drain']:
            kb, cross = [], [dict(drain_entry(st, 'inside the same-mesh bracket in the source record'), independent=False)]
        else:
            kb, cross = [drain_entry(st, 'C10 E10-1: drain-side current violates the same-mesh Rayleigh bracket (cancellation)')], []
        related = []
        if a['deck'] == 'W' and a['na'] == 1e18 and a['tox'] == 2.0 and a['L'] == 1.0 and a['hw'] == 0.5 \
                and a.get('psibi_shift', 0.0) == 0.0 and a['level'] == 0 and a['vg'] in A_W_VG:
            related.append(case_id('AW', a))
        cases.append(zb_case(
            'B-W', cid, st,
            [src(row['file'], 'endprobe state (geometry, mesh, network conductances)'),
             src(c10, 'C10 state index and E10-1 drain/source bracket flags')],
            not_established('bracketed tier: no independent path and no refinement pair exported for this state; scored by the same-mesh bracket only'),
            kb, cross, [ZB_NOTE, 'Ward-truncated deck W: standard TCAD cannot set up these boundary data; internal tier.'], related))
    return cases


def export_aw():
    eq = 'results/phase3/exploration/equivalence_audit/equivalence_summary.json'
    rows = {round(r['vg_V'], 6): r for r in load(eq)['rows']}
    fv = {0.488174: 'results/phase3/C11/indep/fd2d_VG0.488174.json',
          0.621433: 'results/phase3/exploration/equivalence_na1e18/fd2d_VG0.621433.json',
          0.7835: 'results/phase3/C11/indep/fd2d_VG0.783500.json'}
    lev1 = {0.7835: 'results/phase3/C10/runs/e7_W_L1_lev1_VG0.7835_tox2_NA1e+18_hw0.5_ov0_sg0_dpsi0.json'}
    c10_lev0 = {vg: f'results/phase3/C10/runs/e7_W_L1_lev0_VG{fmt(vg)}_tox2_NA1e+18_hw0.5_ov0_sg0_dpsi0.json' for vg in A_W_VG}
    cases = []
    for vg in A_W_VG:
        path = f'results/phase3/exploration/equivalence_energy_final/e7_W_L1_lev0_VG{fmt(vg)}_tox2_NA1e+18_hw0.5_ov0_sg0_dpsi0.json'
        st = zb_state(path)
        r, a, row = st['rec'], st['args'], rows[vg]
        L = st['L_nm']
        assert abs(row['G_source_S_per_cm'] / st['G_src'] - 1) < 1e-12, vg
        e = r['energy_identity']['constant_mobility']
        G_edge = e['edge_energy_G']
        dL_edge = dl_nm(G_edge, st['G_gca'], L)
        f = load(fv[vg])
        assert abs(f['args']['vg'] - vg) < 1e-9 and f['args']['na'] == a['na'] and f['args']['tox_nm'] == a['tox']
        dL_fv = f['dL_net_nm']
        assert abs(dL_fv - row['dL_independent_FV_nm']) < 1e-12
        c10 = zb_state(c10_lev0[vg])
        u_cross = abs(dL_fv - st['dL_src'])
        comps = dict(cross_path_dL_nm=u_cross)
        basis = ['discretization: |dL_independent_FV - dL_ref| (independent mesh, Newton and network assembly; C11 solver)']
        cross = [
            dict(name='edge_energy', same_mesh=True, independent=False, G0_S_per_cm=G_edge, dL_nm=dL_edge,
                 note='phi^T L phi of the same network; algebraically equal to the source current'),
            dict(name='independent_finite_volume', same_mesh=False, independent=True, dL_nm=dL_fv, mesh=dict(nx=f['nx'], ny=f['ny']),
                 note='phase3/c11/fd2d_indep.py: separate Poisson solve, mesh and network; half device by symmetry. Its own dL_slice/dL_rows '
                      'are trapezoid integrals, not exact bounds of its network, and are not exported as a bracket.',
                 fv_dL_slice_quadrature_nm=f['dL_slice_nm'], fv_dL_rows_quadrature_nm=f['dL_rows_nm']),
            dict(name='c10_original_run', same_mesh=True, independent=False, G0_S_per_cm=c10['G_src'], dL_nm=c10['dL_src'],
                 note='same arguments, earlier run (C10); identical mesh')]
        sources = [src(path, 'endprobe re-run with --devsim and edge-energy identity (reference state)'),
                   src(eq, 'equivalence audit rows (source, edge-energy, drain, G_E, G_M, independent FV)'),
                   src(fv[vg], 'independent finite-volume solver result'),
                   src(c10_lev0[vg], 'C10 original level-0 run (same arguments)')]
        if vg in lev1:
            s1 = zb_state(lev1[vg])
            u_lev = abs(s1['dL_src'] - st['dL_src'])
            comps['level_0_to_1_dL_nm'] = u_lev
            basis.append('discretization: |dL(level 1) - dL(level 0)| from C10')
            cross.append(dict(name='level_1_refinement', same_mesh=False, independent=False, G0_S_per_cm=s1['G_src'], dL_nm=s1['dL_src'],
                              mesh=dict(nx=s1['rec']['mesh']['nx'], ny=s1['rec']['mesh']['ny'])))
            sources.append(src(lev1[vg], 'C10 level-1 run (refinement pair)'))
        u_disc = max(comps.values())
        u_round = abs(dL_edge - st['dL_src'])
        u_solver = 0.5 * (L - 2 * st['dL_src']) * st['net_rel_residual']
        u_tot = u_disc + u_round + u_solver
        tol_dl = COVERAGE_K * u_tot
        unc = dict(status='established', coverage_factor=COVERAGE_K,
                   components=dict(discretization_dL_nm=u_disc, discretization_inputs=comps, roundoff_dL_nm=u_round, solver_dL_nm=u_solver,
                                   combined_dL_nm=u_tot),
                   tolerance=dict(dL_nm=tol_dl, G0_rel=g_tol_rel(tol_dl, st['dL_src'], L)),
                   basis='; '.join(basis) + '; round-off: |dL(edge energy) - dL(source)|; solver: network relative residual mapped to dL; '
                         'combined linearly (conservative), tolerance = k * combined with k = 2 (convention, not a statistical coverage). '
                         'The cross-path spread is itself part of the estimate, so passing it is a consistency check, not independent validation.')
        kb = [dict(name='coupled_small_signal_G_E', same_mesh=True, G0_S_per_cm=row['G_E_equilibrated_S_per_cm'],
                   rel_vs_reference=row['rel_G_E_minus_source'], kappa_G=row['kappa_G'], E_dwr=row['E_dwr'], E_round=row['E_round'],
                   evidence='DEVSIM full coupled Jacobian, row/column-equilibrated SuperLU solve; backward residual ~1e-16 yet forward error large (equivalence audit 2026-11-01)'),
              dict(name='coupled_small_signal_G_M', same_mesh=True, G0_S_per_cm=row['G_M_UMFPACK_S_per_cm'],
                   rel_vs_reference=row['rel_G_M_minus_source'],
                   evidence='same Newton linear system solved with UMFPACK (equivalence audit 2026-11-01)')]
        if st['bounds_ok_drain']:
            cross.append(dict(drain_entry(st, 'inside the same-mesh bracket'), independent=False))
        else:
            kb.append(drain_entry(st, 'drain-side current violates the same-mesh Rayleigh bracket (classical flux-cancellation failure)'))
        notes = [ZB_NOTE,
                 'Ward-truncated deck W: standard TCAD cannot set up these boundary data; internal reference tier.',
                 f"Realized well depth h_w = {r['geometry']['h_w_cm'] * 1e4:.6f} um (mesh line nearest to the requested 0.5 um)."]
        cases.append(zb_case('A-W', case_id('AW', a), st, sources, unc, kb, cross, notes, [case_id('BW', a)]))
    return cases


def export_c_rdeck():
    paths = sorted(rel(p) for p in glob.glob(os.path.join(ROOT, 'results', '**', 'e7_R_*.json'), recursive=True))
    groups = collections.defaultdict(list)
    for p in paths:
        groups[os.path.basename(p)].append(p)
    states, dup_info = {}, {}
    for name, ps in sorted(groups.items()):
        canon = [p for p in ps if '/C7/dev/' not in p] or ps
        assert len(canon) == 1, ps
        states[name] = canon[0]
        others = [p for p in ps if p != canon[0]]
        if others:
            dup_info[name] = others
    # refinement pairs: identical args except level
    def key(a):
        return (a['L'], a['vg'], a['tox'], a['na'], a['hw'], a['lw'], a['ov'], a['sigma'])
    loaded = {n: zb_state(p) for n, p in states.items()}
    by_key = collections.defaultdict(dict)
    for n, st in loaded.items():
        by_key[key(st['args'])][st['args']['level']] = n
    cases = []
    for name in sorted(loaded):
        st = loaded[name]
        a = st['args']
        sources = [src(states[name], 'endprobe state (deck R, resolved wells, ohmic contacts)')]
        cross, notes = [], [ZB_NOTE, 'Near/above-threshold resolved device (V_G >= 1.0 V); stress tier.']
        cross.append(dict(drain_entry(st, 'drain-side current of the same network'), independent=False))
        pair = by_key[key(a)]
        if len(pair) >= 2:
            assert sorted(pair) == [0, 1], pair
            other_name = pair[1 - a['level']]
            other = loaded[other_name]
            u_lev = abs(other['dL_src'] - st['dL_src'])
            tol = COVERAGE_K * u_lev
            unc = dict(status='established', coverage_factor=COVERAGE_K, components=dict(discretization_dL_nm=u_lev, discretization_inputs=dict(level_0_to_1_dL_nm=u_lev)),
                       tolerance=dict(dL_nm=tol, G0_rel=g_tol_rel(tol, st['dL_src'], st['L_nm'])),
                       basis='single refinement pair (level 0 <-> 1) only; no independent path; tolerance = 2 |dL(level 1) - dL(level 0)|')
            cross.append(dict(name=f"level_{other['args']['level']}_refinement", same_mesh=False, independent=False,
                              G0_S_per_cm=other['G_src'], dL_nm=other['dL_src']))
            sources.append(src(states[other_name], 'refinement partner'))
        else:
            unc = not_established('no refinement pair and no independent path for this state')
        for d in dup_info.get(name, []):
            dd = zb_state(d)
            same = all(dd[k] == st[k] for k in ('G_src', 'G_drn', 'G_slice', 'G_rows', 'G_gca'))
            entry = dict(name='duplicate_run', same_mesh=True, independent=False, path=d, identical_network_values=same)
            dv = dd['rec'].get('devsim')
            if dv:
                entry['devsim_coupled_jacobian'] = dict(G_E_S_per_cm=dv['G_E'], G_M_S_per_cm=dv['G_M'], kappa_G=dv['kappa_G'],
                                                        rel_G_E_vs_reference=dv['G_E'] / st['G_src'] - 1, rel_G_M_vs_reference=dv['G_M'] / st['G_src'] - 1)
            cross.append(entry)
            sources.append(src(d, 'duplicate file name in another directory'))
            notes.append(f'Duplicate file name also at {d}; network values identical: {same}.')
        cases.append(zb_case('C', case_id('C', a), st, sources, unc, [], cross, notes))
    return cases, dup_info


def export_c_welldepth():
    rows_csv = 'results/phase3/direction_reset/well_depth/well_depth_rows.csv'
    check_blob(rows_csv)
    cases = []
    for hw in ('0.1', '0.25', '0.5', '1'):
        p = f'results/phase3/direction_reset/well_depth/well_depth_L2_lev0_VG1.25_hw{hw}.json'
        r = load(p)
        d = r['decomposition']
        a = dict(r['args'], deck='R', ov=0.0)
        L_nm = a['L'] * 1e3
        G_src = -d['G_net_source_current']
        st = dict(G_src=G_src, G_drn=d['G_net'], G_slice=d['G_slice'], G_rows=d['G_rows'], G_gca=d['G_GCA'], net_rel_residual=d['net_rel_residual'],
                  dL_src=dl_nm(G_src, d['G_GCA'], L_nm), dL_drn=d['dL_net_nm'], dL_slice=d['dL_slice_nm'], dL_rows=d['dL_rows_nm'],
                  bounds_ok_drain=bool(d['bounds_ok']))
        assert abs(dl_nm(d['G_net'], d['G_GCA'], L_nm) - d['dL_net_nm']) < 1e-9
        sources = [src(p, 'well-depth probe state (deck R, abrupt wells)'), src(rows_csv, 'well-depth summary table')]
        cross = [dict(drain_entry(st, 'drain-side current of the same network'), independent=False)]
        dv_path = f'results/phase3/direction_reset/well_depth_devsim/well_depth_L2_lev0_VG1.25_hw{hw}.json'
        if os.path.exists(os.path.join(ROOT, dv_path)):
            dv = load(dv_path)['devsim']
            cross.append(dict(name='devsim_coupled_jacobian', same_mesh=True, independent=False, G_E_S_per_cm=dv['G_E'], G_M_S_per_cm=dv['G_M'],
                              kappa_G=dv['kappa_G'], E_dwr=dv['E_dwr'], E_round=dv['E_round'],
                              rel_G_E_vs_reference=dv['G_E'] / G_src - 1, rel_G_M_vs_reference=dv['G_M'] / G_src - 1))
            sources.append(src(dv_path, 'DEVSIM coupled-Jacobian cross-check run'))
        cid = f"C__Rwd_L{fmt(a['L'])}_NA{a['na']:.0e}_tox{fmt(a['tox'])}_hw{hw}_sg{fmt(a['sigma'])}_VG{fmt(a['vg'])}_lev{a['level']}".replace('+', '')
        obs = zb_observables(st)
        cases.append(dict(
            schema_version=SCHEMA_VERSION, id=cid, tier='C', deck='R', family='well_depth', related_cases=[],
            device=dict(geometry=geometry(r, 'R'), doping=doping(dict(args=a), 'R'), contacts=contacts('R'), physics=physics(psi_b0_of(a['na'])),
                        bias=dict(V_G_V=a['vg'], V_D_V=0.0, V_S_V=0.0)),
            mesh=mesh(r, 'phase3/direction_reset/well_depth_probe.py'), observables=obs,
            uncertainty=not_established('single mesh; bottom-boundary margin confounds the depth trend (Direction_Reset_Correction_and_Next_Investment_2026-11-01.md)'),
            known_bad_outputs=[], cross_paths=cross,
            provenance=dict(base_commit=BASE_COMMIT, sources=sources, producing_script='phase3/direction_reset/well_depth_probe.py',
                            command=(f"python phase3/direction_reset/well_depth_probe.py --L {fmt(a['L'])} --level {a['level']} --vg {fmt(a['vg'])} --hw {fmt(a['hw'])} "
                                     f"--na {a['na']:.0e} --tox {fmt(a['tox'])} --sigma {fmt(a['sigma'])} --lw {fmt(a['lw'])} --out <dir> --no-devsim   (reconstructed)"),
                            exporter=EXPORTER),
            notes=[ZB_NOTE, 'Bottom boundary fixed at y = 1.2 um; at h_w = 1 um the margin is 0.2 um (< W_dep), so depth and bottom-boundary effects are not separated.']))
    return cases


def export_c_finite_bias():
    groups = collections.defaultdict(list)
    for p in sorted(glob.glob(os.path.join(ROOT, 'results', 'phase3', 'next_question', 'physical_r2_*', 'r2_finite_*.json'))):
        rp = rel(p)
        r = load(rp)
        a = r['args']
        groups[(a['L'], a['vg'], a['s_couple'] if 's_couple' in a else 0.0, a['level'])].append((rp, r))
    cases, dup_report = [], []
    for k in sorted(groups):
        items = groups[k]
        by_vd = {}
        dup_max = 0.0
        for rp, r in items:
            for s in r['finite_states']:
                vd = round(s['vd_V'], 12)
                row = dict(V_D_V=s['vd_V'], I_drain_A_per_cm=s['I_drain_A_per_cm'], I_source_A_per_cm=s['I_source_A_per_cm'],
                           I_circuit_A_per_cm=s['I_circuit_A_per_cm'], G_eff_S_per_cm=s['G_eff_S_per_cm'], kcl_rel=s['kcl_rel'],
                           terminal_vs_circuit_rel=s['terminal_vs_circuit_rel'], cross_section_abs_spread_rel=s['cross_section_abs_spread_rel'],
                           source=rp)
                if vd in by_vd:
                    d = abs(row['I_drain_A_per_cm'] / by_vd[vd]['I_drain_A_per_cm'] - 1)
                    dup_max = max(dup_max, d)
                    by_vd[vd].setdefault('also_in', []).append(rp)
                else:
                    by_vd[vd] = row
        rp0, r0 = items[0]
        a = dict(r0['args'], deck='R')
        L, vg, sc, lev = k
        if any(len(v.get('also_in', [])) for v in by_vd.values()):
            dup_report.append(dict(group=f'L{fmt(L)}_VG{fmt(vg)}_s{fmt(sc)}', files=[i[0] for i in items], max_rel_I_drain_difference=dup_max))
        cid = f"C__R2fb_L{fmt(L)}_NA{a['na']:.0e}_tox{fmt(a['tox'])}_hw{fmt(a['hw'])}_sg{fmt(a['sigma'])}_VG{fmt(vg)}_s{fmt(sc)}_lev{lev}".replace('+', '')
        psi_b0 = psi_b0_of(a['na'])
        cases.append(dict(
            schema_version=SCHEMA_VERSION, id=cid, tier='C', deck='R', family='finite_bias_dd', related_cases=[],
            device=dict(geometry=geometry(r0, 'R'), doping=doping(dict(args=a), 'R'), contacts=contacts('R', s_couple=sc), physics=physics(psi_b0),
                        bias=dict(V_G_V=vg, V_S_V=0.0, V_D_V=[v['V_D_V'] for _, v in sorted(by_vd.items())])),
            mesh=mesh(r0, 'phase3/next_question/physical_r2_finite_bias.py'),
            observables=dict(finite_states=[v for _, v in sorted(by_vd.items())],
                             definition='full drift-diffusion at finite V_D with Newton continuation from zero; I_drain/I_source are terminal sums, I_circuit the circuit-node current'),
            uncertainty=not_established('single mesh (level 0) and single code path; stress tier'),
            known_bad_outputs=[], cross_paths=[],
            provenance=dict(base_commit=BASE_COMMIT, sources=[src(i[0], 'finite-bias run') for i in items],
                            producing_script='phase3/next_question/physical_r2_finite_bias.py',
                            command=(f"python phase3/next_question/physical_r2_finite_bias.py --L {fmt(L)} --level {lev} --vg {fmt(vg)} --vd <V_D list> "
                                     f"--tox {fmt(a['tox'])} --na {a['na']:.0e} --sigma {fmt(a['sigma'])} --ov {fmt(a['ov'])} --lw {fmt(a['lw'])} --hw {fmt(a['hw'])} "
                                     f"--s-couple {fmt(sc)} --out <dir>   (reconstructed)"),
                            exporter=EXPORTER),
            notes=['Finite-bias drift-diffusion stress case; not a reference. Rows from several files are merged by V_D; '
                   'duplicated V_D values are cross-checked (max relative I_drain difference recorded in the WP1 report).']))
    return cases, dup_report


# ----------------------------------------------------------------------------------------------- main
def dump(obj):
    return json.dumps(obj, sort_keys=True, indent=1, ensure_ascii=False, allow_nan=False) + '\n'


def main():
    import export_wp2   # WP2 tiers (schema 1.1.0); WP1 cases above are exported unchanged (schema 1.0.0)
    os.makedirs(CASES, exist_ok=True)
    for p in glob.glob(os.path.join(CASES, '*.json')):
        os.remove(p)
    aw = export_aw()
    bw = export_bw()
    cr, dup_r = export_c_rdeck()
    cw = export_c_welldepth()
    cf, dup_fb = export_c_finite_bias()
    wp2, mesh_files, wp2_info = export_wp2.export(sys.modules[__name__])
    allc = aw + bw + cr + cw + cf + wp2
    ids = [c['id'] for c in allc]
    assert len(ids) == len(set(ids)), [i for i, n in collections.Counter(ids).items() if n > 1]
    for c in allc:
        with open(os.path.join(CASES, c['id'] + '.json'), 'w', encoding='utf-8', newline='\n') as f:
            f.write(dump(c))
    meshes = os.path.join(BENCH, 'meshes')
    os.makedirs(meshes, exist_ok=True)
    for p in glob.glob(os.path.join(meshes, '*.json')):
        os.remove(p)
    for p, b in sorted(mesh_files.items()):
        with open(os.path.join(ROOT, p), 'wb') as f:
            f.write(b)
    index = dict(schema_version='1.1.0', schema_versions_used=sorted({c['schema_version'] for c in allc}), base_commit=BASE_COMMIT,
                 wp2=wp2_info,
                 counts=dict(sorted(collections.Counter(c['tier'] + '/' + c['family'] for c in allc).items())),
                 cases=[dict(id=c['id'], tier=c['tier'], deck=c['deck'], family=c['family'], uncertainty=c['uncertainty']['status'])
                        for c in sorted(allc, key=lambda c: c['id'])],
                 inconsistencies=dict(duplicate_rdeck_file_names={k: v for k, v in sorted(dup_r.items())}, finite_bias_duplicate_vd=dup_fb))
    with open(os.path.join(BENCH, 'INDEX.json'), 'w', encoding='utf-8', newline='\n') as f:
        f.write(dump(index))
    files = sorted(glob.glob(os.path.join(CASES, '*.json'))) + [os.path.join(BENCH, 'INDEX.json')] + \
        sorted(glob.glob(os.path.join(BENCH, 'schema', '*.json'))) + sorted(glob.glob(os.path.join(meshes, '*.json')))
    with open(os.path.join(BENCH, 'MANIFEST.sha256'), 'w', encoding='utf-8', newline='\n') as f:
        for p in files:
            f.write(f'{sha256(rel(p))}  {rel(p)}\n')
    print(json.dumps(index['counts']), f'total={len(allc)}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
