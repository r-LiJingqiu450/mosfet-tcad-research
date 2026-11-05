"""WP2 part of the exporter (schema 1.1.0): resolved-well tier A-R, held-out tier HO and B0 companions B-R.

Called by benchmark/tools/export_cases.py; not a stand-alone script. Reads only committed WP2 result files
(benchmark/wp2/results/..., byte-identical to WP2_COMMIT) and runs no simulation.

Definitions follow the frozen pre-registration benchmark/specs/WP3_Preregistration_AR_and_HeldOut_2026-11-04.md:
  * B1 (ohmic body) three-terminal observables G_DD, G_SD (scored) and G_DB (reported), section 3;
  * reference = P1 at level 2, tolerance = 2 (u_disc + u_round + u_solver), section 7, recomputed here and
    asserted equal to benchmark/wp2/results/evaluation.json (written by benchmark/wp2/evaluate_wp2.py);
  * tier A-R only for states that pass H1 (section 8); a state that fails H1 is exported as tier C;
  * same-mesh B1 scoring with tau_num = 1e-6 (section 7) against the P1 value on the declared mesh;
  * B0 companions are scored by the same-mesh bracket only (section 10, step 6) -> tier B-R, one case per mesh;
  * held-out states (tier HO): level 1 only, same-mesh scoring only.
Every known output path of the project's own codes (direct drain flux sum, DEVSIM contact currents, circuit-node
and weighting-function currents, coupled small-signal G_E / G_M, the independent FV solver) is attached to the
case it estimates and classified against that case's own acceptance rule: |rel| <= tau_num on the same mesh,
|rel| <= tolerance on another mesh (bracket for B-R). Classification is mechanical; nothing is hand-assigned.
"""
from __future__ import annotations

import json
import math
import os

WP2_COMMIT = '454db70e71660a896b23ca378d2f5c3413808204'     # commit that contains the WP2 production results, ledger, evaluation and meshes
SCHEMA_VERSION = '1.1.0'
RES = 'benchmark/wp2/results'
PROD = RES + '/production'
MESH_SRC = RES + '/meshes'
MESH_PUB = 'benchmark/meshes'
TAU_NUM = 1e-6
PREREG = 'benchmark/specs/WP3_Preregistration_AR_and_HeldOut_2026-11-04.md'
# WP4 (benchmark/specs/WP4_Preregistration_External_Paths_2026-11-04.md): when WP4_COMMIT is set, the A-R cases also carry
# the WP4 output paths (DEVSIM extended precision; DEVSIM built-in AC analysis), all on the level-1 mesh; those cases then
# cite WP4_COMMIT as their base commit (it contains every WP2 source unchanged).
WP4_COMMIT = None
WP4_PROD = 'benchmark/wp4/results/production'
PREREG_WP4 = 'benchmark/specs/WP4_Preregistration_External_Paths_2026-11-04.md'
TOL_FORMS = 1e-9
RICHARDSON = 'G_R = (100 G(+-0.1 mV) - G(+-1 mV)) / 99; removes the physical (V_D/V_T)^2/6 term of the symmetric difference (residual ~2e-10)'
GEOM_OF = {('dev1', 1.0): 'dev1_L1', ('dev2', 2.0): 'dev2_L2', ('dev1', 2.0): 'dev1_L2'}

ZB3_NOTE = ('Zero-bias linear response of the three-terminal electron network (ohmic body B1): at equilibrium grad(phi_n) = 0, so the '
            'O(V_D) electron current is that of sigma0 = q mu_n n0(x, y) (C6 theory, P4). Excite phi_n = 1 at the drain and 0 at source '
            'and body: G_DD = I_D, G_SD = -I_S, G_DB = -I_B, G_DD = G_SD + G_DB. Any drift-diffusion code can be scored on these '
            'limits (V_D -> 0); report the ELECTRON component of each terminal current (pre-registration section 2.5).')


class Ctx:
    """Committed-file access: every source must be byte-identical to its file at `commit` (default WP2_COMMIT)."""

    def __init__(self, E, commit=None):
        self.E, self.ok, self.commit = E, {}, commit or WP2_COMMIT
        if not (len(self.commit) == 40 and all(ch in '0123456789abcdef' for ch in self.commit)):
            raise SystemExit('export commit is not set to a commit hash')

    def check(self, path):
        if path not in self.ok:
            at = self.E.git('rev-parse', f'{self.commit}:{path}')
            here = self.E.git('hash-object', path)
            if at != here:
                raise SystemExit(f'source file differs from {self.commit[:7]}: {path}')
            self.ok[path] = True

    def load(self, path):
        self.check(path)
        with open(os.path.join(self.E.ROOT, path), encoding='utf-8') as f:
            return json.load(f)

    def src(self, path, role):
        self.check(path)
        return dict(path=path, sha256=self.E.sha256(path), role=role)


def close(a, b, rtol=1e-12):
    return abs(a - b) <= rtol * max(abs(a), abs(b), 1e-300)


def order(q0, q1, q2):
    d1, d2 = q1 - q0, q2 - q1
    if d1 == 0 or d2 == 0 or (d1 > 0) != (d2 > 0):
        return None
    return math.log(abs(d1) / abs(d2)) / math.log(math.sqrt(2.0))


def p1(sid, body, lev):
    return f'{PROD}/P1_{sid}_{body}_lev{lev}.json'


def p2(sid, body, lev):
    return f'{PROD}/P2_{sid}_{body}_lev{lev}.json'


def p4(sid, mesh):
    return f'{PROD}/P4_{sid}_{mesh}.json'


def mesh_id(rec):
    a = rec['args']
    dev = 'dev1' if a['na'] == 1e18 else 'dev2'
    return f"{GEOM_OF[(dev, float(a['L']))]}_lev{a['level']}"


def fmt(x):
    return f'{x:g}'


def state_tag(a):
    return f"R_L{fmt(a['L'])}_NA{a['na']:.0e}_tox{fmt(a['tox'])}_VG{fmt(a['vg'])}".replace('+', '')


# ------------------------------------------------------------------------------------------- device description
def geometry(rec):
    a, g = rec['args'], rec['geometry']
    return dict(L_um=a['L'], t_ox_nm=a['tox'], si_x_range_um=[-0.3, a['L'] + 0.3], si_depth_um=1.2, l_w_um=0.3,
                gate_overlap_nm=0.0, gate_x_range_um=[0.0, a['L']], h_w_um=0.5,
                contact_extent_um_this_mesh=g['contact_extent_cm'] * 1e4,
                contact_extent_note='source/drain contacts cover 0 <= y <= h_w; on a tensor mesh the realized extent is the depth line '
                                    'nearest h_w (this mesh: contact_extent_um_this_mesh). Its effect on G_DD/G_SD/G_DB is <= 6.5e-9 '
                                    'relative (benchmark/wp2/results/contact_extent_check.json, amendment A4).',
                coordinates='x along the channel (source side first), y depth into Si (y = 0 at the Si/SiO2 interface, oxide at -t_ox <= y <= 0)')


def doping(na):
    return dict(N_A_cm3=na, N_D_peak_cm3=1e20, sigma_nm=15.0, sigma_y_nm=15.0, h_w_um=0.5,
                profile=('N_A uniform plus donors N_D(x, y) = 1e20 f_x(x) f_y(y): f_x = 0.5 erfc((x - x_c)/(sqrt(2) sigma)) + '
                         '0.5 erfc(((L - x) - x_c)/(sqrt(2) sigma)), x_c = -sqrt(2) sigma erfcinv(2 N_A/N_D) (N_D = N_A at x = 0 and x = L); '
                         'f_y = 0.5 erfc((y - h_w)/(sqrt(2) sigma_y)); fully analytic, no doping step on a mesh line'))


def contacts(body):
    c = dict(
        gate='oxide top y = -t_ox over 0 <= x <= L: psi = V_G + psi_b0 (V_FB = 0)',
        source='ohmic, far edge x = -l_w, 0 <= y <= h_w: psi = V_T asinh(N_net/(2 n_i)); phi_n = 0',
        drain='ohmic, far edge x = L + l_w, 0 <= y <= h_w: psi = V_T asinh(N_net/(2 n_i)) + V_D; phi_n = V_D',
        faces_below_well='far edges for y > h_w: homogeneous Neumann for psi and zero normal electron flux',
        oxide_top_outside_gate='homogeneous Neumann',
        oxide_triple_points='the Si/SiO2 nodes at x = -l_w and x = L + l_w are tied to the adjacent contact psi')
    if body == 'B1':
        c['body'] = 'B1 (ohmic, primary): Si bottom y = 1.2 um, full width: psi = psi_b0, phi_n = 0 (n = n_i exp(psi_b0/V_T)); holes at equilibrium'
    else:
        c['body'] = ('B0 (electron-blocking companion, as all existing decks): Si bottom y = 1.2 um: psi = psi_b0 and ZERO normal electron '
                     'flux (no electron contact equation). In a code with ohmic-only contacts: zero electron surface recombination velocity.')
    return c


def physics(E, psi_b0, na):
    p = E.physics(psi_b0)
    wf = 1.0362 if na == 1e18 else 0.9172
    p.update(holes='not part of the observable: report the electron component of each terminal current, or hold holes at equilibrium '
                   'p = n_i exp(-psi/V_T) (pre-registration section 2.5); the body is ohmic for holes in B0 and B1',
             gate_work_function=f'phi_m = chi_Si + {wf} V (equals the p-substrate work function, so V_FB = 0); E_c - E_i = 0.56 eV')
    return p


def mesh_obj(rec, files):
    mid = mesh_id(rec)
    g, m = rec['geometry'], rec['mesh']
    pub = f'{MESH_PUB}/mesh_{mid}.json'
    return dict(level=g['level'], mesh_id=mid, nx=m['nx'], ny=m['ny'], n_lateral=g['n_lateral'], n_depth_si=g['n_depth_si'], n_ox=g['n_ox'],
                dx_min_nm=g['dx_min_nm'], dx_max_nm=g['dx_max_nm'], max_adjacent_ratio=g['max_adjacent_ratio'],
                node_lines_file=pub, node_lines_file_sha256=files[pub],
                generator=('endprobe deck-R maps (phase3.c6.rprobe.lateral_nodes, phase3.c7.endprobe.depth_map, oxide_map(t_ox, 5)); '
                           'refinement factor sqrt(2) per level; each rectangle split into two triangles (node_lines_file has the exact '
                           'node lines and rules); box integration, Scharfetter-Gummel edge currents (DEVSIM 2.11.2)'))


# ------------------------------------------------------------------------------------------- classification
def judge(entry, q, ref, same_mesh, tol_rel):
    """Mechanical classification of one output path: True = consistent with the reference under the case's own rule.
    Paths computed on a published DEVSIM mesh are compared with the reference-path value on THAT mesh (tau_num), so that
    discretization differences between levels are never counted as failures; the independent P4 solver (own mesh) is
    compared with the level-2 reference under the established tolerance."""
    v = entry.get(f'{q}_S_per_cm')
    if v is None or not v > 0:
        entry[f'rel_{q}_vs_reference'] = None if v is None else (v / ref - 1.0)
        entry['rule'] = 'sign: a zero-bias conductance is positive'
        return False
    r = v / ref - 1.0
    entry[f'rel_{q}_vs_reference'] = r
    if same_mesh:
        entry['rule'] = f'same mesh: |rel| <= tau_num = {TAU_NUM:g}'
        return abs(r) <= TAU_NUM
    if tol_rel is None:
        entry['rule'] = 'other mesh, no established tolerance: reported only'
        return None
    entry['rule'] = f'other mesh: |rel| <= tolerance {tol_rel:.3e}'
    return abs(r) <= tol_rel


def place(entry, verdicts, known_bad, cross):
    """Known bad if any judged observable is inconsistent; cross path if all judged ones are consistent."""
    vs = [v for v in verdicts if v is not None]
    if vs and not all(vs):
        rels = ', '.join(f"{k[4:-13]} rel = {v:+.3e}" if v is not None else f"{k[4:-13]} not positive"
                         for k, v in entry.items() if k.startswith('rel_G_') and k.endswith('_vs_reference'))
        entry['evidence'] = f"fails the case rule ({entry.get('rule', 'sign')}): {rels}"
        known_bad.append(entry)
    else:
        entry.setdefault('independent', False)
        cross.append(entry)


def p2_paths(d, refs, mid, suffix='', path_label='P2'):
    """Output paths of one P2 drift-diffusion run (DEVSIM contact currents, circuit node, weighting function), judged against the
    reference-path values `refs` on the same mesh `mid`."""
    lev = d['args']['level']
    s = d['symmetric']
    out = []
    names = (('I_drain', 'devsim_drain_contact_current', 'G_DD', 'DEVSIM get_contact_current at the drain (default post-solve current)'),
             ('I_circuit', 'circuit_node_current', 'G_DD', '-x[V1.I] of the drain circuit node (Newton linear-system current)'),
             ('I_weighted', 'weighting_function_current', 'G_DD', 'sum_e J_e (w_i - w_j) with w the zero-bias network solution (P1)'))
    for k, nm, q, how in names:
        g1, g01 = s['1mV'][k]['G'], s['0.1mV'][k]['G']
        gr = (100.0 * g01 - g1) / 99.0
        e = dict(name=f'{nm}{suffix}_richardson_lev{lev}', path=path_label, level=lev, same_mesh=True, mesh_id=mid, evaluation=how, definition=RICHARDSON,
                 **{f'{q}_S_per_cm': gr}, symmetric_pm1mV_S_per_cm=g1, symmetric_pm0p1mV_S_per_cm=g01,
                 rel_pm1mV=g1 / refs[q] - 1.0, rel_pm0p1mV=g01 / refs[q] - 1.0,
                 note_pm='the +-1 mV and +-0.1 mV symmetric differences contain the physical (V_D/V_T)^2/6 term (2.5e-4 and 2.5e-6)')
        out.append((e, [judge(e, q, refs[q], True, None)]))
    # DEVSIM source contact current: G_SD estimate (two-terminal B0: also G0)
    I_S = {st['V_D']: st['I_source'] for st in d['states']}
    gs1 = -(I_S[1e-3] - I_S[-1e-3]) / 2e-3
    gs01 = -(I_S[1e-4] - I_S[-1e-4]) / 2e-4
    gsr = (100.0 * gs01 - gs1) / 99.0
    e = dict(name=f'devsim_source_contact_current{suffix}_richardson_lev{lev}', path=path_label, level=lev, same_mesh=True, mesh_id=mid,
             evaluation='-DEVSIM get_contact_current at the source', definition=RICHARDSON, G_SD_S_per_cm=gsr,
             symmetric_pm1mV_S_per_cm=gs1, symmetric_pm0p1mV_S_per_cm=gs01,
             raw_source_current_A_per_cm={f'{v:+g}': I_S[v] for v in sorted(I_S)})
    out.append((e, [judge(e, 'G_SD', refs['G_SD'], True, None)]))
    return out


# ------------------------------------------------------------------------------------------- B1 three-terminal cases
def three_terminal_case(X, E, sid, files, ev, ledger, kind):
    """kind 'AR' (reference level 2, tolerance per section 7) or 'HO' (level 1, same-mesh scoring only)."""
    levels = (0, 1, 2) if kind == 'AR' else (1,)
    lref = levels[-1]
    P = {lv: X.load(p1(sid, 'B1', lv)) for lv in levels}
    B0 = {lv: X.load(p1(sid, 'B0', lv)) for lv in levels}
    r = P[lref]
    a, n = r['args'], {lv: P[lv]['network'] for lv in levels}
    L_nm = a['L'] * 1e3
    sources = [X.src(p1(sid, 'B1', lv), f'P1 network, body B1, level {lv}' + (' (reference)' if lv == lref else '')) for lv in levels]
    refs = {q: n[lref][q] for q in ('G_DD', 'G_SD', 'G_DB')}
    G_gca = n[lref]['G_GCA_1COL']
    dL_ref = 0.5 * L_nm * (1.0 - G_gca / refs['G_SD'])
    assert close(dL_ref, n[lref]['dL_SD_nm'], 1e-10)
    notes = [ZB3_NOTE]
    tols, comps, unc_status = {}, {}, 'not_established'
    if kind == 'AR':
        h1 = ev['H1'][sid]
        P4 = {m: X.load(p4(sid, m)) for m in ('refined', 'base')}
        for m in ('refined', 'base'):
            sources.append(X.src(p4(sid, m), f'P4 independent finite-volume solver ({m} mesh), bodies B1 and B0'))
        for q in ('G_DD', 'G_SD'):
            q4 = P4['refined']['results']['B1'][q]
            cross, refine = q4 / refs[q] - 1.0, refs[q] / n[1][q] - 1.0
            u_disc = max(abs(cross), abs(refine))
            u_round = abs(n[2][f'{q}_rel_forms'])
            u_solver = n[2]['residual_drain_excited']
            tol = 2.0 * (u_disc + u_round + u_solver)
            e = h1['reference'][q]
            assert close(e['value'], refs[q]) and close(e['tol_rel'], tol), (sid, q)
            comps[q] = dict(u_disc=u_disc, u_round=u_round, u_solver=u_solver, combined=u_disc + u_round + u_solver,
                            inputs=dict(cross_P4_refined_vs_L2=cross, refine_L2_vs_L1=refine,
                                        P4_refined_vs_P4_base=q4 / P4['base']['results']['B1'][q] - 1.0,
                                        observed_order_levels_0_1_2=order(n[0][q], n[1][q], n[2][q])))
            tols[q] = tol
        dl_tol = (0.5 * L_nm - dL_ref) * tols['G_SD']
        assert close(dl_tol, h1['dL_tol_nm']) and close(dL_ref, h1['dL_SD_nm'])
        if h1['pass']:
            tier, unc_status = 'A-R', 'established'
        else:
            tier = 'C'
            failed = [f"{q}: cross {c['cross_P4_vs_L2']:+.2e}, refine {c['refine_L2_vs_L1']:+.2e}" for q, c in h1['checks'].items()
                      if not (c['pass_i'] and c['pass_ii'])]
            failed += [f"level {x['level']}: forms/KCL/bracket" for x in h1['self_consistency'] if not x['ok']]
            notes.append('Failed pre-registered H1 (' + '; '.join(failed) + '): exported as tier C, tolerance not established.')
        cid = f"AR__{state_tag(a)}_B1"
    else:
        h2 = ev['H2'][sid]
        tier = 'HO'
        assert close(h2['G_DD'], refs['G_DD']) and h2['pass'] is True
        cid = f"HO__{state_tag(a)}_B1"
        notes.append('Held-out state (pre-registration section 6): level 1 only, not used to set any threshold or tolerance. '
                     'No refinement pair or independent path, so only same-mesh scoring applies.')
    jt = tols if unc_status == 'established' else {}       # other-mesh classification only against an established tolerance
    # same-mesh references on every mesh the reference path was run on
    smr = {}
    for lv in levels:
        nn = n[lv]
        smr[mesh_id(P[lv])] = dict(level=lv, G_DD_S_per_cm=nn['G_DD'], G_SD_S_per_cm=nn['G_SD'], G_DB_S_per_cm=nn['G_DB'],
                                   G_DD_rel_forms=nn['G_DD_rel_forms'], G_SD_rel_forms=nn['G_SD_rel_forms'], G_DB_rel_forms=nn['G_DB_rel_forms'],
                                   kcl_closure_rel=nn['kcl_closure_rel'], network_relative_residual=nn['residual_drain_excited'],
                                   node_lines_file=f"{MESH_PUB}/mesh_{mesh_id(P[lv])}.json")
    obs = dict(
        G_DD_S_per_cm=dict(value=refs['G_DD'], path='edge energy phi^T L phi of the drain-excited network (cancellation-free)',
                           second_form=n[lref]['G_DD_second'], second_form_definition='-(I_S + I_B), evaluated where phi = 0'),
        G_SD_S_per_cm=dict(value=refs['G_SD'], path='-I_S of the drain-excited network (source side, cancellation-free)',
                           second_form=n[lref]['G_SD_second'], second_form_definition='reciprocity: -I_D with the source excited'),
        G_DB_S_per_cm=dict(value=refs['G_DB'], path='-I_B of the drain-excited network', second_form=n[lref]['G_DB_second'],
                           second_form_definition='G_DD - G_SD (cancellation-prone)', scored=False,
                           note='reported, not scored: the pre-registered tolerance rule (section 7) covers G_DD and G_SD only'),
        dL_nm=dict(value=dL_ref, definition='per-end correction dL = (L/2)(1 - G_GCA/G_SD)'),
        G_GCA_S_per_cm=dict(value=G_gca, definition='long-channel (1-column) sheet conductance / L on the reference depth mesh; the dL reference'),
        psi_hat=dict(value=r['onecol']['psi_hat'], definition='(psi_s,inf - psi_b0)/V_T of the 1-column problem (long-channel surface potential)'),
        same_mesh_references=smr,
        scored=['G_DD_S_per_cm', 'G_SD_S_per_cm', 'dL_nm'])
    # body effect (B1 vs B0 on the reference mesh)
    g0 = B0[lref]['network']['G0']
    be = dict(G0_B0_S_per_cm=g0, Delta_DD=refs['G_DD'] / g0 - 1.0, Delta_SD=refs['G_SD'] / g0 - 1.0, G_DB_over_G0=refs['G_DB'] / g0,
              B0_case=f"BR__{state_tag(a)}_B0_lev{lref}")
    if kind == 'AR':
        hb = ev['HB'][sid]
        assert close(hb['Delta_DD'], be['Delta_DD'], 1e-9) and close(hb['Delta_SD'], be['Delta_SD'], 1e-9)
        be.update(band_DD=hb['band_DD'], band_SD=hb['band_SD'], in_band_DD=hb['in_DD'], in_band_SD=hb['in_SD'],
                  hypothesis='H-B (pre-registration section 8); never a kill')
    obs['body_effect_vs_B0'] = be
    # output paths of the project's own codes
    known_bad, cross = [], []
    mids = {lv: mesh_id(P[lv]) for lv in levels}
    for lv in levels:
        nn = n[lv]
        e = dict(name=f'drain_side_flux_sum_lev{lv}', path='P1', level=lv, same_mesh=True, mesh_id=mids[lv], G_DD_S_per_cm=nn['G_DD_direct_drain_sum'],
                 evaluation='direct sum of the edge currents at the drain contact nodes (phi = 1 there)')
        place(e, [judge(e, 'G_DD', nn['G_DD'], True, None)], known_bad, cross)
        e = dict(name=f'second_forms_lev{lv}', path='P1', level=lv, same_mesh=True, mesh_id=mids[lv], G_DD_S_per_cm=nn['G_DD_second'],
                 G_SD_S_per_cm=nn['G_SD_second'], note='-(I_S + I_B) and source-excited reciprocity on the same network')
        place(e, [judge(e, q, nn[q], True, None) for q in ('G_DD', 'G_SD')], known_bad, cross)
    ss = P[1].get('smallsignal')
    if ss:
        for k, nm, how in (('G_E', 'coupled_small_signal_G_E', 'DEVSIM full coupled Jacobian, row/column-equilibrated SuperLU solve, -x[V1.I]'),
                           ('G_M', 'coupled_small_signal_G_M', 'same Newton linear system solved with UMFPACK')):
            e = dict(name=f'{nm}_lev1', path='P3', level=1, same_mesh=True, mesh_id=mids[1], G_DD_S_per_cm=ss[k], evaluation=how,
                     kappa_G=ss['kappa_G'], E_dwr=ss['E_dwr'], E_round=ss['E_round'], eta=ss['eta'])
            place(e, [judge(e, 'G_DD', n[1]['G_DD'], True, None)], known_bad, cross)
    for lv in ((1, 2) if kind == 'AR' else (1,)):
        if not os.path.exists(os.path.join(E.ROOT, p2(sid, 'B1', lv))):
            notes.append(f'P2 (drift-diffusion) at level {lv} has no successful run: it exceeded the memory limit of the run environment '
                         '(pre-registration amendment A8; failed attempts are in benchmark/wp2/results/ledger.jsonl).')
            continue
        d = X.load(p2(sid, 'B1', lv))
        sources.append(X.src(p2(sid, 'B1', lv), f'P2 drift-diffusion at V_D = +-0.1, +-1 mV, body B1, level {lv}'))
        assert close(d['network_same_run']['G_DD'], n[lv]['G_DD'], 1e-14)
        for e, vs in p2_paths(d, dict(G_DD=n[lv]['G_DD'], G_SD=n[lv]['G_SD']), mids[lv]):
            place(e, vs, known_bad, cross)
    if kind == 'AR':
        for m in ('refined', 'base'):
            R4 = P4[m]
            e = dict(name=f'independent_finite_volume_{m}', path='P4', same_mesh=False, independent=True,
                     G_DD_S_per_cm=R4['results']['B1']['G_DD'], G_SD_S_per_cm=R4['results']['B1']['G_SD'], G_DB_S_per_cm=R4['results']['B1']['G_DB'],
                     mesh=dict(nx=R4['mesh']['nx'], ny=R4['mesh']['ny'], dx_min_nm=R4['mesh']['dx_min_nm'], dy_at_interface_nm=R4['mesh']['dy_at_interface_nm']),
                     note='benchmark/wp2/fv_rdeck.py: own mesh, Poisson-Boltzmann Newton and SG network; imports nothing from phase3/, '
                          'semiconductor_sim/ or DEVSIM. Its refined-mesh spread is part of u_disc, so passing is a consistency check.')
            place(e, [judge(e, q, refs[q], False, jt.get(q)) for q in ('G_DD', 'G_SD')], known_bad, cross)
    base_commit = WP2_COMMIT
    if kind == 'AR' and WP4_COMMIT:
        base_commit = WP4_COMMIT
        X4 = Ctx(E, WP4_COMMIT)
        mid1 = mids[1]
        refs1 = dict(G_DD=n[1]['G_DD'], G_SD=n[1]['G_SD'])
        px = f'{WP4_PROD}/P2X_{sid}_B1_lev1.json'
        d = X4.load(px)
        sources.append(X4.src(px, 'WP4 P2 with DEVSIM extended precision (float128 models/equations/Newton residual), body B1, level 1'))
        for e, vs in p2_paths(d, refs1, mid1, suffix='_extended', path_label='P2X'):
            e['precision'] = 'DEVSIM extended_solver/extended_model/extended_equation (float128 residuals, double LU)'
            place(e, vs, known_bad, cross)
        for prec in ('double', 'extended'):
            pa = f'{WP4_PROD}/AC_{sid}_B1_lev1_{prec}.json'
            d = X4.load(pa)
            sources.append(X4.src(pa, f'WP4 DEVSIM built-in AC analysis at 1e-6 Hz, {prec} precision, body B1, level 1'))
            e = dict(name=f'devsim_ac_small_signal_{prec}_lev1', path='AC', level=1, same_mesh=True, mesh_id=mid1,
                     G_DD_S_per_cm=d['ac']['G_AC'], evaluation=f"devsim.solve(type='ac', frequency=1e-6), -Re(V1.I); {prec} precision",
                     frequency_Hz=d['frequency_Hz'])
            place(e, [judge(e, 'G_DD', refs1['G_DD'], True, None)], known_bad, cross)
        notes.append('WP4 output paths on the level-1 mesh (benchmark/specs/WP4_Preregistration_External_Paths_2026-11-04.md): DEVSIM '
                     'post-solve contact currents with extended precision, and DEVSIM built-in AC analysis in double and extended precision.')
    if unc_status == 'established':
        unc = dict(status='established', coverage_factor=E.COVERAGE_K, tau_num_rel=TAU_NUM,
                   components=dict(G_DD=comps['G_DD'], G_SD=comps['G_SD']),
                   tolerance=dict(G_DD_rel=tols['G_DD'], G_SD_rel=tols['G_SD'], dL_nm=dl_tol),
                   basis=('pre-registration section 7: reference = P1 at level 2; u_disc = max(|Q_L2/Q_L1 - 1|, |Q_P4/Q_L2 - 1|); '
                          'u_round = |reference form / second form - 1|; u_solver = network relative residual; tolerance = 2 (sum), '
                          'relative, applied separately to G_DD and G_SD; dL tolerance = (L/2 - dL) * tol(G_SD). Same-mesh B1 submissions: '
                          'tau_num = 1e-6 against the P1 value on the declared mesh (observables.same_mesh_references). '
                          'The P4 spread is part of the estimate, so passing it is a consistency check, not independent validation.'))
    else:
        basis = ('held-out state: level 1 only, no refinement pair or independent path; same-mesh scoring only (tau_num = 1e-6)'
                 if kind == 'HO' else 'failed pre-registered H1; same-mesh scoring only (tau_num = 1e-6)')
        unc = dict(status='not_established', components={}, coverage_factor=E.COVERAGE_K, tolerance='not_established', tau_num_rel=TAU_NUM,
                   basis=basis)
        if kind == 'AR':
            unc['components_if_rule_applied'] = dict(G_DD=comps['G_DD'], G_SD=comps['G_SD'])
    sources.append(X.src(f'{RES}/evaluation.json', 'WP2 evaluation against the frozen pre-registration (H1, H-B, H2, P-DD, P-SS, E-CIRC)'))
    sources.append(X.src(p1(sid, 'B0', lref), 'P1 network, body B0, same mesh (body effect)'))
    job = f'P1_{sid}_B1_lev{lref}'
    return dict(
        schema_version=SCHEMA_VERSION, id=cid, tier=tier, deck='R', family='zero_bias_3terminal',
        related_cases=sorted(f"BR__{state_tag(a)}_B0_lev{lv}" for lv in levels),
        device=dict(geometry=geometry(r), doping=doping(a['na']), contacts=contacts('B1'), physics=physics(E, r['bulk']['psi_b0'], a['na']),
                    bias=dict(V_G_V=a['vg'], V_D_V=0.0, V_S_V=0.0, V_B_V=0.0), state_id=sid, body='B1'),
        mesh=mesh_obj(r, files), observables=obs, uncertainty=unc,
        known_bad_outputs=known_bad, cross_paths=cross,
        provenance=dict(base_commit=base_commit, sources=sources, producing_script='benchmark/wp2/rdeck_probe.py',
                        command=ledger[job]['cmd'] + '   (benchmark/wp2/results/ledger.jsonl)', exporter='benchmark/tools/export_cases.py',
                        pre_registration=PREREG if base_commit == WP2_COMMIT else f'{PREREG}; {PREREG_WP4}'),
        notes=notes)


# ------------------------------------------------------------------------------------------- B0 companions (tier B-R)
def b0_case(X, E, sid, lev, files, ledger):
    path = p1(sid, 'B0', lev)
    r = X.load(path)
    a, nn = r['args'], r['network']
    L_nm = a['L'] * 1e3
    G0, G_gca = nn['G0'], nn['G_GCA_1COL']
    dl = lambda G: 0.5 * L_nm * (1.0 - G_gca / G)  # noqa: E731
    lo, hi = nn['G_rows'] * (1 - E.BRACKET_RTOL), nn['G_slice'] * (1 + E.BRACKET_RTOL)
    inside = lambda G: G is not None and lo <= G <= hi  # noqa: E731
    obs = dict(
        G0_S_per_cm=dict(value=G0, path='-I_S of the two-terminal network (source side, cancellation-free); G_DD = G_SD = G0 with body B0'),
        dL_nm=dict(value=dl(G0), definition='per-end correction dL = (L/2)(1 - G_GCA/G0)'),
        G_GCA_S_per_cm=dict(value=G_gca, definition='long-channel (1-column) sheet conductance / L on this depth mesh'),
        bracket=dict(G_rows_S_per_cm=nn['G_rows'], G_slice_S_per_cm=nn['G_slice'], dL_rows_nm=dl(nn['G_rows']), dL_slice_nm=dl(nn['G_slice']),
                     rule=f'G_rows (1 - {E.BRACKET_RTOL:g}) <= G <= G_slice (1 + {E.BRACKET_RTOL:g})',
                     scope='Rayleigh bounds of THIS discrete network on THIS mesh; valid only for same-mesh submissions'),
        network_relative_residual=nn['residual_drain_excited'],
        two_terminal_aliases=dict(G_DD_S_per_cm='G0_S_per_cm', G_SD_S_per_cm='G0_S_per_cm'),
        psi_hat=dict(value=r['onecol']['psi_hat']))
    known_bad, cross = [], []

    def put(e, G):
        e['G0_S_per_cm'] = G
        e['inside_bracket'] = bool(inside(G)) if G is not None and G > 0 else False
        e['rel_vs_reference'] = G / G0 - 1.0 if G is not None else None
        e['rule'] = 'same mesh: inside the bracket' if e['same_mesh'] else 'other mesh: reported only (bracket not applicable)'
        if e['same_mesh'] and not e['inside_bracket']:
            e['evidence'] = (f"outside the same-mesh bracket [{lo:.6e}, {hi:.6e}] S/cm: G = {G!r} S/cm"
                             + (f", rel vs reference {e['rel_vs_reference']:+.3e}" if e['rel_vs_reference'] is not None else ''))
            known_bad.append(e)
        else:
            e.setdefault('independent', False)
            cross.append(e)
    put(dict(name='edge_energy', path='P1', same_mesh=True, note='phi^T L phi of the same network (G_DD reference form)'), nn['G_DD'])
    put(dict(name='drain_side_flux_sum', path='P1', same_mesh=True, evaluation='direct sum of edge currents at the drain contact (phi = 1)'),
        nn['G_DD_direct_drain_sum'])
    sources = [X.src(path, 'P1 network, body B0 (reference)')]
    ss = r.get('smallsignal')
    if ss:
        for k, nm in (('G_E', 'coupled_small_signal_G_E'), ('G_M', 'coupled_small_signal_G_M')):
            put(dict(name=nm, path='P3', same_mesh=True, kappa_G=ss['kappa_G'], E_dwr=ss['E_dwr'], E_round=ss['E_round']), ss[k])
    p2p = p2(sid, 'B0', lev)
    if os.path.exists(os.path.join(E.ROOT, p2p)):
        d = X.load(p2p)
        sources.append(X.src(p2p, f'P2 drift-diffusion at V_D = +-0.1, +-1 mV, body B0, level {lev}'))
        for e, _ in p2_paths(d, dict(G_DD=G0, G_SD=G0), mesh_id(r)):
            G = e.pop('G_DD_S_per_cm', None)
            if G is None:
                G = e.pop('G_SD_S_per_cm')
            for k in [k for k in e if k.startswith('rel_G_') or k == 'rule']:
                e.pop(k)
            put(e, G)
    if lev == 2 and sid.startswith('AR'):
        for m in ('refined', 'base'):
            R4 = X.load(p4(sid, m))
            sources.append(X.src(p4(sid, m), f'P4 independent finite-volume solver ({m} mesh)'))
            put(dict(name=f'independent_finite_volume_{m}', path='P4', same_mesh=False, independent=True,
                     mesh=dict(nx=R4['mesh']['nx'], ny=R4['mesh']['ny'])), R4['results']['B0']['G_SD'])
    tag = state_tag(a)
    job = f'P1_{sid}_B0_lev{lev}'
    return dict(
        schema_version=SCHEMA_VERSION, id=f'BR__{tag}_B0_lev{lev}', tier='B-R', deck='R', family='zero_bias_network',
        related_cases=[f"{'AR' if sid.startswith('AR') else 'HO'}__{tag}_B1"],
        device=dict(geometry=geometry(r), doping=doping(a['na']), contacts=contacts('B0'), physics=physics(E, r['bulk']['psi_b0'], a['na']),
                    bias=dict(V_G_V=a['vg'], V_D_V=0.0, V_S_V=0.0), state_id=sid, body='B0'),
        mesh=mesh_obj(r, files), observables=obs,
        uncertainty=E.not_established('B0 companion (pre-registration sections 2.3 and 10.6): scored by the same-mesh bracket only'),
        known_bad_outputs=known_bad, cross_paths=cross,
        provenance=dict(base_commit=WP2_COMMIT, sources=sources, producing_script='benchmark/wp2/rdeck_probe.py',
                        command=ledger[job]['cmd'] + '   (benchmark/wp2/results/ledger.jsonl)', exporter='benchmark/tools/export_cases.py',
                        pre_registration=PREREG),
        notes=['Two-terminal companion of the three-terminal case (electron-blocking body B0, the boundary data of every existing deck). '
               'Not separately pre-registered for reference quality; no tolerance is established, so only same-mesh (bracket) scoring applies.'])


# ------------------------------------------------------------------------------------------- driver
def export(E):
    """Returns (cases, published mesh files {published path: bytes}, index info)."""
    X = Ctx(E)
    ev = X.load(f'{RES}/evaluation.json')
    X.check(f'{RES}/ledger.jsonl')
    ledger = {}
    with open(os.path.join(E.ROOT, RES, 'ledger.jsonl'), encoding='utf-8') as f:
        for line in f:
            rec = json.loads(line)
            if rec['exit'] == 0:
                ledger[rec['job']] = rec
    # meshes: byte copies of the committed node-line files
    mesh_files = {}
    for name in sorted(os.listdir(os.path.join(E.ROOT, MESH_SRC))):
        if name.endswith('.json'):
            sp = f'{MESH_SRC}/{name}'
            X.check(sp)
            with open(os.path.join(E.ROOT, sp), 'rb') as f:
                mesh_files[f'{MESH_PUB}/{name}'] = f.read()
    import hashlib
    files = {p: hashlib.sha256(b).hexdigest() for p, b in mesh_files.items()}
    h1s, h2s = ev['H1_summary'], ev['H2_summary']
    if h2s['kill'] or not h2s['complete']:
        raise SystemExit('H2 failed or incomplete: release stops (pre-registration section 8)')
    cases = []
    ar_ids = sorted(ev['H1'])
    if h1s['kill']:
        info = 'H1 kill: fewer than 4 of 6 A-R states passed; A-R tier not released (states exported as tier C)'
    else:
        info = f"H1: {h1s['n_pass']} of 6 A-R states passed"
    for sid in ar_ids:
        c = three_terminal_case(X, E, sid, files, ev, ledger, 'AR')
        if h1s['kill'] and c['tier'] == 'A-R':
            c['tier'] = 'C'
            c['notes'].append('H1 kill criterion met for the tier as a whole: not released as A-R.')
        cases.append(c)
        for lv in (0, 1, 2):
            cases.append(b0_case(X, E, sid, lv, files, ledger))
    for sid in sorted(ev['H2']):
        cases.append(three_terminal_case(X, E, sid, files, ev, ledger, 'HO'))
        cases.append(b0_case(X, E, sid, 1, files, ledger))
    return cases, mesh_files, dict(wp2_commit=WP2_COMMIT, h1=info, h2=f"H2: {h2s['n_pass']} of {h2s['n']} held-out states passed",
                                   pre_registration=PREREG)
