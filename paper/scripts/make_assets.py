#!/usr/bin/env python3
"""Make every figure, table and quoted number of the manuscript from the benchmark's released files.

Nothing in the manuscript's results is typed by hand: each number is a LaTeX macro in paper/generated/numbers.tex,
each table a file paper/generated/tab_*.tex, each figure paper/figures/fig*.pdf, all written here from the files
listed in INPUTS. Their SHA-256 are written to paper/generated/inputs.sha256.

Usage (repository root):  python paper/scripts/make_assets.py
"""
from __future__ import annotations

import collections
import glob
import hashlib
import json
import math
import os
import statistics

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from scipy.special import erfc, erfcinv  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
ROOT = os.path.dirname(PAPER)
B = os.path.join(ROOT, 'benchmark')
GEN = os.path.join(PAPER, 'generated')
FIG = os.path.join(PAPER, 'figures')

INPUTS = {
    'index': 'benchmark/INDEX.json',
    'wp2': 'benchmark/wp2/results/evaluation.json',
    'atlas': 'benchmark/wp2/results/atlas_observations.json',
    'wp4': 'benchmark/wp4/results/evaluation_wp4.json',
    'wp6': 'benchmark/wp6/results/evaluation_wp6_stage1.json',
    'wp6_sens': 'benchmark/wp6/results/sensitivity_wp6_stage1.json',
    'wp6_expl': 'benchmark/wp6/results/exploratory_wp6_stage1.json',
    'wp6_methods': 'benchmark/wp6/results/method_coding_stage1.json',
    'wp6_verif': 'benchmark/wp6/results/verification_coding_stage1.json',
    'wp8': 'benchmark/wp8/results/evaluation_wp8.json',
    'wp8x': 'benchmark/wp8/results/exploratory/summary_xsens.json',
    'wp8x4_ar1': 'benchmark/wp8/results/exploratory/X4_AR-1_B1_lev1.json',
    'wp8x4_ar4': 'benchmark/wp8/results/exploratory/X4_AR-4_B1_lev1.json',
    'wp6_tasks': 'benchmark/wp6/tasks/tasks.json',
}

# palette: reference categorical slots 1-3 (validated all-pairs, light surface) + neutral context gray
BLUE, ORANGE, AQUA, GRAY = '#2a78d6', '#eb6834', '#1baf7a', '#7d7c78'
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e6e5e0'
plt.rcParams.update({
    'font.size': 8, 'axes.labelsize': 8, 'axes.titlesize': 8, 'legend.fontsize': 7, 'xtick.labelsize': 7,
    'ytick.labelsize': 7, 'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
    'text.color': INK, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.minor.width': 0.4, 'ytick.minor.width': 0.4, 'lines.linewidth': 1.0, 'legend.frameon': False,
    'pdf.fonttype': 42, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02, 'axes.spines.top': False,
    'axes.spines.right': False, 'mathtext.default': 'regular'})


def load(key):
    return json.load(open(os.path.join(ROOT, INPUTS[key])))


def sha(path):
    with open(path, 'rb') as f:
        return hashlib.sha256(f.read()).hexdigest()


def num(x, sig=2):
    """siunitx \\num{} with `sig` significant digits."""
    if x == 0:
        return r'\num{0}'
    e = int(math.floor(math.log10(abs(x))))
    m = round(x / 10 ** e, sig - 1)
    if abs(m) >= 10:
        m, e = m / 10, e + 1
    if -2 <= e <= 1:
        return r'\num{' + f'{m * 10 ** e:.{max(sig - 1 - e, 0)}f}' + '}'
    ms = f'{m:.{sig - 1}f}'
    return r'\num{' + f'{ms}e{e}' + '}'


def pct(k, n):
    return f'{k}/{n}'


# ----------------------------------------------------------------------------- data
def cases():
    out = {}
    for p in sorted(glob.glob(os.path.join(B, 'cases', '*.json'))):
        c = json.load(open(p))
        out[c['id']] = c
        INPUT_FILES.add(os.path.relpath(p, ROOT))
    return out


INPUT_FILES = set(INPUTS.values())


def three_terminal(cs):
    rows = []
    for cid, c in cs.items():
        if c['tier'] not in ('A-R', 'HO'):
            continue
        d, o = c['device'], c['observables']
        L_um = d['geometry']['L_um']
        dev = ('dev1' if L_um == 1.0 else 'dev1-L2') if d['doping']['N_A_cm3'] == 1e18 else 'dev2'
        tol = c['uncertainty'].get('tolerance') if c['tier'] == 'A-R' else None
        rows.append(dict(id=cid, state=d['state_id'], tier=c['tier'], dev=dev, VG=d['bias']['V_G_V'],
                         psi=o['psi_hat']['value'], G_DD=o['G_DD_S_per_cm']['value'], G_SD=o['G_SD_S_per_cm']['value'],
                         G_GCA=o['G_GCA_S_per_cm']['value'], dL=o['dL_nm']['value'], case=c, tol=tol))
    rows.sort(key=lambda r: (r['tier'] != 'A-R', r['state']))
    return rows


def tolerance_of(c, q='G_DD'):
    u = c['uncertainty']
    t = u.get('tolerance') or {}
    if isinstance(t, dict) and q in t:
        v = t[q]
        return v['relative'] if isinstance(v, dict) else v
    comp = u['components'][q]
    return 2 * comp['combined'] if 'combined' in comp else None


# ----------------------------------------------------------------------------- numbers
def write_numbers(cs, tt, wp2, atlas, wp4, wp6, sens, expl, methods):
    N = {}
    idx = load('index')
    cnt = collections.Counter(c['tier'] for c in cs.values())
    N['nCasesTotal'] = str(len(cs))
    for t, k in (('A-R', 'nAR'), ('HO', 'nHO'), ('B-R', 'nBR'), ('A-W', 'nAW'), ('B-W', 'nBW'), ('C', 'nC')):
        N[k] = str(cnt[t])
    assert sum(idx['counts'].values()) == len(cs) if isinstance(idx.get('counts'), dict) else True
    g = [r['G_DD'] for r in tt]
    N['GddMin'], N['GddMax'] = num(min(g), 2), num(max(g), 2)
    N['IoneMin'], N['IoneMax'] = num(min(g) * 1e-3, 2), num(max(g) * 1e-3, 2)
    ar = [r for r in tt if r['tier'] == 'A-R']
    tols = [tolerance_of(r['case']) for r in ar]
    N['tolMin'], N['tolMax'] = num(min(tols), 2), num(max(tols), 2)
    N['psiMin'], N['psiMax'] = f"{min(r['psi'] for r in tt):.1f}", f"{max(r['psi'] for r in tt):.1f}"
    # H1 / H2 / E-CIRC / P-DD / P-SS
    h1 = wp2['H1']
    N['HonePass'] = str(wp2['H1_summary']['n_pass'])
    N['HtwoPass'] = str(wp2['H2_summary']['n_pass'])
    N['maxPfour'] = num(max(abs(h1[s]['checks'][q]['cross_P4_vs_L2']) for s in h1 for q in ('G_DD', 'G_SD')), 2)
    N['maxRefine'] = num(max(abs(h1[s]['checks'][q]['refine_L2_vs_L1']) for s in h1 for q in ('G_DD', 'G_SD')), 2)
    N['ecircWorst'] = num(wp2['ECIRC']['ecirc_worst'], 2)
    N['nPtwoRuns'] = str(len(wp2['P2_rows']))
    pdd = wp2['PDD_by_state']
    N['pddPredicted'] = str(sum(1 for v in pdd.values() if v.get('prediction') in ('fail', 'pass')))
    N['pddConfirmed'] = str(sum(1 for v in pdd.values() if v.get('outcome') == 'confirmed'))
    N['pssConfirmed'] = str(wp2['PSS_summary']['n_confirmed'])
    N['pssPredicted'] = str(wp2['PSS_summary']['n_predicted'])
    b1 = [r for r in wp2['PSS_rows'] if r['body'] == 'B1']
    N['pssMaxGE'] = num(max(abs(r['G_E_rel']) for r in b1), 2)
    N['pssMaxGM'] = num(max(abs(r['G_M_rel']) for r in b1), 2)
    N['kappaMin'], N['kappaMax'] = num(min(r['kappa_G'] for r in b1), 2), num(max(r['kappa_G'] for r in b1), 2)
    N['etaMax'] = num(max(r['eta'] for r in b1), 1)
    nf = atlas['noise_floor_by_mesh']
    errs = [x for m in nf.values() for k in ('max_abs_err_drain_range', 'max_abs_err_source_range') for x in m[k]]
    N['floorMin'], N['floorMax'] = num(min(errs), 2), num(max(errs), 2)
    ident = [x for x in atlas['contact_current_identity'] if x['mesh'] == 'dev1_L1_lev1' and x['body'] == 'B1'
             and x['V_D'] == 0.0001 and x['current'] == 'I_drain'][0]
    assert ident['n_distinct'] == 1
    N['identValue'] = num(ident['modal_value_A_per_cm'], 3)
    N['identN'] = str(ident['n_runs'])
    N['identGlo'], N['identGhi'] = num(ident['G_DD_range_of_sharing_states'][0], 2), num(ident['G_DD_range_of_sharing_states'][1], 2)
    N['identStates'] = ', '.join(ident['states_sharing_modal_value'])
    for mid, key in (('dev1_L1_lev1', 'meshDevOneLevOne'), ('dev2_L2_lev1', 'meshDevTwoLevOne'),
                     ('dev1_L1_lev2', 'meshDevOneLevTwo'), ('dev2_L2_lev2', 'meshDevTwoLevTwo')):
        mp = os.path.join(B, 'meshes', f'mesh_{mid}.json')
        INPUT_FILES.add(os.path.relpath(mp, ROOT))
        m = json.load(open(mp))
        N[key] = rf"${m['nx']} \times {m['ny']}$"
    # H3 and WP4
    h3 = wp4['H3']
    n_fail = sum(v.get('fail', 0) for v in h3['by_tier'].values())
    n_all = sum(sum(v.values()) for v in h3['by_tier'].values())
    N['hthreeFail'], N['hthreeAll'] = str(n_fail), str(n_all)
    B1r, B2r = wp4['B1'], wp4['B2']
    N['extPass'] = str(sum(1 for r in B1r if r['X1'] and r['X2'] and r['X3']))
    N['extN'] = str(len(B1r))
    N['extMaxRich'] = num(max(abs(r['default_richardson']) for r in B1r), 2)
    N['extMaxErr'] = num(max(max(r['e_D'], r['e_S']) for r in B1r), 2)
    cost = [r['wall_s'] / r['wall_s_double'] for r in B1r]
    N['extCostMin'], N['extCostMax'] = f'{min(cost):.1f}', f'{max(cost):.1f}'
    N['acMaxDouble'] = num(max(abs(r['double']['rel']) for r in B2r), 2)
    N['acMaxExt'] = num(max(abs(r['extended']['rel']) for r in B2r), 2)
    N['acFailDouble'] = str(sum(1 for r in B2r if abs(r['double']['rel']) > 1e-6))
    N['acFailExt'] = str(sum(1 for r in B2r if abs(r['extended']['rel']) > 1e-6))
    N['acN'] = str(len(B2r))
    N['acBitwise'] = str(sum(1 for r in B2r if r['double']['rel'] == r['P3_G_M_rel_WP2']))
    dmax = [abs(r['default_richardson_double']) for r in B1r]
    N['defMaxDouble'] = num(max(dmax), 2)
    # agent study
    s = wp6['summary']['A']
    N['agentRuns'], N['agentValid'] = str(s['n_runs']), str(s['n_valid'])
    for k in ('SE', 'CS', 'FE', 'FC', 'NR'):
        N['agent' + k] = str(s['counts'][k])
    N['agentSEhi'] = f"{s['SE_CI95'][1]:.3f}"
    N['agentCSlo'], N['agentCShi'] = f"{s['CS_CI95'][0]:.2f}", f"{s['CS_CI95'][1]:.2f}"
    N['agentCSrate'] = f"{s['CS_rate']:.2f}"
    N['agentMedianWall'] = f"{s['median_wall_s']:.0f}"
    N['agentCost'] = f"{s['total_cost_usd']:.2f}"
    runs = wp6['runs']
    rel = sorted(abs(r['rel_error']) for r in runs if r['rel_error'] is not None)
    N['agentMaxRel'], N['agentMedianRel'] = num(max(rel), 2), num(statistics.median(rel), 2)
    N['agentNreported'] = str(len(rel))
    N['agentLow'] = str(sum(1 for r in runs if r['task_class'] == 'low'))
    N['agentModerate'] = str(sum(1 for r in runs if r['task_class'] == 'moderate'))
    h = expl['method_hints']
    N['agentCalledContact'], N['agentCalledAC'] = str(len(h['get_contact_current'])), str(len(h['ac_analysis']))
    N['agentCalledCircuit'] = str(len(h['circuit_node_value']))
    N['agentExtended'], N['agentJacobian'] = str(len(h['extended_precision'])), str(len(h['jacobian_export']))
    cats = collections.Counter(v['category'] for v in methods['runs'].values())
    for k in 'ABCD':
        N['method' + k] = str(cats[k])
    vk = collections.Counter(v['kind'] for v in load('wp6_verif')['runs'].values())
    N['verifCurrent'] = str(vk['cuts'] + vk['network'] + vk['contacts'])
    N['verifSheet'] = str(vk['sheet_density'])
    assert vk['cuts'] + vk['network'] + vk['contacts'] + vk['sheet_density'] == len(rel)
    fl = [(r['psi'], kb['rel_G_DD_vs_reference']) for r in tt for kb in r['case']['known_bad_outputs']
          if kb['name'].startswith('drain_side_flux_sum')]
    worst12 = max((v for p_, v in fl if abs(p_ - 12) < 0.1), key=abs)
    assert worst12 < 0
    N['fluxMax'] = num(abs(worst12), 2)
    bc = expl['by_class']
    N['agentWallLow'], N['agentWallModerate'] = f"{bc['low']['median_wall_s']:.0f}", f"{bc['moderate']['median_wall_s']:.0f}"
    v = sens['variants']
    worst = v['worst case, original audit: contaminated = SE']['A']
    N['sensWorstSE'] = str(worst['counts']['SE'])
    N['sensWorstLo'], N['sensWorstHi'] = f"{worst['SE_CI95'][0]:.3f}", f"{worst['SE_CI95'][1]:.3f}"
    N['sensOrigValid'] = str(v['original audit (frozen 710e55a)']['A']['n_valid'])
    N['sensOrigContaminated'] = str(len(sens['runs_whose_validity_differs']))
    rc = expl['refused_calls']
    N['refusedCalls'] = str(sum(rc['by_mechanism'].values()))
    N['refusedRuns'] = str(len(rc['runs_with_refusals']))
    N['refusedClassifier'] = str(sum(v for k, v in rc['by_mechanism'].items() if k.startswith('classifier')))
    write_numbers_wp8(N, load('wp8'), load('wp8x'), {r['state']: r['dev'] for r in tt})
    with open(os.path.join(GEN, 'numbers.tex'), 'w') as f:
        f.write('% generated by paper/scripts/make_assets.py from the benchmark release files; do not edit\n')
        for k in sorted(N):
            f.write(f'\\newcommand{{\\{k}}}{{{N[k]}}}\n')
    return N


def write_numbers_wp8(N, w8, w8x, devof):
    """WP8: mechanism of the small-signal failure (M0-M2), grid convergence of the A-R references, and the exploratory
    first-order sensitivity analysis X (pre-registration amendment 1)."""
    m, ms = w8['mechanism']['rows'], w8['mechanism']['summary']
    g, gs = w8['grid']['rows'], w8['grid']['summary']
    gate = [x for x in m if x['in_gate']]
    N['ssStates'] = str(len(m))
    N['ssRepro'] = str(ms['n_reproducing'])
    N['ssGate'] = str(ms['n_in_gate'])
    N['ssPMoneOk'], N['ssPMtwoOk'] = str(ms['P_M1']['n_ok']), str(ms['P_M2']['n_ok'])
    N['ssConverged'] = str(sum(x['M1_converged'] for x in m))
    N['ssIterMax'] = str(max(x['M1_iterations'] for x in m))
    N['ssOmegaMax'] = num(max(x['omega_dd'] for x in m), 2)
    N['ssMbMax'] = num(max(abs(x['M1b_rel_vs_G_IR']) for x in m), 1) if all(x['M1b_rel_vs_G_IR'] is not None for x in m) \
        and max(abs(x['M1b_rel_vs_G_IR']) for x in m) > 0 else r'\num{0}'
    N['ssExactMin'], N['ssExactMax'] = num(min(x['e_IR'] for x in m), 2), num(max(x['e_IR'] for x in m), 2)
    N['ssSpreadMin'], N['ssSpreadMax'] = num(min(x['S'] for x in m), 2), num(max(x['S'] for x in m), 2)
    ratio = [x['S'] / x['e_IR'] for x in m if x['e_IR'] > 0]
    N['ssRatioMin'], N['ssRatioMax'] = f'{min(ratio):.1f}', f'{max(ratio):.1f}'
    absd = {'dev1': [], 'dev2': [], 'dev1-L2': []}
    for x in m:
        absd[devof[x['state']]].append(x['e_IR'] * x['G_net'])
    for k, key in (('dev1', 'One'), ('dev2', 'Two'), ('dev1-L2', 'Three')):
        N[f'ssAbs{key}Min'], N[f'ssAbs{key}Max'] = num(min(absd[k]), 2), num(max(absd[k]), 2)
    share = {d: [x['e_IR'] / x['e_dbl'] for x in gate if devof[x['state']] == d] for d in absd}
    def p2(v):                                    # percent, two significant digits
        return f'{float(f"{100 * v:.2g}"):g}'
    N['ssShareOneMin'], N['ssShareOneMax'] = p2(min(share['dev1'])), p2(max(share['dev1']))
    N['ssShareTwoMin'], N['ssShareTwoMax'] = p2(min(share['dev2'])), p2(max(share['dev2']))
    N['ssShareThree'] = p2(share['dev1-L2'][0])
    sr = [x['S'] / x['e_dbl'] for x in gate]
    N['ssSdblMin'], N['ssSdblMax'] = p2(min(sr)), p2(max(sr))
    N['ssEirSix'] = num(next(x['e_IR'] for x in m if x['state'] == 'AR-6'), 2)
    N['ssResolved'] = str(ms['n_resolved_1e6'])
    N['ssGateShareMin'] = num(min(x['e_IR'] / x['e_dbl'] for x in gate), 2)
    mono = [x for x in g if x['monotone']]
    N['gciPairs'], N['gciMonotone'] = str(len(g)), str(len(mono))
    N['gciPmin'], N['gciPmax'] = f"{min(x['p'] for x in mono):.2f}", f"{max(x['p'] for x in mono):.2f}"
    N['gciAsymMin'] = f"{min(x['asymptotic_ratio'] for x in mono):.3f}"
    N['gciAsymMax'] = f"{max(x['asymptotic_ratio'] for x in mono):.3f}"
    N['gciFineWithin'] = str(sum(x['GCI_fine'] <= x['tolerance'] for x in mono))
    N['gciRelMin'], N['gciRelMax'] = num(min(abs(x['rel_f1_vs_ext']) for x in mono), 2), num(max(abs(x['rel_f1_vs_ext']) for x in mono), 2)
    N['gciFineMin'], N['gciFineMax'] = num(min(x['GCI_fine'] for x in mono), 2), num(max(x['GCI_fine'] for x in mono), 2)
    tr = [x['tolerance'] / abs(x['rel_f1_vs_ext']) for x in mono]
    N['gciTolRatioMin'], N['gciTolRatioMax'] = f'{min(tr):.1f}', f'{max(tr):.1f}'
    N['gciPfourMax'] = num(max(abs(x['P4_rel_vs_ext']) for x in mono if x['P4_rel_vs_ext'] is not None), 2)
    N['gciWithin'] = str(gs['P_G3']['n_within'])
    dd = [x for x in mono if x['quantity'] == 'G_DD' and x['P4_rel_vs_ext'] is not None]
    N['gciPfourCloser'] = str(sum(abs(x['P4_rel_vs_ext']) < abs(x['rel_f1_vs_ext']) for x in dd))
    N['gciStates'] = str(len(dd))
    # exploratory analysis X
    N['xsStates'] = str(len(w8x))
    N['xsKappaMin'], N['xsKappaMax'] = num(min(x['kappa_A'] for x in w8x), 2), num(max(x['kappa_A'] for x in w8x), 2)
    fl = {}
    for x in w8x:
        fl.setdefault(devof[x['state']], []).append(x['abs_rms_u'])
    for d, key in (('dev1', 'One'), ('dev2', 'Two'), ('dev1-L2', 'Three')):
        N[f'xsFloor{key}'] = num(sum(fl[d]) / len(fl[d]), 2)
    N['xsFloorSpread'] = num(max(max(v) / min(v) - 1 for v in fl.values()), 2)
    means = [sum(v) / len(v) for v in fl.values()]
    N['xsFloorMin'], N['xsFloorMax'] = num(min(means), 2), num(max(means), 2)
    N['xsLinearMax'] = num(max(abs(r) for x in w8x for r in x['rel_dG1']), 1)
    N['xsSecondMax'] = num(max(abs(a / b) for x in w8x for a, b in zip(x['dG2'], x['shifts'])), 1)
    N['xsSameMtwo'] = str(sum(all(x['same_as_M2']) for x in w8x))
    N['xsDrainMin'] = f"{100 * min(x['by_block_and_class']['si_electrons|n_drain_side']['share'] for x in w8x):.2f}"
    x4 = [load('wp8x4_ar1'), load('wp8x4_ar4')]
    N['xfRowSumMax'] = num(max(r['rel_rowsum_psi']['max'] for r in x4), 2)
    N['xfLamDev'] = num(max(max(abs(r['lambda_region']['min'] - 1), abs(r['lambda_region']['max'] - 1)) for r in x4), 1)
    N['xfRows'] = ', '.join(str(r['region_rows']) for r in x4)
    N['xsRowsEleMin'] = f"{100 * min(x['by_row_block']['si_electrons'] for x in w8x):.1f}"
    N['xsColsPotMin'] = f"{100 * min(x['by_column_block']['si_potential'] for x in w8x):.2f}"
    N['xsUnrefMax'] = num(max(abs(x['unrefined_abs_rms_u'] / x['abs_rms_u'] - 1) for x in w8x), 2)
    act = [abs(x['abs_err_IR']) / x['abs_rms_u'] for x in w8x]
    N['xsActMin'], N['xsActMax'] = num(min(act), 2), num(max(act), 2)
    N['xsLamDevMax'] = num(max(abs(t['lambda_i'] - 1) for x in w8x for t in x['top']), 1)
    N['xsEntryMin'] = num(min(t['max_abs_entry'] for x in w8x for t in x['top']), 1)
    N['xsEntryMax'] = num(max(t['max_abs_entry'] for x in w8x for t in x['top']), 2)
    N['xsWpTwoMin'] = num(min(x['kappa_G_wp2'] for x in w8x), 1)
    N['xsWpTwoMax'] = num(max(x['kappa_G_wp2'] for x in w8x), 1)


def write_table_gci(w8):
    L = [r'\begin{tabular}{@{}lrrrrr@{}}', r'\toprule',
         r'State & $p$ & $|G_{L2}/G_\mathrm{ext} - 1|$ & $\mathrm{GCI}_\mathrm{fine}$ & Tolerance & $|G_{P4}/G_\mathrm{ext} - 1|$ \\',
         r'\midrule']
    for x in w8['grid']['rows']:
        if x['quantity'] != 'G_DD':
            continue
        L.append(f"{x['state']} & {x['p']:.2f} & {num(abs(x['rel_f1_vs_ext']), 2)} & {num(x['GCI_fine'], 2)} & "
                 f"{num(x['tolerance'], 2)} & {num(abs(x['P4_rel_vs_ext']), 2)} \\\\")
    L += [r'\bottomrule', r'\end{tabular}']
    open(os.path.join(GEN, 'tab_gci.tex'), 'w').write('\n'.join(L) + '\n')


# ----------------------------------------------------------------------------- tables
def write_tables(cs, tt, wp6, expl):
    # states
    L = [r'\begin{tabular}{@{}llrrrrrr@{}}', r'\toprule',
         r'State & Device & $V_G$ (V) & $\hat\psi$ & $G_{DD}$ (S/cm) & $G_{SD}$ (S/cm) & $\Delta L$ (nm) & Tolerance \\',
         r'\midrule']
    for r in tt:
        tol = num(tolerance_of(r['case']), 2) if r['tier'] == 'A-R' else r'$\tau_\mathrm{num}$ only'
        devl = {'dev1': 'dev1', 'dev2': 'dev2', 'dev1-L2': r'dev1, $L = 2\,\mu$m'}[r['dev']]
        L.append(f"{r['state']} & {devl} & {r['VG']:.6g} & {r['psi']:.1f} & {num(r['G_DD'], 5)} & {num(r['G_SD'], 5)} & "
                 f"{r['dL']:.2f} & {tol} \\\\")
        if r['state'] == 'AR-6':
            L.append(r'\midrule')
    L += [r'\bottomrule', r'\end{tabular}']
    open(os.path.join(GEN, 'tab_states.tex'), 'w').write('\n'.join(L) + '\n')

    # scored examples
    names = [
        ('three_terminal_reference_same_mesh', 'Reference forms (A-R, HO), same mesh'),
        ('three_terminal_circuit_node_current_same_mesh', r'DEVSIM circuit-node current $-x[\mathrm{V1.I}]$, Richardson, same mesh'),
        ('three_terminal_independent_finite_volume_other_mesh', 'Independent finite volume (A-R), other mesh'),
        ('three_terminal_devsim_drain_contact_current_same_mesh', r'DEVSIM \texttt{get\_contact\_current}, double, Richardson, same mesh'),
        ('three_terminal_devsim_contact_currents_extended_level1_mesh', r'DEVSIM \texttt{get\_contact\_current}, extended precision (A-R)'),
        ('three_terminal_coupled_small_signal_G_E_level1_mesh', r'Coupled small-signal solve $G_E$ (full Jacobian), level-1 mesh'),
        ('three_terminal_devsim_ac_small_signal_double_level1_mesh', r"DEVSIM AC analysis, double (A-R)"),
        ('three_terminal_devsim_ac_small_signal_extended_level1_mesh', r"DEVSIM AC analysis, extended precision (A-R)"),
        ('source_side_reference_same_mesh', 'Two-terminal source-side reference, same mesh'),
        ('drain_side_flux_sum_same_mesh', 'Two-terminal drain-side flux sum, same mesh'),
        ('coupled_small_signal_G_E_same_mesh', r'Two-terminal coupled small-signal $G_E$, same mesh'),
        ('independent_finite_volume_other_mesh', 'Two-terminal independent finite volume (A-W), other mesh'),
    ]
    L = [r'\begin{tabular}{@{}lrr@{}}', r'\toprule', r'Submitted evaluation path & Pass & Fail \\', r'\midrule']
    for i, (fn, lab) in enumerate(names):
        p = os.path.join(B, 'examples', fn + '.report.json')
        INPUT_FILES.add(os.path.relpath(p, ROOT))
        rep = json.load(open(p))
        res = rep['results']
        items = res.values() if isinstance(res, dict) else res
        c = collections.Counter(x['status'] for x in items)
        L.append(f"{lab} & {c.get('pass', 0)} & {c.get('fail', 0)} \\\\")
        if i == 7:
            L.append(r'\midrule')
    L += [r'\bottomrule', r'\end{tabular}']
    open(os.path.join(GEN, 'tab_examples.tex'), 'w').write('\n'.join(L) + '\n')

    # agent study per task
    L = [r'\begin{tabular}{@{}lllrlr@{}}', r'\toprule',
         r'Task & State & Class & $G_\mathrm{ref}$ (S/cm) & Outcomes & max $|\epsilon|$ \\', r'\midrule']
    for x in expl['per_task']:
        m = num(x['max_abs_rel_error'], 2) if x['max_abs_rel_error'] is not None else '--'
        L.append(f"{x['task']} & {x['state_id']} & {x['task_class']} & {num(x['G_ref'], 4)} & {' '.join(x['outcomes'])} & {m} \\\\")
    L += [r'\bottomrule', r'\end{tabular}']
    open(os.path.join(GEN, 'tab_agent.tex'), 'w').write('\n'.join(L) + '\n')


# ----------------------------------------------------------------------------- figures
def fig_device(cs):
    c = cs['AR__R_L1_NA1e18_tox2_VG0.488174_B1']
    d = c['device']
    L, lw, H = d['geometry']['L_um'], d['geometry']['l_w_um'], d['geometry']['si_depth_um']
    hw, NA, ND = d['doping']['h_w_um'], d['doping']['N_A_cm3'], d['doping']['N_D_peak_cm3']
    sg, sgy = d['doping']['sigma_nm'] * 1e-3, d['doping']['sigma_y_nm'] * 1e-3
    xc = -math.sqrt(2) * sg * erfcinv(2 * NA / ND)
    x = np.linspace(-lw, L + lw, 801)
    y = np.linspace(0, H, 601)
    X, Y = np.meshgrid(x, y)
    fx = 0.5 * erfc((X - xc) / (math.sqrt(2) * sg)) + 0.5 * erfc(((L - X) - xc) / (math.sqrt(2) * sg))
    fy = 0.5 * erfc((Y - hw) / (math.sqrt(2) * sgy))
    nd = ND * fx * fy
    net = np.log10(np.maximum(nd, 1.0))
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    im = ax.pcolormesh(X, Y, net, cmap='Blues', vmin=17, vmax=20.2, shading='auto', rasterized=True)
    ax.contour(X, Y, nd - NA, levels=[0], colors=[INK], linewidths=0.8)
    tox = 0.05   # drawn thickness of the oxide (not to scale)
    ax.add_patch(plt.Rectangle((-lw, -tox), L + 2 * lw, tox, facecolor='#f3e9c6', edgecolor=INK2, lw=0.5))
    ax.plot([0, L], [-tox, -tox], color=ORANGE, lw=3, solid_capstyle='butt')
    for xx, lab in ((-lw, 'source'), (L + lw, 'drain')):
        ax.plot([xx, xx], [0, hw], color=ORANGE, lw=3, solid_capstyle='butt')
        ax.text(xx + (-0.04 if lab == 'source' else 0.04), hw / 2, lab, rotation=90, va='center',
                ha='right' if lab == 'source' else 'left', fontsize=7)
    ax.plot([-lw, L + lw], [H, H], color=ORANGE, lw=3, solid_capstyle='butt')
    ax.text(L / 2, H + 0.035, 'body contact (B1)', ha='center', va='top', fontsize=7)
    ax.text(L / 2, -tox - 0.035, 'gate', ha='center', va='bottom', fontsize=7)
    ax.text(-0.15, 0.25, r'n$^+$', ha='center', va='center', fontsize=8)
    ax.text(L + 0.15, 0.25, r'n$^+$', ha='center', va='center', fontsize=8)
    ax.text(L / 2, 0.8, rf'p, $N_A = 10^{{18}}$ cm$^{{-3}}$', ha='center', va='center', fontsize=7)
    ax.annotate('oxide (2 nm; not to scale)', xy=(0.2, -tox / 2), xytext=(0.2, -0.22), fontsize=6.5,
                ha='center', arrowprops=dict(arrowstyle='-', lw=0.5, color=INK2), color=INK2)
    ax.set_xlim(-lw - 0.13, L + lw + 0.13)
    ax.set_ylim(H + 0.12, -0.3)
    ax.set_xlabel(r'$x$ ($\mu$m)')
    ax.set_ylabel(r'$y$ ($\mu$m)')
    ax.set_aspect('equal')
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label(r'$\log_{10} N_D$ (cm$^{-3}$)', fontsize=7)
    cb.ax.tick_params(labelsize=6.5)
    cb.outline.set_linewidth(0.4)
    fig.savefig(os.path.join(FIG, 'fig1_device.pdf'), dpi=300, metadata={'CreationDate': None})
    plt.close(fig)


def fig_references(tt):
    fig, ax = plt.subplots(figsize=(3.3, 2.3))
    for dev, col in (('dev1', BLUE), ('dev2', ORANGE)):
        for tier, fill in (('A-R', True), ('HO', False)):
            rr = sorted([r for r in tt if r['dev'] == dev and r['tier'] == tier], key=lambda r: r['psi'])
            ax.plot([r['psi'] for r in rr], [r['G_DD'] for r in rr], ls='none', marker='o', ms=5,
                    mfc=col if fill else 'white', mec=col, mew=1.0)
        rr = sorted([r for r in tt if r['dev'] == dev], key=lambda r: r['psi'])
        ax.plot([r['psi'] for r in rr], [r['G_DD'] for r in rr], color=col, lw=0.8, zorder=0)
    v = [r for r in tt if r['dev'] == 'dev1-L2']
    ax.plot([r['psi'] for r in v], [r['G_DD'] for r in v], ls='none', marker='D', ms=4.5, mfc='white', mec=BLUE, mew=1.0)
    for r in v:
        ax.annotate(rf"{r['state']}: dev1, $L = 2\,\mu$m", xy=(r['psi'], r['G_DD']), xytext=(r['psi'] + 1.6, r['G_DD'] * 0.6),
                    fontsize=6.3, color=INK2, va='center', ha='left', arrowprops=dict(arrowstyle='-', lw=0.5, color=INK2))
    ax.set_yscale('log')
    ax.set_xlabel(r'$\hat\psi = (\psi_{s,\infty} - \psi_{b0})/V_T$')
    ax.set_ylabel(r'reference $G_{DD}$ (S/cm)')
    ax.set_xlim(10.5, 22)
    ax.set_ylim(3e-13, 1e-4)
    ax.grid(True, which='major', color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    h = [Line2D([], [], color=BLUE, lw=0.8, label='dev1'), Line2D([], [], color=ORANGE, lw=0.8, label='dev2'),
         Line2D([], [], ls='none', marker='o', mfc=INK2, mec=INK2, ms=5, label='A-R'),
         Line2D([], [], ls='none', marker='o', mfc='white', mec=INK2, ms=5, label='HO')]
    ax.legend(handles=h, loc='upper left', fontsize=6.5, ncol=2, columnspacing=1.0)
    fig.savefig(os.path.join(FIG, 'fig2_references.pdf'), metadata={'CreationDate': None})
    plt.close(fig)


def fig_atlas(cs, tt, wp2, wp4):
    G = {r['state']: r['G_DD'] for r in tt}
    p2 = [r for r in wp2['P2_rows'] if r['level'] == 1 and r['body'] == 'B1']
    flux = []
    for r in tt:
        for kb in r['case']['known_bad_outputs']:
            if kb['name'] == 'drain_side_flux_sum_lev1':
                flux.append((G[r['state']] * 1e-3, abs(kb['rel_G_DD_vs_reference'])))
    floor = 1e-11

    def clip(v):
        return max(abs(v), floor)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.6, 2.7), sharey=True)
    a1.plot([x for x, _ in flux], [clip(v) for _, v in flux], ls='none', marker='x', ms=5, mew=1.0, color=GRAY,
            label='flux sum at drain (network)')
    a1.plot([G[r['state']] * 1e-3 for r in p2], [clip(r['default_richardson']) for r in p2], ls='none', marker='v', ms=5,
            mfc=ORANGE, mec='white', mew=0.5, label='DEVSIM contact, double')
    a1.plot([G[r['state']] * 1e-3 for r in wp4['B1']], [clip(r['default_richardson']) for r in wp4['B1']], ls='none',
            marker='^', ms=6.5, mfc='none', mec=AQUA, mew=1.1, label='DEVSIM contact, extended')
    a1.plot([G[r['state']] * 1e-3 for r in p2], [clip(r['circuit_richardson']) for r in p2], ls='none', marker='o', ms=4.5,
            mfc=BLUE, mec='white', mew=0.5, label=r'DEVSIM V1.I, double')
    a1.set_title('(a) post-solve drain current', loc='left')
    pss = [r for r in wp2['PSS_rows'] if r['body'] == 'B1']
    a2.plot([G[r['state']] * 1e-3 for r in pss], [clip(r['G_M_rel']) for r in pss], ls='none', marker='o', ms=6.5, mfc='none',
            mec=BLUE, mew=1.0, label=r'full Jacobian $G_M$, double')
    a2.plot([G[r['state']] * 1e-3 for r in wp4['B2']], [clip(r['double']['rel']) for r in wp4['B2']], ls='none', marker='v',
            ms=4.5, mfc=ORANGE, mec='white', mew=0.5, label=r"DEVSIM AC, double")
    a2.plot([G[r['state']] * 1e-3 for r in wp4['B2']], [clip(r['extended']['rel']) for r in wp4['B2']], ls='none', marker='^',
            ms=4.5, mfc=AQUA, mec='white', mew=0.5, label='DEVSIM AC, extended')
    w8 = load('wp8')['mechanism']['rows']
    a2.plot([x['G_net'] * 1e-3 for x in w8], [clip(x['e_IR']) for x in w8], ls='none', marker='D', ms=4.0, mfc=GRAY,
            mec='white', mew=0.5, label='exact solve of the double system')
    devof = {r['state']: r['dev'] for r in tt}
    fl = {}
    for x in load('wp8x'):
        fl.setdefault(devof[x['state']], []).append(x['abs_rms_u'])
    Ix = np.logspace(np.log10(8e-16), np.log10(5e-8), 50)
    F1, F3 = sum(fl['dev1']) / len(fl['dev1']), sum(fl['dev1-L2']) / len(fl['dev1-L2'])
    assert abs(F3 / F1 - 1) < 1e-3                          # HO-5 shares the dev1 drain region and floor
    floor_lines = []
    for d, ls, lab in (('dev1', '-', 'dev1 and HO-5'), ('dev2', (0, (4, 2)), 'dev2')):
        F = sum(fl[d]) / len(fl[d])
        floor_lines += a2.plot(Ix, F * 1e-3 / Ix, color=GRAY, lw=0.8, ls=ls, zorder=1, label=lab)
    leg = a2.legend(handles=floor_lines, loc='lower left', fontsize=6.3, handlelength=2.6, title='one-ulp floor $u\\|c\\|_2/G$',
                    title_fontsize=6.3)
    leg._legend_box.align = 'left'
    a2.add_artist(leg)
    a2.set_title('(b) coupled small-signal solve', loc='left')
    for ax, xa in ((a1, 0.44), (a2, 0.58)):
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylim(5e-12, 3e4)
        ax.set_xlim(8e-16, 5e-8)
        ax.axhline(1e-6, color=INK2, lw=0.6, ls='--', zorder=0)
        ax.axhline(1e-2, color=INK2, lw=0.6, ls=':', zorder=0)
        ax.text(1.1e-15, 1.4e-6, r'$\tau_\mathrm{num} = 10^{-6}$', fontsize=6.5, color=INK2, va='bottom')
        ax.text(1.1e-15, 1.4e-2, '1%', fontsize=6.5, color=INK2, va='bottom')
        ax.grid(True, which='major', color=GRID, lw=0.5)
        ax.set_axisbelow(True)
        ax.set_xlabel(r'reference current at 1 mV, $G_{DD}\times 1$ mV (A/cm)')
        hs = [h for h in ax.get_lines() if not h.get_label().startswith('_') and h not in floor_lines]
        ax.legend(handles=hs, loc='upper center', bbox_to_anchor=(xa, -0.2), ncol=2, fontsize=6.3, handletextpad=0.3,
                  columnspacing=1.0)
    a1.set_ylabel(r'$|G/G_\mathrm{ref} - 1|$ (same mesh, level 1)')
    fig.savefig(os.path.join(FIG, 'fig3_atlas.pdf'), metadata={'CreationDate': None})
    plt.close(fig)


def fig_agents(wp6, methods, tasks):
    tref = {t['task']: t['reference_G_DD_S_per_cm'] for t in tasks['tasks']}
    cat = {k: v['category'] for k, v in methods['runs'].items()}
    style = {'A': (BLUE, 'o', 'central differences, double'), 'B': (ORANGE, 'v', 'central differences, extended'),
             'C': (AQUA, 's', 'Scharfetter-Gummel network')}
    fig, ax = plt.subplots(figsize=(3.3, 2.4))
    seen = collections.defaultdict(int)
    for r in sorted(wp6['runs'], key=lambda r: r['run_id']):
        x0 = tref[r['task']] * 1e-3
        j = seen[r['task']]
        seen[r['task']] += 1
        x = x0 * 10 ** ((j - 1) * 0.07)
        if r['rel_error'] is None:
            ax.plot([x], [2e-2], ls='none', marker='x', ms=5, mew=1.0, color=GRAY)
            ax.annotate('no report\n(time limit)', xy=(x, 2e-2), xytext=(x * 0.08, 2.6e-3), fontsize=6, color=INK2,
                        ha='center', arrowprops=dict(arrowstyle='-', lw=0.5, color=INK2))
            continue
        col, mk, _ = style[cat[r['run_id']]]
        ax.plot([x], [max(abs(r['rel_error']), 1e-8)], ls='none', marker=mk, ms=4.5, mfc=col, mec='white', mew=0.5)
    ax.axhline(1e-2, color=INK2, lw=0.6, ls=':')
    ax.text(1.2e-15, 1.25e-2, 'required accuracy (1%)', fontsize=6.5, color=INK2, va='bottom')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(3e-7, 6e-2)
    ax.set_xlim(8e-16, 5e-8)
    ax.set_xlabel(r'task: reference current at 1 mV (A/cm)')
    ax.set_ylabel(r'reported $|G/G_\mathrm{ref} - 1|$')
    ax.grid(True, which='major', color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    h = [Line2D([], [], ls='none', marker=m, mfc=c, mec='white', ms=5, label=lab) for c, m, lab in style.values()]
    ax.legend(handles=h, loc='upper left', bbox_to_anchor=(0.0, 0.86), fontsize=6.3, title='final method',
              title_fontsize=6.3, alignment='left')
    fig.savefig(os.path.join(FIG, 'fig4_agents.pdf'), metadata={'CreationDate': None})
    plt.close(fig)


def main():
    os.makedirs(GEN, exist_ok=True)
    os.makedirs(FIG, exist_ok=True)
    cs = cases()
    tt = three_terminal(cs)
    wp2, atlas, wp4 = load('wp2'), load('atlas'), load('wp4')
    wp6, sens, expl, methods, tasks = load('wp6'), load('wp6_sens'), load('wp6_expl'), load('wp6_methods'), load('wp6_tasks')
    N = write_numbers(cs, tt, wp2, atlas, wp4, wp6, sens, expl, methods)
    write_tables(cs, tt, wp6, expl)
    write_table_gci(load('wp8'))
    fig_device(cs)
    fig_references(tt)
    fig_atlas(cs, tt, wp2, wp4)
    fig_agents(wp6, methods, tasks)
    with open(os.path.join(GEN, 'inputs.sha256'), 'w') as f:
        for p in sorted(INPUT_FILES):
            f.write(f'{sha(os.path.join(ROOT, p))}  {p}\n')
    print(json.dumps(dict(numbers=len(N), inputs=len(INPUT_FILES)), indent=0))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
