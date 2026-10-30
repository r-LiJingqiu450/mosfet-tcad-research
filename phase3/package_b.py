"""Package B driver (Research_05 v2.1 B-1..B-5), reconstruction run 2026-10-30.

Usage: python phase3/package_b.py --out results/phase3/B_rebuild_2026-10-30/<env_tag>

The 2026-10-29 Package B report is used as an acceptance reference only; every
number written here is recomputed. Production files are read, never written.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from phase3 import m0, runlog
from semiconductor_sim.stage3_mosfet import reference as prod

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF_CSV = os.path.join(ROOT, 'semiconductor_sim/results/phase1/independent_reference.csv')
BG_CSV = os.path.join(ROOT, 'semiconductor_sim/results/phase1/neutral_background.csv')
P1_JSON = os.path.join(ROOT, 'semiconductor_sim/results/phase1/phase1_results.json')
GOLDEN = os.path.join(ROOT, 'semiconductor_sim/results/golden/stage3.json')
GOLDEN_SHA = '17ca32745b69c34355113347eb0f51f95772cf8a76f76b1decdd85925ef67828'

# Acceptance reference values from claude/Phase3_Package_B_Report.md (2026-10-29)
PRIOR = dict(B1_current=0.0, B2=0.0, B3_max_U_rel=1.70e-8, B5=0.0,
             vt_cc_direct=1.0428740, vt_pchip_minus_direct_mV=-0.114,
             vt_phase1=1.0432740, m0_minus_phase1_mV=-0.400, linearity_max=0.0191,
             rii_window_vg=0.50)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        h.update(fh.read())
    return h.hexdigest()


def read_csv(path):
    with open(path, encoding='utf-8-sig') as fh:
        return list(csv.DictReader(fh))


def write_csv(path, rows, header_comment):
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        fh.write(f'# {header_comment}\n')
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    out = os.path.abspath(a.out); os.makedirs(out, exist_ok=False)
    P = m0.M0Params.f0(); TOL = m0.Tolerances(); TIGHT = TOL.tightened(100.)
    rec = runlog.header(dict(params=P.to_dict(), tol=TOL.__dict__, tight=TIGHT.__dict__,
                             prior_report=PRIOR))
    rec['inputs'] = {p: sha256(p) for p in (REF_CSV, BG_CSV, P1_JSON, GOLDEN)}
    assert rec['inputs'][GOLDEN] == GOLDEN_SHA, 'golden changed before run'
    t0 = time.time(); res = {}

    # ---------------- B-1 / B-2 / B-3 on the 42 stored points
    rows = []
    for r in read_csv(REF_CSV):
        vg, vd, Is = float(r['vg']), float(r['vd']), float(r['reference_A'])
        pr = prod.reference_current(vg, vd)
        mm = m0.reference_current(vg, vd, P, TOL)
        mt = m0.reference_current(vg, vd, P, TIGHT)
        dtight = abs(mt['ID'] - mm['ID'])
        U = max(dtight, mm['quadrature_abs_estimate'])
        rows.append(dict(VG_V=vg, VD_V=vd, I_stored_A=Is, I_prod_A=pr['ID'], I_m0_A=mm['ID'],
                         I_m0_tight_A=mt['ID'],
                         e_rep_prod=abs(pr['ID'] - Is) / abs(Is),
                         e_rep_m0=abs(mm['ID'] - Is) / abs(Is),
                         e_b2=abs(mm['ID'] - pr['ID']) / abs(pr['ID']),
                         tight_minus_default_rel=dtight / abs(mm['ID']),
                         quad_est_rel=mm['quadrature_abs_estimate'] / abs(mm['ID']),
                         U_GCA_rel=U / abs(mm['ID']),
                         U_dominant='tight_diff' if dtight >= mm['quadrature_abs_estimate'] else 'quad_estimate'))
    write_csv(os.path.join(out, 'b1_b3_42points.csv'), rows,
              'currents in A for W=1um (= A/um); relative quantities dimensionless; U_GCA is an estimate')
    res['B1'] = dict(n=len(rows), max_e_rep_prod=max(r['e_rep_prod'] for r in rows),
                     max_e_rep_m0=max(r['e_rep_m0'] for r in rows))
    res['B2'] = dict(max_rel=max(r['e_b2'] for r in rows))
    res['B3'] = dict(n=len(rows), max_U_GCA_rel=max(r['U_GCA_rel'] for r in rows),
                     max_tight_minus_default_rel=max(r['tight_minus_default_rel'] for r in rows),
                     dominant=sorted({r['U_dominant'] for r in rows}))

    # ---------------- B-4 derived quantities (F0, L = 1 um)
    vgs = np.round(np.arange(1, 41) * .05, 10)
    vds_mV = [0, .25, .5, 1, 2, 4, 8, 50]
    grid = []
    for vg in vgs:
        s0 = m0.reference_slice(vg, 0., P, TOL)
        G0 = m0.conductance0(vg, P, TOL)
        for vdm in vds_mV:
            vd = vdm * 1e-3
            if vd == 0:
                I = Ie = 0.; fbg = s0['ns_background'] / s0['ns']; ratio = 1.
            else:
                c = m0.reference_current(vg, vd, P, TOL)
                I = c['ID'] / m0.w_um(P); Ie = c['ID_excess'] / m0.w_um(P)
                fbg = c['background_fraction']; ratio = I / (vd * G0)
            grid.append(dict(VG_V=vg, VD_mV=vdm, I_A_per_um=I, I_excess_A_per_um=Ie, f_bg=fbg,
                             Ns_phin0_cm2=s0['ns'] * 1e-4, Ns_inv_phin0_cm2=s0['ns_excess'] * 1e-4,
                             phi_s_bulkref_V=s0['phi_s'], psi_s_intrinsicref_V=s0['phi_s'] + P.psi_b0,
                             G0_S_per_um=G0, I_over_VD_G0=ratio))
    write_csv(os.path.join(out, 'm0_F0.csv'), grid,
              'F0 (NA=1e16 cm-3, tox=20 nm, xmax=1.2 um, L=W=1 um); VD=0 rows: I=0, f_bg from slice; '
              'phi_s bulk-referenced, psi_s = phi_s + psi_b0 intrinsic-referenced')

    ith = m0.ith_per_um(P)
    vt_d, st_d = m0.vt_cc_direct(.05, P, tol=TOL)
    i50 = [g['I_A_per_um'] for g in grid if g['VD_mV'] == 50]
    vt_t, st_t = m0.vt_cc_tabulated(vgs, i50, ith)
    fine = np.round(np.arange(math.floor(vt_d * 100) / 100 - .05, vt_d + .06, .01), 10)
    vt_f, st_f = m0.vt_cc_tabulated(fine, [m0.current_per_um(v, .05, P, TOL) for v in fine], ith)
    vt_p1 = json.load(open(P1_JSON))['threshold']['vt_cc']
    lin = [abs(g['I_over_VD_G0'] - 1) for g in grid if 0 < g['VD_mV'] <= 1]
    res['B4'] = dict(ith_A_per_um=ith, vt_cc_direct_V=vt_d, status=st_d,
                     vt_pchip_005_V=vt_t, pchip_005_minus_direct_mV=(vt_t - vt_d) * 1e3,
                     vt_pchip_001_V=vt_f, pchip_001_minus_direct_mV=(vt_f - vt_d) * 1e3,
                     vt_phase1_table_V=vt_p1, m0_minus_phase1_mV=(vt_d - vt_p1) * 1e3,
                     linearity_max_VD_le_1mV=max(lin), linearity_ref_1mV_over_2VT=1e-3 / (2 * P.vt),
                     psi_b0_V=P.psi_b0)

    # ---------------- B-5 background shares
    b5 = []
    for r in read_csv(BG_CSV):
        vg, vd = float(r['vg']), float(r['vd'])
        f = m0.reference_current(vg, vd, P, TOL)['background_fraction']
        fs = float(r['neutral_background_fraction'])
        b5.append(dict(VG_V=vg, VD_V=vd, f_bg_stored=fs, f_bg_m0=f, rel_diff=abs(f - fs) / fs))
    write_csv(os.path.join(out, 'b5_background.csv'), b5, 'dimensionless background shares')
    spec = {(.05, .05): .962459, (.25, .05): .265584, (.5, .05): .000575826}
    spec_dev = {f'{k}': abs(next(x['f_bg_m0'] for x in b5 if (x['VG_V'], x['VD_V']) == k) - v) / v
                for k, v in spec.items()}
    window = [g['VG_V'] for g in grid if g['VD_mV'] == 50 and g['f_bg'] < 1e-3]
    res['B5'] = dict(n=len(b5), max_rel=max(x['rel_diff'] for x in b5),
                     spec_values_rel_dev_6digit=spec_dev, rii_window_min_VG_V=min(window))

    # ---------------- figure
    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(10, 4))
        for vdm in (1, 50):
            ax[0].semilogy(vgs, [g['I_A_per_um'] for g in grid if g['VD_mV'] == vdm], label=f'VD={vdm} mV')
        ax[0].axhline(ith, ls='--', c='k', lw=.8); ax[0].axvline(vt_d, ls=':', c='k', lw=.8)
        ax[0].set(xlabel='V_G (V)', ylabel="I' (A/um)", title='M0 transfer, F0, L=1 um'); ax[0].legend()
        ax[1].semilogy(vgs, [g['f_bg'] for g in grid if g['VD_mV'] == 50])
        ax[1].axhline(1e-3, ls='--', c='k', lw=.8)
        ax[1].set(xlabel='V_G (V)', ylabel='f_bg', title='Background share, VD=50 mV')
        fig.tight_layout(); fig.savefig(os.path.join(out, 'B_m0_F0.png'), dpi=150)
        fig.savefig(os.path.join(out, 'B_m0_F0.svg'))
    except Exception as exc:
        res['figure_error'] = repr(exc)

    # ---------------- acceptance
    pinned = rec['versions']['numpy'] == '2.3.5' and rec['versions']['scipy'] == '1.17.0'
    acc = dict(
        B1=dict(criterion='<=1e-12' if pinned else '<=1e-9',
                pass_=res['B1']['max_e_rep_prod'] <= (1e-12 if pinned else 1e-9)
                and res['B1']['max_e_rep_m0'] <= (1e-12 if pinned else 1e-9)),
        B2=dict(criterion='<=1e-13', pass_=res['B2']['max_rel'] <= 1e-13),
        B3=dict(criterion='U_GCA,rel<=1e-6 (estimate)', pass_=res['B3']['max_U_GCA_rel'] <= 1e-6),
        # The R05 B-5 values are quoted to 6 significant digits, so they are checked
        # to their rounding half-unit; the 1e-6 criterion uses the full-precision CSV.
        B5=dict(criterion='<=1e-6 vs neutral_background.csv; spec values within rounding',
                pass_=res['B5']['max_rel'] <= 1e-6
                and all(v <= 5e-6 for v in spec_dev.values())),
    )
    rec.update(environment_class='pinned' if pinned else 'current', results=res, acceptance=acc,
               end_utc=runlog.now(), wall_s=time.time() - t0,
               golden_sha_after=sha256(GOLDEN))
    rec['golden_unchanged'] = rec['golden_sha_after'] == GOLDEN_SHA
    runlog.write_json(os.path.join(out, 'package_b_run.json'), rec)
    print(json.dumps(dict(acceptance=acc, results=res, wall_s=rec['wall_s'],
                          golden_unchanged=rec['golden_unchanged']), indent=1, default=float))
    return 0 if all(v['pass_'] for v in acc.values()) and rec['golden_unchanged'] else 2


if __name__ == '__main__':
    sys.exit(main())
