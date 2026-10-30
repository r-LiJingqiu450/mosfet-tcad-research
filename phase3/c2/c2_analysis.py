"""C-2 analysis: T0/T2 against M0, numerical uncertainty, P/F equivalence, figures.

Inputs: results/phase3/C2/onecol/{P,F}_L{0..3}_tau1e-10.json, P_L3_tau1e-12.json,
P_L{0..3}_tau1e-10_vtlong.json, P_L3_tau1e-12_vtlong.json.

Tests (both published forms are evaluated; they differ by a factor 2):
  T0  R05 C-2 : |psi_s,1COL,finest - psi_s,M0| <= U_disc(psi_s) + 1e-6 V
      R04 8.1 : ... <= 2 U_disc(psi_s) + 1e-6 V
  T2  R05 C-2 : |N_s,1COL,finest/N_s,M0 - 1| <= U_disc,rel + U_GCA,rel
      R04 8.1 : ... <= 2 (U_disc,rel + U_GCA,rel)
U_disc per R04 7.4 on levels 3 (finest) .. 0, r = sqrt(2), no order assumed.
Verdict per point:
  PASS          criterion met and the 4-level sequence is classified asymptotic
                (or, for T0, |diff| <= 1e-6 V, which does not depend on U);
  FAIL          criterion violated with an asymptotic classification;
  INCONCLUSIVE  the uncertainty estimate is not reliable (non-asymptotic /
                oscillatory) and the decision depends on it.
U_iter = |f(tau) - f(tau/100)| on the finest level (T5: U_iter <= 0.1 U_disc).
U_GCA  = M0 quadrature/root uncertainty: |f(tol) - f(tol/100)| combined with the
         `quad` error estimate (N_s) or the brentq xtol (psi_s). All U are estimates.
"""
from __future__ import annotations

import csv
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np

from phase3 import m0, runlog
from common import constants as C

R = math.sqrt(2.0)
BASE = 'results/phase3/C2/onecol'
OUT = 'results/phase3/C2/t0t2'
PSI_B0 = m0.M0Params.f0().psi_b0


def load(name):
    return json.load(open(os.path.join(BASE, name)))


def key(s):
    return (round(s['VG_V'], 10), round(s['phin_V'], 10))


def index(rec):
    return {key(s): s for s in rec['states']}


def classify(f):
    """R04 7.4. f = [f1 (finest), f2, f3, f4]."""
    f1, f2, f3, f4 = f
    d12, d23, d34 = f1 - f2, f2 - f3, f3 - f4
    out = dict(f=f, R123=None, R234=None, p123=None, p234=None)
    if d23 == 0 or d34 == 0:
        if d12 == 0 and d23 == 0:
            return dict(out, cls='converged_to_roundoff', U_disc=0.0, f_rich=f1)
        return dict(out, cls='oscillatory', U_disc=3 * (max(f[:3]) - min(f[:3])) / 2, f_rich=None)
    R123, R234 = d12 / d23, d23 / d34
    out.update(R123=R123, R234=R234)
    if R123 < 0 or R234 < 0:
        return dict(out, cls='oscillatory', U_disc=3 * (max(f[:3]) - min(f[:3])) / 2, f_rich=None)
    p123 = math.log(abs(d23 / d12)) / math.log(R) if d12 != 0 else float('inf')
    p234 = math.log(abs(d34 / d23)) / math.log(R)
    out.update(p123=p123, p234=p234)
    if 0 < R123 < 1 and 0 < R234 < 1 and 0.5 <= p123 <= 3 and abs(p123 - p234) <= 0.25 * p123:
        U = 1.25 * abs(d12) / (R ** p123 - 1)
        return dict(out, cls='asymptotic', U_disc=U, f_rich=f1 + d12 / (R ** p123 - 1))
    # monotone, non-asymptotic: least-squares f0 + a h^p (p free) and f0 + a1 h + a2 h^2
    h = np.array([1.0, R, R ** 2, R ** 3]); y = np.array(f)
    best = None
    for p in np.linspace(0.3, 4.0, 371):
        A = np.vstack([np.ones(4), h ** p]).T; c, res, *_ = np.linalg.lstsq(A, y, rcond=None)
        rr = float(np.sum((A @ c - y) ** 2))
        if best is None or rr < best[0]:
            best = (rr, c[0], p)
    A2 = np.vstack([np.ones(4), h, h ** 2]).T; c2, *_ = np.linalg.lstsq(A2, y, rcond=None)
    dM = max(abs(np.max(y) - np.min(y)), 0.0)
    U = 3 * max(abs(f1 - best[1]), abs(f1 - c2[0]), dM)
    return dict(out, cls='monotone_nonasymptotic', U_disc=U, f_rich=None, p_fit=best[2])


def m0_point(vg, phin):
    P = m0.M0Params.f0()
    if vg == 0.0:
        # flat band, exactly: V_FB + phi_s - Q_s(phi_s)/C_ox = 0 has the root phi_s = 0
        # (Q_s(0) = 0); M0's reference_slice refuses V_G < 0.05 V by design.
        return dict(psi_s=PSI_B0, U_psi=0.0, Ns_cm2=None, U_Ns=None, note='analytic flat band')
    a = m0.reference_slice(vg, phin, P); b = m0.reference_slice(vg, phin, P, m0.DEFAULT_TOL.tightened(100.))
    U_psi = abs(a['phi_s'] - b['phi_s']) + 2e-14
    U_ns = max(abs(a['ns'] - b['ns']), a['quadrature_abs_estimate'])
    return dict(psi_s=a['phi_s'] + PSI_B0, U_psi=U_psi, Ns_cm2=a['ns'] * 1e-4, U_Ns=U_ns * 1e-4,
                bulk_shift_V=a['bulk_shift'], note='m0.reference_slice (phi_s + psi_b0)')


def main():
    os.makedirs(OUT, exist_ok=True)
    rec = runlog.header(dict(analysis='C-2 T0/T2', base=BASE))
    P = [index(load(f'P_L{k}_tau1e-10.json')) for k in range(4)]
    Pv = [index(load(f'P_L{k}_tau1e-10_vtlong.json')) for k in range(4)]
    F = [index(load(f'F_L{k}_tau1e-10.json')) for k in range(4)]
    Pt = index(load('P_L3_tau1e-12.json')); Ptv = index(load('P_L3_tau1e-12_vtlong.json'))
    for k in range(4):
        P[k].update(Pv[k])
    Pt.update(Ptv)
    vtl = float(C.VT_LONG)

    rows = []
    t0_pts = [(0.0, 0.0), (round(vtl, 10), 0.0), (2.0, 0.0)]
    t2_pts = [(vg, ph) for vg in (0.8, 1.25, 2.0) for ph in (0.0, 0.05, 0.5)]
    for test, pts, qty in (('T0', t0_pts, 'psi_s_V'), ('T2', t2_pts, 'Ns_FV_cm2')):
        for kp in pts:
            f = [P[3][kp][qty], P[2][kp][qty], P[1][kp][qty], P[0][kp][qty]]
            c = classify(f)
            U_iter = abs(P[3][kp][qty] - Pt[kp][qty])
            ref = m0_point(*kp)
            if test == 'T0':
                diff = f[0] - ref['psi_s']; ad = abs(diff)
                lim05 = c['U_disc'] + 1e-6; lim04 = 2 * c['U_disc'] + 1e-6
                if ad <= 1e-6:
                    v05 = v04 = 'PASS'
                elif c['cls'] in ('asymptotic', 'converged_to_roundoff'):
                    v05 = 'PASS' if ad <= lim05 else 'FAIL'; v04 = 'PASS' if ad <= lim04 else 'FAIL'
                else:
                    v05 = v04 = 'INCONCLUSIVE'
                rows.append(dict(test='T0', VG_V=kp[0], phin_V=kp[1], quantity='psi_s [V, intrinsic ref]',
                                 M0=ref['psi_s'], U_GCA=ref['U_psi'], L3=f[0], L2=f[1], L1=f[2], L0=f[3],
                                 diff_signed=diff, diff_abs=ad, rel_diff=None, cls=c['cls'], R123=c['R123'],
                                 R234=c['R234'], p123=c['p123'], p234=c['p234'], U_disc=c['U_disc'],
                                 U_iter=U_iter, T5_ok=(U_iter <= 0.1 * c['U_disc']) if c['U_disc'] else U_iter < 1e-14,
                                 f_richardson_diag=c['f_rich'],
                                 rich_minus_M0=(c['f_rich'] - ref['psi_s']) if c['f_rich'] is not None else None,
                                 limit_R05=lim05, verdict_R05=v05, limit_R04=lim04, verdict_R04=v04, note=ref['note']))
            else:
                diff = f[0] - ref['Ns_cm2']; rel = f[0] / ref['Ns_cm2'] - 1
                Ud = c['U_disc'] / abs(f[0]); Ug = ref['U_Ns'] / ref['Ns_cm2']
                lim05 = Ud + Ug; lim04 = 2 * (Ud + Ug)
                if c['cls'] in ('asymptotic', 'converged_to_roundoff'):
                    v05 = 'PASS' if abs(rel) <= lim05 else 'FAIL'; v04 = 'PASS' if abs(rel) <= lim04 else 'FAIL'
                else:
                    v05 = v04 = 'INCONCLUSIVE'
                rows.append(dict(test='T2', VG_V=kp[0], phin_V=kp[1], quantity='N_s [cm^-2]',
                                 M0=ref['Ns_cm2'], U_GCA=ref['U_Ns'], L3=f[0], L2=f[1], L1=f[2], L0=f[3],
                                 diff_signed=diff, diff_abs=abs(diff), rel_diff=rel, cls=c['cls'], R123=c['R123'],
                                 R234=c['R234'], p123=c['p123'], p234=c['p234'], U_disc=c['U_disc'],
                                 U_disc_rel=Ud, U_GCA_rel=Ug, U_iter=U_iter,
                                 T5_ok=(U_iter <= 0.1 * c['U_disc']) if c['U_disc'] else U_iter < 1e-14 * abs(f[0]),
                                 f_richardson_diag=c['f_rich'],
                                 rich_rel_minus_M0=(c['f_rich'] / ref['Ns_cm2'] - 1) if c['f_rich'] is not None else None,
                                 limit_R05=lim05, verdict_R05=v05, limit_R04=lim04, verdict_R04=v04, note=ref['note']))
    keys = sorted({k for r in rows for k in r})
    with open(os.path.join(OUT, 't0t2.csv'), 'w', newline='') as fh:
        fh.write('# C-2 T0/T2; psi in V (intrinsic-referenced), N_s in cm^-2; U values are estimates\n')
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
    rec['rows'] = rows

    # ---- P vs F equivalence (all common states, every level)
    eq = []
    for k in range(4):
        common = sorted(set(P[k]) & set(F[k]))
        dpsi = max(abs(P[k][c]['psi_s_V'] - F[k][c]['psi_s_V']) for c in common)
        dns = max(abs(P[k][c]['Ns_FV_cm2'] / F[k][c]['Ns_FV_cm2'] - 1) for c in common)
        eq.append(dict(level=k, n_states=len(common), max_abs_dpsi_s_V=dpsi, max_rel_dNs=dns,
                       F_continuity_residual_rel_max=max(F[k][c]['continuity_residual_rel'] for c in common),
                       F_sg_flux_rel_max=max(F[k][c]['sg_flux_rel_to_term_scale'] for c in common),
                       F_n_relation_max_rel=max(F[k][c]['n_relation_max_rel'] for c in common),
                       F_lateral_psi_spread_max_V=max(F[k][c]['lateral_psi_spread_V'] for c in common),
                       P_sg_flux_rel_max=max(P[k][c]['sg_flux_rel_to_term_scale'] for c in common),
                       P_continuity_residual_rel_max=max(P[k][c]['continuity_residual_rel'] for c in common),
                       max_gate_bc_err_V=max(abs(P[k][c]['gate_bc_err_V']) for c in common),
                       max_bottom_bc_err_V=max(abs(P[k][c]['bottom_bc_err_V']) for c in common),
                       gauss_rel_max_Qs_nonzero=max(P[k][c]['gauss_rel'] for c in common if c[0] > 0)))
    rec['equivalence_P_vs_F'] = eq

    # ---- full V_G sweep at phi_n = 0, finest level (diagnostic, no criterion)
    sweep = []
    for vg in np.round(np.arange(0.05, 2.0 + 1e-9, 0.05), 10):
        kp = (vg, 0.0); ref = m0_point(vg, 0.0); s = P[3][kp]
        sweep.append(dict(VG_V=vg, psi_s_L3=s['psi_s_V'], psi_s_M0=ref['psi_s'], dpsi_V=s['psi_s_V'] - ref['psi_s'],
                          Ns_L3=s['Ns_FV_cm2'], Ns_M0=ref['Ns_cm2'], rel_dNs=s['Ns_FV_cm2'] / ref['Ns_cm2'] - 1,
                          Qs_FV_C_cm2=s['Qs_FV_C_cm2'], gauss_rel=s['gauss_rel']))
    with open(os.path.join(OUT, 'sweep_L3_vs_M0.csv'), 'w', newline='') as fh:
        fh.write('# finest level (L3) vs M0, phi_n = 0; psi V, N_s cm^-2, Q_s C/cm^2 (diagnostic only)\n')
        w = csv.DictWriter(fh, fieldnames=list(sweep[0].keys())); w.writeheader(); w.writerows(sweep)
    rec['sweep_summary'] = dict(max_abs_dpsi_V=max(abs(s['dpsi_V']) for s in sweep),
                                max_abs_rel_dNs=max(abs(s['rel_dNs']) for s in sweep))

    # ---- figures
    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        hs = np.array([R ** -3, R ** -2, R ** -1, 1.0])          # relative h, L3..L0
        for r in rows:
            err = np.abs(np.array([r['L3'], r['L2'], r['L1'], r['L0']]) - r['M0'])
            if r['test'] == 'T2':
                err = err / abs(r['M0'])
            a = ax[0] if r['test'] == 'T0' else ax[1]
            if np.all(err > 0):
                a.loglog(hs, err, 'o-', ms=3, lw=1, label=f"VG={r['VG_V']:.3g}, phin={r['phin_V']:.2g}")
        for a, t in zip(ax, ('T0: |psi_s,h - psi_s,M0| (V)', 'T2: |N_s,h/N_s,M0 - 1|')):
            a.loglog(hs, 1e-6 * hs ** 2 / hs[-1] ** 2 * 0.1, 'k:', lw=.8, label='slope 2 (guide)')
            a.set_xlabel('relative mesh size h/h0'); a.set_title(t); a.legend(fontsize=6)
        fig.tight_layout(); fig.savefig(os.path.join(OUT, 'C2_convergence.png'), dpi=150)
        fig.savefig(os.path.join(OUT, 'C2_convergence.svg'))
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        vg = [s['VG_V'] for s in sweep]
        ax[0].semilogy(vg, [abs(s['dpsi_V']) + 1e-18 for s in sweep], '.-'); ax[0].set_title('|psi_s,L3 - psi_s,M0| (V), phi_n = 0')
        ax[1].semilogy(vg, [abs(s['rel_dNs']) + 1e-18 for s in sweep], '.-'); ax[1].set_title('|N_s,L3/N_s,M0 - 1|, phi_n = 0')
        for a in ax: a.set_xlabel('V_G (V)')
        fig.tight_layout(); fig.savefig(os.path.join(OUT, 'C2_sweep_L3_vs_M0.png'), dpi=150)
        fig.savefig(os.path.join(OUT, 'C2_sweep_L3_vs_M0.svg'))
    except Exception as exc:
        rec['figure_error'] = repr(exc)
    rec['end_utc'] = runlog.now()
    runlog.write_json(os.path.join(OUT, 'c2_t0t2.json'), rec)
    for r in rows:
        print('%s VG=%-6.4g phin=%-4.2g cls=%-22s diff=%+.3e rel=%s Udisc=%.2e Uiter=%.1e UGCA=%.1e R05=%s R04=%s p=%s' % (
            r['test'], r['VG_V'], r['phin_V'], r['cls'], r['diff_signed'],
            '%+.2e' % r['rel_diff'] if r['rel_diff'] is not None else '   -   ', r['U_disc'], r['U_iter'],
            r['U_GCA'], r['verdict_R05'], r['verdict_R04'], ('%.2f/%.2f' % (r['p123'], r['p234'])) if r['p123'] else '-'))
    print('EQUIV', json.dumps(eq, default=float))
    print('SWEEP', rec['sweep_summary'])
    return 0


if __name__ == '__main__':
    sys.exit(main())
