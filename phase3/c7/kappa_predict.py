"""Pre-registered out-of-sample test of the kappa(alpha_inf) scaling law (C-7 Addendum B).
kappa curve = piecewise-linear interpolation in ln(alpha) of the pooled points in kappa_table.csv
(all deck-W abrupt states of C-6/C-7 incl. Addendum A; equal-alpha points averaged). Predictions:
dL = l0(V_G) * kappa(alpha_inf(V_G)) for two never-simulated devices at four alpha targets."""
import csv, json, math, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
from scipy.optimize import brentq
from phase3.c7.model_m7 import Device, VT

T = list(csv.DictReader(open('results/phase3/C7/analysis/kappa_table.csv')))
pts = {}
for r in T:
    a = round(float(r['alpha_inf']), 2); pts.setdefault(a, []).append(float(r['kappa']))
A = np.array(sorted(pts)); K = np.array([np.mean(pts[a]) for a in A])
def kappa(alpha):
    return float(np.interp(math.log(alpha), np.log(A), K))
out = dict(curve=dict(alpha=A.tolist(), kappa=K.tolist(), n_points=len(T)), rule='SUPPORTED if every |dL_2D/dL_pred - 1| <= 0.12; REFUTED if any > 0.20; else PARTIAL', states=[])
jobs = []
for na, tox in ((3e17, 5.0), (1e17, 10.0)):
    d = Device(na=na, tox_nm=tox)
    def alpha(vg):
        _, p, _ = d.far_field(vg); return d.ainf * math.sqrt(d.ll) * d.lam ** ((1 - p / (VT * d.ll)) / 2)
    for at in (3.0, 1.2, 0.6, 0.2):
        vg = round(brentq(lambda v: alpha(v) - at, 0.2, 6.0, xtol=1e-7), 4)
        _, p, Wd = d.far_field(vg); l0 = d.l0(Wd) * 1e7; al = alpha(vg)
        out['states'].append(dict(NA=na, tox_nm=tox, VG=vg, alpha_inf=al, l0_nm=l0, kappa_pred=kappa(al), dL_pred_nm=l0 * kappa(al)))
        jobs.append(dict(script='phase3/c7/endprobe.py', args=['--deck', 'W', '--L', '1', '--level', '0', '--vg', f'{vg:g}', '--tox', f'{tox:g}', '--na', f'{na:g}'],
                         expect=f'results/phase3/C7/kappa_test/e7_W_L1_lev0_VG{vg:g}_tox{tox:g}_NA{na:.0e}_hw0.5_ov0_sg0_dpsi0.json'))
json.dump(out, open('results/phase3/C7/kappa_test/predictions.json', 'w'), indent=1)
json.dump(jobs, open('results/phase3/C7/jobs_kappa_test.json', 'w'), indent=1)
for s in out['states']: print(s)
