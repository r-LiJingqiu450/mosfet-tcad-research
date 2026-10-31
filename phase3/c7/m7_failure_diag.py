"""C-7 post-hoc failure diagnosis of M7 (NOT a validation): replace M7's quasi-2-D surface-potential profile by the
2-D surface potential of each test state and keep everything else (Ward charge with alpha_0l). If this 'hybrid'
reproduces the 2-D dL, M7's failure is located in its profile; if not, in the charge relation."""
import glob, json, math, sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
from phase3.c7.model_m7 import Device, VT

rows = []
for p in sorted(glob.glob('results/phase3/C7/m7/e7_W_*.json')) + sorted(glob.glob('results/phase3/C7/realism/e7_W_*dpsi0.json')):
    r = json.load(open(p)); a = r['args']
    d = Device(na=a['na'], tox_nm=a['tox'])
    s = r['result']['surface']; x = np.array(s['x_cm']); ps = np.array(s['psi_V'])
    L = x[-1]; h = x <= L / 2 + 1e-15; xm, pm = x[h], ps[h]
    pc = ps[np.argmin(np.abs(x - L / 2))]
    qinf = d.qward(pc, d.ainf)
    q = np.array([d.qward(p_, d.a0(t / d.epsL)) for p_, t in zip(pm, xm)])
    hyb = float(np.trapezoid(1 - qinf / q, xm)) * 1e7
    m7 = d.dL(a['vg'])
    _, pinf, Wd = d.far_field(a['vg'])
    rows.append(dict(tox=a['tox'], NA=a['na'], VG=a['vg'], dL_2D=r['result']['const_mu']['dL_net_nm'],
                     dL_slice=r['result']['const_mu']['dL_slice_nm'], hybrid=hyb, M7=m7['dL_nm'], l0=m7['l0_nm'],
                     alpha=m7['alpha_inf'], ratio_2D_over_l0=r['result']['const_mu']['dL_net_nm'] / m7['l0_nm']))
os.makedirs('results/phase3/C7/analysis', exist_ok=True)
json.dump(rows, open('results/phase3/C7/analysis/m7_failure_diag.json', 'w'), indent=1)
for z in rows:
    print(f"tox {z['tox']:4g} NA {z['NA']:.0e} VG {z['VG']:.4f} a {z['alpha']:5.2f} | 2D {z['dL_2D']:6.1f} slice {z['dL_slice']:6.1f} hybrid {z['hybrid']:6.1f} ({z['hybrid']/z['dL_slice']-1:+.2f}) M7 {z['M7']:6.1f} | dL/l0 {z['ratio_2D_over_l0']:.2f}")
