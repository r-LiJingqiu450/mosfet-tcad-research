"""E8-1 evaluation (pre-registered rule): H-sat vs H-log for the weak-inversion branch of kappa = dL/l0."""
import glob, json, math, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from phase3.c7.model_m7 import Device, VT
P = json.load(open('results/phase3/C8/subth/preregistration_data.json'))
res = []
for s in P['states']:
    f = glob.glob(f"results/phase3/C8/subth/e7_W_L{s['L']:g}_lev0_VG{s['VG']:g}_tox{s['tox']:g}_NA{s['NA']:.0e}_*.json")[0]
    r = json.load(open(f)); c = r['result']['const_mu']
    res.append(dict(**s, dL_net=c['dL_net_nm'], dL_slice=c['dL_slice_nm'], kappa=c['dL_net_nm'] / s['l0_nm'],
                    bounds=c['bounds_ok'], Ns=r['result']['Ns_cm2'], newton_final=r['equilibrium']['final']['last'][-1]))
out = dict(states=res, verdict={})
for na in (1e16, 1e18):
    e = sorted([z for z in res if z['NA'] == na], key=lambda z: z['alpha'])
    out['verdict'][f'{na:.0e}'] = dict(dkappa=e[-1]['kappa'] - e[0]['kappa'], kappas=[z['kappa'] for z in e])
dk = [v['dkappa'] for v in out['verdict'].values()]
out['verdict']['E8-1'] = 'H-sat SUPPORTED' if max(dk) <= 0.9 else ('H-log SUPPORTED' if min(dk) >= 1.8 else 'INCONCLUSIVE')
json.dump(out, open('results/phase3/C8/subth/e8_1_result.json', 'w'), indent=1)
for z in res:
    print(f"NA {z['NA']:.0e} alpha {z['alpha']:6.1f} VG {z['VG']} Ns {z['Ns']:.2e} dL {z['dL_net']:.2f} kappa {z['kappa']:.3f} (H-log {z['kappa_Hlog']:.2f}) bounds {z['bounds']}")
print(json.dumps(out['verdict'], indent=1))
