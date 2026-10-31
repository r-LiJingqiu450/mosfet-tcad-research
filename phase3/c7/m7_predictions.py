"""Pre-registered predictions of model M7 (and of Ward's original Sec. 6.4) for UNSEEN states.
Also records M7's in-sample (training) comparison with the C-6 2-D data, labelled as such.
Output: results/phase3/C7/m7/m7_predictions.json and jobs_m7_test.json (2-D test jobs, deck W, L = 2 um)."""
from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from phase3.c7.model_m7 import Device, dL_ward, ns_1d, vg_for_ns

OUT = 'results/phase3/C7/m7'
TARGETS = {'threshold': 6.24e9, 'moderate': 1.607e11, 'strong': 8.967e11}   # C-6 inversion levels (1.0/1.25/2.0 V)
STRUCT = {'S1_tox10': dict(tox=10.0, na=1e16), 'S2_tox40': dict(tox=40.0, na=1e16), 'S3_NA4e16': dict(tox=20.0, na=4e16)}


def main():
    os.makedirs(OUT, exist_ok=True)
    rec = dict(note='M7 frozen at commit of this file; predictions made before any 2-D run of these states.',
               training_C6=[], predictions=[])
    d0 = Device()
    c6 = {1.0: 369.873, 1.25: 153.687, 1.5: 92.878, 2.0: 53.385, 3.0: 29.909}   # C-6 W1 dL_net, L = 4 um, level 0
    for vg, v in c6.items():
        m = d0.dL(vg)
        rec['training_C6'].append(dict(VG=vg, dL_2D_net=v, dL_M7=m['dL_nm'], ratio=v / m['dL_nm'], dL_Ward=dL_ward(d0, vg),
                                       alpha_inf=m['alpha_inf'], l0_nm=m['l0_nm']))
    states = [('S0_original', 20.0, 1e16, 1.1, 'unseen bias')]
    for sname, p in STRUCT.items():
        dv = Device(na=p['na'], tox_nm=p['tox'])
        for lvl, ns in TARGETS.items():
            states.append((sname, p['tox'], p['na'], round(vg_for_ns(dv, ns), 4), lvl))
    jobs = []
    for sname, tox, na, vg, lvl in states:
        dv = Device(na=na, tox_nm=tox)
        m = dv.dL(vg)
        _, pinf, _ = dv.far_field(vg)
        rec['predictions'].append(dict(state=sname, level_label=lvl, tox_nm=tox, N_A=na, VG=vg, Ns_1D=ns_1d(dv, pinf),
                                       dL_M7_nm=m['dL_nm'], dL_Ward_nm=dL_ward(dv, vg), alpha_inf=m['alpha_inf'],
                                       l0_nm=m['l0_nm'], Wd_nm=m['Wd_nm']))
        jobs.append(dict(script='phase3/c7/endprobe.py',
                         args=['--deck', 'W', '--L', '2', '--level', '0', '--vg', f'{vg:g}', '--tox', f'{tox:g}', '--na', f'{na:g}'],
                         expect=f'results/phase3/C7/m7/e7_W_L2_lev0_VG{vg:g}_tox{tox:g}_NA{na:.0e}_hw0.5_ov0_sg0_dpsi0.json'))
    json.dump(rec, open(os.path.join(OUT, 'm7_predictions.json'), 'w'), indent=1)
    json.dump(jobs, open('results/phase3/C7/jobs_m7_test.json', 'w'), indent=1)
    for t in rec['training_C6']:
        print('train', t['VG'], round(t['dL_M7'], 1), t['dL_2D_net'], round(t['ratio'], 3), 'Ward', round(t['dL_Ward'], 1))
    for p in rec['predictions']:
        print('pred', p['state'], p['level_label'], p['VG'], f"{p['Ns_1D']:.3e}", 'M7', round(p['dL_M7_nm'], 1),
              'Ward', round(p['dL_Ward_nm'], 1), 'alpha', round(p['alpha_inf'], 2), 'l0', round(p['l0_nm'], 1))


if __name__ == '__main__':
    main()
