"""C-10 design numbers (computed BEFORE any C-10 2-D run): gate voltages for target normalised far-field band
bending psi_hat = (psi_inf - psi_b0)/V_T, similarity variables, eigen-lengths l0, and the job lists.

Similarity map (C-10 note, Sec. 2): in the electron-free limit the zero-bias end correction obeys
dL/L_D = F(psi_hat, psi_hat_bi, t_hat, h_hat, ...), t_hat = (eps_si/eps_ox) t_ox/L_D, h_hat = h_w/L_D,
psi_hat_bi = (psi_bi - psi_b0)/V_T.  Device 1 = (1e16 cm^-3, 20 nm, h_w = 0.5 um) and device 2 = (1e18, 2 nm,
h_w = 0.5 um) share t_hat; device 2' = (1e18, 2 nm, h_w = 0.05 um, psi_bi shifted by V_T ln(1e16/1e18)) is the
similarity image of device 1 (h_hat and psi_hat_bi matched as well).  With the decks' V_FB = 0 convention the
normalised gate drive is V_G/V_T, so device 1 and device 2' are compared at the SAME V_G.
Writes results/phase3/C10/design.json and results/phase3/C10/jobs_*.json.
"""
from __future__ import annotations

import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from scipy.optimize import brentq  # noqa: E402

from phase3.c7 import endprobe as E7  # noqa: E402  (mesh maps only; no DEVSIM solve)
from phase3.c7.model_m7 import Device, VT, Q, ESI, EOX  # noqa: E402

OUT = 'results/phase3/C10'
SHIFT = VT * math.log(1e16 / 1e18)          # psi_bi shift that matches psi_hat_bi of device 2' to device 1


def vg_for(dev, psi_hat):
    f = lambda vg: (dev.far_field(vg)[1] - dev.psib) / VT - psi_hat
    return brentq(f, 1e-3, 6.0, xtol=1e-12)


def info(na, tox, vg, hw_um=0.5, shift=0.0):
    dev = Device(na, tox)
    psig, pinf, Wd = dev.far_field(vg)
    LD = math.sqrt(ESI * VT / (Q * na))
    return dict(na=na, tox=tox, vg=vg, hw_um=hw_um, shift=shift, psi_hat=(pinf - dev.psib) / VT,
                psibi_hat=(dev.psibi + shift - dev.psib) / VT, t_hat=(ESI / EOX) * tox * 1e-7 / LD, h_hat=hw_um * 1e-4 / LD,
                LD_nm=LD * 1e7, Wd_nm=Wd * 1e7, l0_nm=dev.l0(Wd) * 1e7,
                alpha=dev.ainf * math.sqrt(dev.ll) * dev.lam ** ((1 - pinf / (VT * dev.ll)) / 2))


def snapped_hw(na, hw_um, level=0):
    ys = E7.depth_map(na).nodes(level)
    j = min(range(len(ys)), key=lambda i: abs(ys[i] - hw_um * 1e-4))
    return float(ys[j]) * 1e4


def job(deck, L, vg, tox, na, hw=0.5, shift=0.0, level=0, out=OUT + '/runs'):
    args = ['--deck', deck, '--L', f'{L:g}', '--level', str(level), '--vg', f'{vg:.6f}', '--tox', f'{tox:g}', '--na', f'{na:.0e}',
            '--hw', f'{hw:g}', '--psibi-shift', f'{shift:.6f}']
    name = (f"e7_{deck}_L{L:g}_lev{level}_VG{float(f'{vg:.6f}'):g}_tox{tox:g}_NA{na:.0e}_hw{hw:g}_ov0_sg0_dpsi{float(f'{shift:.6f}'):g}.json")
    return dict(script='phase3/c7/endprobe.py', args=args, expect=f'{out}/{name}')


def main():
    os.makedirs(OUT, exist_ok=True)
    d1, d2 = Device(1e16, 20.0), Device(1e18, 2.0)
    vg_a150_d1 = 0.7835                                    # existing C-8 state: device 1, alpha = 150, dL = 416.30 nm
    psi_a150_d1 = (d1.far_field(vg_a150_d1)[1] - d1.psib) / VT
    psi_deep = 16.0
    vg_deep_d1 = vg_for(d1, psi_deep)
    sweep_hats = [27.0, 24.0, 18.0, 15.0, 12.0]
    # h_w is snapped to an existing depth node by the deck; the shallow device-2' junction is placed on the level-0 node
    # near 50 nm (47.27 nm) and device 1 of the similarity pair on its node nearest 10x that depth, so h_hat matches.
    hw_s = round(snapped_hw(1e18, 0.05), 5)
    hw_1s = round(snapped_hw(1e16, 10 * hw_s), 5)
    design = dict(shift_V=SHIFT, psi_hat_device1_alpha150=psi_a150_d1, vg_deep_device1=vg_deep_d1, hw_shallow_um=hw_s,
                  hw_device1_similarity_um=hw_1s, h_hat_device2prime=hw_s * 1e-4 / math.sqrt(ESI * VT / (Q * 1e18)),
                  h_hat_device1_pair=hw_1s * 1e-4 / math.sqrt(ESI * VT / (Q * 1e16)), hw_device1_snapped_um=snapped_hw(1e16, 0.5))
    jobs_sim, jobs_sweep, jobs_tox, jobs_mesh = [], [], [], []
    # similarity block at V_G = 0.7835 (matched to device 1, alpha = 150) and at the deeper matched bias
    blocks = []
    for vg in (vg_a150_d1, vg_deep_d1):
        for hw, sh, tag in ((0.5, 0.0, 'A: device 2 as is'), (0.5, SHIFT, 'B: psi_bi matched'), (hw_s, 0.0, 'C: h_w matched'),
                            (hw_s, SHIFT, "D: device 2' (full similarity image)")):
            if vg == vg_deep_d1 and tag[0] in 'BC':
                continue
            jobs_sim.append(job('W', 1.0, vg, 2.0, 1e18, hw, sh)); blocks.append(dict(tag=tag, **info(1e18, 2.0, vg, hw, sh)))
    for vg in (vg_a150_d1, vg_deep_d1):
        jobs_sim.append(job('W', 2.0, vg, 20.0, 1e16, hw_1s))
        blocks.append(dict(tag='device 1, similarity pair (h_w matched to 10x device 2\')', **info(1e16, 20.0, vg, hw_1s)))
    design['similarity_states'] = blocks
    # bias-range sweeps (device 2: L = 1 um; device 1: L = 2 um), deck W, h_w = 0.5 um
    sw = []
    for dev, na, tox, L in ((d2, 1e18, 2.0, 1.0), (d1, 1e16, 20.0, 2.0)):
        for ph in sweep_hats:
            if na == 1e16 and ph > 21.0:
                continue
            vg = vg_for(dev, ph)
            jobs_sweep.append(job('W', L, vg, tox, na)); sw.append(info(na, tox, vg))
    design['sweep_states'] = sw
    # oxide-thickness similarity breaking (device 2 geometry), at psi_hat of the existing alpha = 150 state and at 21
    psi_a150_d2 = (d2.far_field(1.08)[1] - d2.psib) / VT
    ox = []
    for tox in (6.0, 1.0):
        dv = Device(1e18, tox)
        for ph in (psi_a150_d2, psi_a150_d1):
            vg = vg_for(dv, ph)
            jobs_tox.append(job('W', 1.0, vg, tox, 1e18)); ox.append(info(1e18, tox, vg))
    design['tox_states'] = ox
    design['l0_ratio_pred'] = {f'{s["tox"]:g}nm_psi{s["psi_hat"]:.2f}': s['l0_nm'] / info(1e18, 2.0, vg_for(d2, s['psi_hat']))['l0_nm'] for s in ox}
    # mesh check of the key similarity comparison: device 1 at alpha = 150 and device 2' at the same V_G, level 1
    jobs_mesh.append(job('W', 2.0, vg_a150_d1, 20.0, 1e16, level=1))
    jobs_mesh.append(job('W', 1.0, vg_a150_d1, 2.0, 1e18, 0.5, 0.0, level=1))
    for n, j in (('sim', jobs_sim), ('sweep', jobs_sweep), ('tox', jobs_tox), ('mesh', jobs_mesh)):
        json.dump(j, open(f'{OUT}/jobs_{n}.json', 'w'), indent=1)
    json.dump(design, open(f'{OUT}/design.json', 'w'), indent=1)
    print(json.dumps({k: v for k, v in design.items() if not isinstance(v, list)}, indent=1))
    for k in ('similarity_states', 'sweep_states', 'tox_states'):
        print(k)
        for s in design[k]:
            print('  ', {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in s.items()})
    print('l0 ratio predictions', {k: round(v, 4) for k, v in design['l0_ratio_pred'].items()})
    print('jobs', {n: len(j) for n, j in (('sim', jobs_sim), ('sweep', jobs_sweep), ('tox', jobs_tox), ('mesh', jobs_mesh))})


if __name__ == '__main__':
    main()
