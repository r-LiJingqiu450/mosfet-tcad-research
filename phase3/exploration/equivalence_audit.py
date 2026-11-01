#!/usr/bin/env python3
"""Summarise corrected equivalence tests for zero-bias conductance.

The endpoint conductance reported by endprobe is a raw sum over the drain
contact.  At low current that sum is cancellation-prone.  This audit treats
(1) the source-side sum and (2) the positive edge-energy identity as the stable
network quantity, then compares it with the independent finite-volume result.
No production solver is changed.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math, subprocess
from datetime import datetime, timezone
from pathlib import Path


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''): h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[2])
    ap.add_argument('--input-dir', type=Path, default=None)
    ap.add_argument('--indep-dir', type=Path, default=None)
    ap.add_argument('--out', type=Path, default=None)
    a = ap.parse_args()
    root = a.root.resolve()
    inp = (a.input_dir or root / 'results/phase3/exploration/equivalence_energy_final').resolve()
    indep = (a.indep_dir or root / 'results/phase3/exploration/equivalence_na1e18').resolve()
    out = (a.out or root / 'results/phase3/exploration/equivalence_audit').resolve()
    out.mkdir(parents=True, exist_ok=True)
    files = sorted(inp.glob('e7_*.json'))
    if len(files) != 3:
        raise RuntimeError(f'expected 3 endprobe JSON files in {inp}, found {len(files)}')
    rows = []
    indep_map = {}
    indep_files = sorted(set(indep.glob('fd2d_VG*.json')) | set((root / 'results/phase3/C11/indep').glob('fd2d_VG*.json')))
    for p in indep_files:
        d = json.loads(p.read_text()); indep_map[round(float(d['args']['vg']), 6)] = d
    for p in files:
        d = json.loads(p.read_text())
        vg = float(d['args']['vg']); c = d['result']['const_mu']; e = d['energy_identity']['constant_mobility']; dv = d['devsim']
        gca = float(d['result']['G_GCA_1COL'])
        gsrc = -float(c['src_current']); gdr = float(c['G_net']); gedge = float(e['edge_energy_G'])
        dl_edge = 0.5 * 1e-4 * (1.0 - gca / gedge) * 1e7
        ref = indep_map.get(round(vg, 6))
        dl_ref = None if ref is None else float(ref['dL_net_nm'])
        rows.append({
            'vg_V': vg, 'G_GCA_S_per_cm': gca, 'G_source_S_per_cm': gsrc,
            'G_edge_energy_S_per_cm': gedge, 'G_drain_raw_S_per_cm': gdr,
            'G_E_equilibrated_S_per_cm': dv['G_E'], 'G_M_UMFPACK_S_per_cm': dv['G_M'],
            'rel_edge_minus_source': gedge / gsrc - 1.0,
            'rel_drain_minus_source': gdr / gsrc - 1.0,
            'rel_G_E_minus_source': dv['G_E'] / gsrc - 1.0,
            'rel_G_M_minus_source': dv['G_M'] / gsrc - 1.0,
            'dL_edge_nm': dl_edge, 'dL_drain_raw_nm': c['dL_net_nm'],
            'dL_independent_FV_nm': dl_ref,
            'rel_dL_edge_vs_FV': None if dl_ref is None else dl_edge / dl_ref - 1.0,
            'network_relative_residual': c['net_rel_residual'],
            'free_node_residual_max': e['free_residual_max'],
            'kappa_G': dv['kappa_G'], 'E_dwr': dv['E_dwr'], 'E_round': dv['E_round'],
            'n_unknowns': dv['n_unknowns'],
        })
    fields = list(rows[0])
    with (out / 'equivalence_rows.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)
    summary = {
        'status': 'CORRECTED_EQUIVALENCE_AUDIT',
        'generated_utc': datetime.now(timezone.utc).isoformat(),
        'git_head': subprocess.check_output(['git','rev-parse','HEAD'], cwd=root, text=True).strip(),
        'inputs_sha256': {str(p.relative_to(root)): sha256(p) for p in [*files, *indep_files]},
        'rows': rows,
        'max_abs_rel_edge_source': max(abs(r['rel_edge_minus_source']) for r in rows),
        'max_abs_rel_dL_edge_vs_FV': max(abs(r['rel_dL_edge_vs_FV']) for r in rows if r['rel_dL_edge_vs_FV'] is not None),
        'max_abs_rel_drain_source': max(abs(r['rel_drain_minus_source']) for r in rows),
        'max_abs_rel_G_E_source': max(abs(r['rel_G_E_minus_source']) for r in rows),
        'max_abs_rel_G_M_source': max(abs(r['rel_G_M_minus_source']) for r in rows),
        'interpretation': {
            'mathematical': 'For a symmetric conductance Laplacian with exact interior solve, edge dissipation phi^T L phi equals the endpoint current. The equality is checked numerically here using the positive edge sum.',
            'numerical': 'The raw drain endpoint sum is cancellation-prone at low current; source-side and edge-energy values agree while the drain sum does not.',
            'physical': 'The equilibrium network and the independent finite-volume solver represent the same zero-bias constant-mobility linear-response quantity for this frozen geometry and boundary condition.',
            'novelty': 'The audit protocol combines standard identities and established conditioning diagnostics; current evidence does not establish a new publishable algorithm.',
        },
    }
    (out / 'equivalence_summary.json').write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k: summary[k] for k in ('status','git_head','max_abs_rel_edge_source','max_abs_rel_dL_edge_vs_FV','max_abs_rel_drain_source','max_abs_rel_G_E_source','max_abs_rel_G_M_source')}, indent=2))

if __name__ == '__main__':
    main()
