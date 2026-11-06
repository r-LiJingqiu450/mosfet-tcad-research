"""C-1(b): is get_matrix_and_rhs(format='csc') the Newton Jacobian, incl. contact rows?

Research_05 v2.1 C-1(b). Test protocol (no damping in any compared update):
  1. Converge the R-I(a)-type deck at V_G = 1.25 V, V_D = 0 (circuit source V1 = 0).
  2. Re-declare all bulk equations with variable_update='default' (x <- x + dx).
  3. Set V1 = delta (state now not converged; Newton step is non-trivial).
  4. At this state X0: export (A_exp, b_exp) = get_matrix_and_rhs('csc').
     Check X0 is unchanged by the export call.
  5. One DEVSIM Newton iteration with a recording solver callback: captures the
     exact (A_nwt, b_nwt) factored by Newton; dx_cb = A_nwt^-1 b_nwt (scipy).
     Record the state X1; dx_state = X1 - X0.
  6. Restore X0; one DEVSIM Newton iteration with DEVSIM's own default direct
     solver (UMFPACK via its stock callback); dx_builtin = X1' - X0.
  7. Compare dx_ext = A_exp^-1 b_exp (scipy SuperLU) with dx_builtin.
  8. Linearization check: central-difference columns of b_exp (b = -F, since
     Newton.hh sets rhssign = -1) vs columns of A_exp at sampled unknowns,
     including contact, interface, oxide, gate and circuit unknowns.
Outputs: results/phase3/C1/c1b/<run>/c1b.json, matrices (.npz), log.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import devsim as ds

from phase3 import runlog
from phase3.c1 import deck as D

REC = {}


def recording_callback(**kw):
    a = kw['action']
    if a == 'init':
        REC['init'] = dict(n=kw['n'], transpose=kw['transpose'])
        return dict(solver_object=object(), matrix_format='csc', status=True, message='')
    if a == 'factor':
        n = REC['init']['n']
        A = sp.csc_matrix((np.array(kw['Ax'], float), np.array(kw['Ai']), np.array(kw['Ap'])),
                          shape=(n, n))
        REC.setdefault('factor', []).append(dict(A=A, is_complex=kw['is_complex']))
        # equation numbering is only valid during a solve: capture it here
        REC['eqnums'] = equation_numbers()
        return dict(status=True, message='')
    if a == 'solve':
        A = REC['factor'][-1]['A']
        b = np.array(kw['b'], float)
        M = A.T.tocsc() if REC['init']['transpose'] else A
        x = spla.splu(M).solve(b)
        REC.setdefault('solve', []).append(dict(b=b, x=x))
        return dict(status=True, message='', x=list(x))
    raise RuntimeError(a)


def equation_numbers():
    out = {}
    for r, var in ((D.SI, 'Potential'), (D.SI, 'Electrons'), (D.OX, 'Potential')):
        out[(r, var)] = np.array(ds.get_equation_numbers(device=D.DEV, region=r, variable=var))
    for node in ds.get_circuit_node_list():
        try:
            out[('circuit', node)] = np.array([ds.get_circuit_equation_number(node=node)])
        except ds.error:
            pass
    return out


def snapshot():
    s = {}
    for r, var in ((D.SI, 'Potential'), (D.SI, 'Electrons'), (D.OX, 'Potential')):
        s[(r, var)] = np.array(ds.get_node_model_values(device=D.DEV, region=r, name=var))
    for node in ds.get_circuit_node_list():
        s[('circuit', node)] = np.array([ds.get_circuit_node_value(solution='dcop', node=node)])
    return s


def restore(s):
    for (r, var), v in s.items():
        if r == 'circuit':
            ds.set_circuit_node_value(solution='dcop', node=var, value=float(v[0]))
        else:
            ds.set_node_values(device=D.DEV, region=r, name=var, values=list(v))


def to_vector(s, eqn, n):
    """Map a state dict into the global unknown vector using equation numbers."""
    v = np.full(n, np.nan)
    for k, rows in eqn.items():
        if k in s:
            valid = rows >= 0
            v[rows[valid]] = s[k][valid]
    return v


def export():
    m = ds.get_matrix_and_rhs(format='csc')
    st = m['static']; n = len(st['rhs'])
    A = sp.csc_matrix((np.array(st['av'], float), np.array(st['ai']), np.array(st['ap'])), shape=(n, n))
    dyn = m['dynamic']
    Ad = sp.csc_matrix((np.array(dyn['av'], float), np.array(dyn['ai']), np.array(dyn['ap'])), shape=(n, n))
    perm = m['permutation']
    return A, np.array(st['rhs'], float), Ad, np.array(dyn['rhs'], float), \
        np.array(perm['row']), np.array(perm['copy']), sorted(m.keys())


def block_rel(dx_a, dx_b, eqn):
    """Per-unknown-block relative difference max|a-b| / max|b|."""
    out = {}
    for k, rows in eqn.items():
        rows = rows[rows >= 0]
        if len(rows) == 0:
            continue
        a, b = dx_a[rows], dx_b[rows]
        den = np.max(np.abs(b))
        out[f'{k[0]}:{k[1]}'] = dict(max_abs_diff=float(np.max(np.abs(a - b))), max_abs_ref=float(den),
                                     rel=float(np.max(np.abs(a - b)) / den) if den > 0 else None)
    return out


def run_state(delta, outdir, tag, n_fd_per_block=6):
    res = dict(delta_V=delta)
    stock_cb = ds.get_parameter(name='solver_callback')
    ds.circuit_alter(name='V1', value=delta)
    X0 = snapshot()
    A, b, Ad, bd, prow, pcopy, keys = export()
    X0b = snapshot()
    res['export_keys'] = keys
    res['export_changes_state_maxabs'] = {f'{k[0]}:{k[1]}': float(np.max(np.abs(X0[k] - X0b[k]))) for k in X0}
    n = A.shape[0]

    # ---- one Newton iteration with recording callback
    REC.clear()
    ds.set_parameter(name='direct_solver', value='custom')
    ds.set_parameter(name='solver_callback', value=recording_callback)
    try:
        info_cb = D.solve(abs_err=1e300, rel_err=1e300, it=1)
    finally:
        ds.set_parameter(name='solver_callback', value=stock_cb)
    X1 = snapshot()
    eqn = REC['eqnums']
    A_n = REC['factor'][0]['A']; b_n = REC['solve'][0]['b']; x_cb = REC['solve'][0]['x']
    res['callback'] = dict(n=REC['init']['n'], transpose=REC['init']['transpose'],
                           n_factor_calls=len(REC['factor']), n_solve_calls=len(REC['solve']),
                           converged_flag=info_cb['converged'])
    # ordering coverage
    allrows = np.concatenate([r[r >= 0] for r in eqn.values()])
    res['ordering'] = dict(n_unknowns=n, n_mapped=int(len(allrows)), n_unique=int(len(np.unique(allrows))),
                           covers_0_to_n_minus_1=bool(np.array_equal(np.sort(np.unique(allrows)), np.arange(n))),
                           blocks={f'{k[0]}:{k[1]}': [int(r.min()), int(r.max()), int(len(r))]
                                   for k, r in eqn.items() if len(r[r >= 0])},
                           unnumbered_nodes={f'{k[0]}:{k[1]}': int(np.sum(r < 0)) for k, r in eqn.items()})
    # matrix / rhs identity: export vs Newton
    dA = (A - A_n)
    res['A_export_vs_newton'] = dict(
        same_pattern=bool(np.array_equal(A.indptr, A_n.indptr) and np.array_equal(A.indices, A_n.indices)),
        max_abs_diff=float(abs(dA).max()) if dA.nnz else 0.0, max_abs_A=float(abs(A).max()), nnz=int(A.nnz))
    res['rhs_export_vs_newton'] = dict(max_abs_diff=float(np.max(np.abs(b - b_n))),
                                       max_abs_rhs=float(np.max(np.abs(b))))
    dx_state_cb = to_vector({k: X1[k] - X0[k] for k in X0}, eqn, n)
    res['dx_state_vs_callback_solution'] = block_rel(dx_state_cb, x_cb, eqn)

    # ---- restore and do one iteration with DEVSIM's stock direct solver
    restore(X0)
    A2, b2, *_ = export()
    res['restore_check'] = dict(A_max_abs_diff=float(abs(A2 - A).max()) if (A2 - A).nnz else 0.0,
                                rhs_max_abs_diff=float(np.max(np.abs(b2 - b))))
    info_bi = D.solve(abs_err=1e300, rel_err=1e300, it=1)
    X1b = snapshot()
    dx_builtin = to_vector({k: X1b[k] - X0[k] for k in X0}, eqn, n)

    # ---- external solve on the exported system
    lu = spla.splu(A.tocsc())
    dx_ext = lu.solve(b)
    res['dx_ext_vs_builtin'] = block_rel(dx_ext, dx_builtin, eqn)
    res['dx_minus_ext_vs_builtin'] = block_rel(-dx_ext, dx_builtin, eqn)
    # componentwise (Oettli-Prager) backward error of each solution
    absA = abs(A)

    def omega(x):
        den = absA @ np.abs(x) + np.abs(b)
        r = np.abs(A @ x - b)
        m = den > 0
        return float(np.max(r[m] / den[m]))
    # NOTE: componentwise omega is degenerate for exactly-decoupled contact rows
    # (one nonzero, b = 0, exact dx = 0): any round-off dx there gives omega = 1.
    res['backward_error_componentwise'] = dict(builtin=omega(dx_builtin), scipy_ext=omega(dx_ext),
                                               callback=omega(x_cb),
                                               note='degenerate on decoupled rows; see normwise_equilibrated')
    # 1-norm condition estimates: raw, and after row/column max-equilibration
    def cond1(M):
        try:
            l = spla.splu(M.tocsc())
            Minv = spla.LinearOperator(M.shape, matvec=l.solve, rmatvec=lambda v: l.solve(v, trans='T'))
            return float(spla.onenormest(M) * spla.onenormest(Minv))
        except Exception as exc:
            return f'failed: {exc!r}'
    Dr = sp.diags(1 / np.asarray(absA.max(axis=1).todense()).ravel())
    As = (Dr @ A).tocsc()
    Dc = sp.diags(1 / np.asarray(abs(As).max(axis=0).todense()).ravel())
    Ae = (As @ Dc).tocsc()
    res['cond1_estimate'] = dict(raw=cond1(A), row_col_equilibrated=cond1(Ae))
    be = Dr @ b
    nAe = float(abs(Ae).sum(axis=1).max())

    def eta(x):           # normwise backward error of the equilibrated system
        xe = x / Dc.diagonal()
        return float(np.max(np.abs(Ae @ xe - be)) / (nAe * np.max(np.abs(xe)) + np.max(np.abs(be))))
    res['backward_error_normwise_equilibrated'] = dict(builtin=eta(dx_builtin), scipy_ext=eta(dx_ext),
                                                       callback=eta(x_cb))
    # expected forward-error scale for two backward-stable solvers: ~ cond * eps
    if isinstance(res['cond1_estimate']['row_col_equilibrated'], float):
        res['forward_error_scale_cond_x_eps'] = res['cond1_estimate']['row_col_equilibrated'] * np.finfo(float).eps

    # ---- contact rows and permutation
    res['permutation'] = dict(n=int(len(prow)), n_moved=int(np.sum(prow != np.arange(len(prow)))),
                              n_dropped=int(np.sum(prow < 0)), n_copy=int(np.sum(pcopy == 1)))
    nodes = D.contact_nodes('drain')
    si_n = eqn[(D.SI, 'Electrons')]; si_p = eqn[(D.SI, 'Potential')]
    nv = np.array(ds.get_node_model_values(device=D.DEV, region=D.SI, name='NodeVolume'))
    cr = []
    Ar = A.tocsr()
    for nd in nodes[:3] + nodes[len(nodes) // 2:len(nodes) // 2 + 1] + nodes[-2:]:
        row = si_n[nd]
        cols = Ar.indices[Ar.indptr[row]:Ar.indptr[row + 1]]; vals = Ar.data[Ar.indptr[row]:Ar.indptr[row + 1]]
        diag = float(vals[cols == row][0]) if np.any(cols == row) else None
        cr.append(dict(node=int(nd), row=int(row), nnz_in_row=int(len(cols)), diag=diag,
                       node_volume=float(nv[nd]), diag_over_node_volume=diag / nv[nd] if diag else None,
                       couples_to_circuit=bool(np.any(cols == eqn[('circuit', 'drain_bias')][0]))
                       if ('circuit', 'drain_bias') in eqn else None,
                       potential_row=int(si_p[nd])))
    res['drain_contact_electron_rows'] = cr
    if ('circuit', 'drain_bias') in eqn:
        crow = int(eqn[('circuit', 'drain_bias')][0])
        res['circuit_rows'] = {k: int(v[0]) for k, v in eqn.items() if k[0] == 'circuit'}
        res['circuit_drain_bias_row_nnz'] = int(Ar.indptr[crow + 1] - Ar.indptr[crow])

    # ---- central-difference Jacobian columns (b = -F)
    restore(X0)
    rng = np.random.default_rng(12345)
    x0, y0 = map(np.array, D.coords())
    samples = []
    special = {'drain_contact': nodes[len(nodes) // 3], 'source_contact': D.contact_nodes('source')[5],
               'interface_si': int(np.where(np.isclose(y0, 0) & np.isclose(x0, 5e-5))[0][0]),
               'body_contact': D.contact_nodes('body')[4]}
    interior = np.where((x0 > 1e-6) & (x0 < 9.9e-5) & (y0 > 1e-7) & (y0 < 1.1e-4))[0]
    for name, nd in list(special.items()) + [(f'interior{i}', int(i)) for i in rng.choice(interior, n_fd_per_block, replace=False)]:
        samples += [(D.SI, 'Potential', nd, name), (D.SI, 'Electrons', nd, name)]
    xo, yo = map(np.array, D.coords(D.OX))
    gate_nd = D.contact_nodes('gate', D.OX)[7]
    ox_int = int(np.where(np.isclose(yo, -1e-6, atol=3e-7) & np.isclose(xo, 5e-5))[0][0])
    samples += [(D.OX, 'Potential', gate_nd, 'gate_contact'), (D.OX, 'Potential', ox_int, 'oxide_interior')]
    if ('circuit', 'drain_bias') in eqn:
        samples += [('circuit', 'drain_bias', 0, 'circuit')]
    fd = []
    for r, var, nd, name in samples:
        col = int(eqn[(r, var)][nd]) if r != 'circuit' else int(eqn[(r, var)][0])
        base = X0[(r, var)][nd]
        h = 1e-6 if var in ('Potential', 'drain_bias') else 1e-6 * abs(base)
        cols = []
        for s in (+1, -1):
            if r == 'circuit':
                ds.set_circuit_node_value(solution='dcop', node=var, value=float(base + s * h))
            else:
                ds.set_node_value(device=D.DEV, region=r, name=var, index=int(nd), value=float(base + s * h))
            cols.append(export()[1])
        if r == 'circuit':
            ds.set_circuit_node_value(solution='dcop', node=var, value=float(base))
        else:
            ds.set_node_value(device=D.DEV, region=r, name=var, index=int(nd), value=float(base))
        fd_col = -(cols[0] - cols[1]) / (2 * h)      # dF/dx_j with F = -b
        a_col = A[:, col].toarray().ravel()
        den = np.max(np.abs(a_col))
        fd.append(dict(unknown=f'{r}:{var}', node=int(nd), where=name, col=col, h=h,
                       max_abs_col=float(den), max_abs_diff=float(np.max(np.abs(fd_col - a_col))),
                       rel=float(np.max(np.abs(fd_col - a_col)) / den) if den else None,
                       rel_if_sign_flipped=float(np.max(np.abs(-fd_col - a_col)) / den) if den else None))
    res['fd_columns'] = fd
    restore(X0)
    np.savez_compressed(os.path.join(outdir, f'{tag}_system.npz'), A_data=A.data, A_indices=A.indices,
                        A_indptr=A.indptr, rhs=b, dx_builtin=dx_builtin, dx_ext=dx_ext)
    # put the circuit back to V1 = 0 and reconverge for the next state
    ds.circuit_alter(name='V1', value=0.0)
    res['reconverge_V1_0'] = D.converge()
    return res


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out', required=True); ap.add_argument('--vg', type=float, default=1.25)
    a = ap.parse_args(); out = os.path.abspath(a.out); os.makedirs(out, exist_ok=False)
    rec = runlog.header(dict(vg=a.vg, test='C-1(b)', pm=D.pm_checklist()))
    t0 = time.time()
    rec['mesh'] = D.build_mesh(); D.build_physics()
    rec['equilibrium'] = D.equilibrium(a.vg)
    D.redeclare_updates('default', 'default')
    rec['undamped_reconverge'] = D.converge()
    rec['devsim_info'] = {k: str(v) for k, v in ds.get_parameter(name='info').items()}
    rec['direct_solver'] = ds.get_parameter(name='direct_solver')
    rec['states'] = {}
    for tag, delta in (('dV1mV', 1e-3), ('dV50mV', 5e-2)):
        rec['states'][tag] = run_state(delta, out, tag)
    rec['wall_s'] = time.time() - t0; rec['end_utc'] = runlog.now()
    runlog.write_json(os.path.join(out, 'c1b.json'), rec)
    print('C1B_DONE', json.dumps({t: dict(A=s['A_export_vs_newton'], rhs=s['rhs_export_vs_newton'],
                                         dx=s['dx_ext_vs_builtin']) for t, s in rec['states'].items()}, default=str))
    return 0


if __name__ == '__main__':
    sys.exit(main())
