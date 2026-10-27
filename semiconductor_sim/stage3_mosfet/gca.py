"""Classical GCA full quasi-Fermi integral (drift and diffusion combined).

The default excludes empirical CLM. At fixed VG, SliceTable's Ns is piecewise
linear in phi_n; integrate that interpolant exactly and invert its cumulative
integral without an intermediate y grid. Both representations are integrated analytically. Exactness for an interpolant does
not imply accuracy for the continuum device. clm=True opts into legacy CLM.
"""
import numpy as np
from common import constants as C


def _areas(knots, values, mode):
    h=np.diff(knots); nl=values[:-1]; nr=values[1:]
    if mode=='linear':
        return .5*(nl+nr)*h
    k=np.log(nr/nl)
    factor=np.ones_like(k)
    np.divide(np.expm1(k),k,out=factor,where=np.abs(k)>1e-12)
    return nl*h*factor


def _evaluate(x, knots, values, mode):
    return (np.exp(np.interp(x,knots,np.log(values)))
            if mode=='log' else np.interp(x,knots,values))


def _integral_between(a, b, knots, values, mode):
    xx=np.r_[a,knots[(knots>a)&(knots<b)],b]
    return float(np.sum(_areas(xx,_evaluate(xx,knots,values,mode),mode)))


def gca_solve(vg, vd, table, y=None, nphi=0, nmin=1e12, clm=False):
    """Solve the classical reduced channel model, not a full 2D DD device.

    id_local is independently finite-differenced and has truncation error.
    id_edges is an integral consistency diagnostic, not independent validation.
    Only nonnegative VD within the table's phi_n domain is supported.
    """
    vd = float(vd)
    if not np.isfinite(vd) or vd < 0 or vd > table.phi_n_grid[-1]:
        raise ValueError('VD outside supported nonnegative table domain')
    y = np.linspace(0., C.L_CH, 161) if y is None else np.asarray(y, float)
    if (y.ndim != 1 or len(y) < 3 or not np.all(np.isfinite(y))
            or np.any(np.diff(y) <= 0) or y[0] < 0 or y[-1] > C.L_CH):
        raise ValueError('y must increase strictly within [0,L]')
    if clm:
        from .legacy_gca import gca_solve as legacy
        r = legacy(vg, vd, table, y=y, nphi=max(int(nphi), 401),
                   nmin=nmin, clm=True)
        r.update(model_mode='legacy_empirical_clm', profile_is_dd_solution=False)
        return r
    knots = np.r_[0., table.phi_n_grid[(table.phi_n_grid > 0)
                                      & (table.phi_n_grid < vd)], vd]
    if int(nphi) > 1:
        knots = np.r_[knots, np.linspace(0., vd, int(nphi))]
    knots = np.unique(knots)
    ns_nodes = table.ns_vs_phi_n(vg, knots)
    if np.any(~np.isfinite(ns_nodes)) or np.any(ns_nodes <= 0):
        raise ValueError('Nonpositive or invalid sheet density')
    mode = getattr(table,"sheet_interpolation","linear")
    if vd == 0:
        phi_n = np.zeros_like(y)
        current = 0.
    else:
        areas = _areas(knots, ns_nodes, mode)
        primitive = np.r_[0., np.cumsum(areas)]
        total = primitive[-1]
        current = C.Q_E*C.W_CH*C.MU_N/C.L_CH*total
        target = total*y/C.L_CH
        idx = np.clip(np.searchsorted(primitive, target, side='right')-1,
                      0, len(knots)-2)
        a = (target-primitive[idx])/np.diff(knots)[idx]
        nl = ns_nodes[idx]
        if mode=='log':
            k=np.log(ns_nodes[idx+1]/nl)
            t=np.zeros_like(a)
            small=np.abs(k)<1e-12
            t[small]=a[small]/nl[small]
            z=k[~small]*a[~small]/nl[~small]
            t[~small]=np.log1p(np.maximum(z,-1+np.finfo(float).eps))/k[~small]
            t=np.clip(t,0.,1.)
        else:
            dn=ns_nodes[idx+1]-nl
            disc=np.maximum(nl*nl+2*dn*a,0.)
            t=np.clip(2*a/(nl+np.sqrt(disc)),0.,1.)
        phi_n = knots[idx]+t*np.diff(knots)[idx]
        phi_n[y == 0.] = 0.
        phi_n[y == C.L_CH] = vd
    ns = _evaluate(phi_n, knots, ns_nodes, mode)
    ps = np.array([table.solve_vg(vg, pn)[0] for pn in phi_n])
    fd = C.Q_E*C.W_CH*C.MU_N*ns*np.gradient(phi_n, y, edge_order=2)
    edges = np.array([C.Q_E*C.W_CH*C.MU_N
                      *_integral_between(a,b,knots,ns_nodes,mode)/dy
                      for a,b,dy in zip(phi_n[:-1],phi_n[1:],np.diff(y))])
    return dict(ID=float(current), vdsat=float('nan'), saturated=False,
                y_p=C.L_CH, y=y, phi_n=phi_n, ns=ns, phi_s=ps,
                id_local=fd, id_edges=edges,
                model_mode='classical_full_integral',
                profile_is_dd_solution=False, integral_knots=len(knots))


def id_vg(table, vg_list, vd, **kw):
    return np.array([gca_solve(vg,vd,table,**kw)['ID'] for vg in vg_list])


def id_vd(table, vg, vd_list, **kw):
    rows = [gca_solve(vg,vd,table,**kw) for vd in vd_list]
    return tuple(np.array([r[k] for r in rows]) for k in ['ID','vdsat','y_p'])


def transfer_characteristics(table, vg_list, vd_lin=.05, vd_sat=1.5, **kw):
    return dict(vg=np.asarray(vg_list,float),
                id_lin=id_vg(table,vg_list,vd_lin,**kw),
                id_sat=id_vg(table,vg_list,vd_sat,**kw))


def output_characteristics(table, vg_list, vd_list, **kw):
    rows = [id_vd(table,vg,vd_list,**kw) for vg in vg_list]
    return dict(vg=np.asarray(vg_list),vd=np.asarray(vd_list),
                **{k:np.array([r[i] for r in rows])
                   for i,k in enumerate(['ids','vdsats','yps'])})
