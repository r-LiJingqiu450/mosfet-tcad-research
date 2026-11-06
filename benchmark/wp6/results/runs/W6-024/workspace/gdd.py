"""Zero-bias drain conductance, three ways:
 (1) linear-response resistor network from the converged equilibrium (exact linearization of the SG
     discretization at phi_n = const: dJ_e = q mu couple/h * n0 B(v) * (dphi0 - dphi1), decoupled from Poisson),
     G evaluated by contact fluxes and by the (variational) dissipation sum;
 (2) DEVSIM small-signal AC at low frequency (real part of V1.I);
 (3) DEVSIM DC central finite differences.
"""
import math, sys, json
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import mosfet_deck as M
import devsim as ds

info = M.build()
print(info, flush=True)
D, R = M.DEV, M.SI

x, y = M.coords(R)
n = np.array(ds.get_node_model_values(device=D, region=R, name='Electrons'))
psi = np.array(ds.get_node_model_values(device=D, region=R, name='Potential'))
phin = psi - M.VT * np.log(n / M.NI)
print('equilibrium phi_n range', phin.min(), phin.max())

ds.edge_from_node_model(device=D, region=R, node_model='node_index')
i0 = np.array(ds.get_edge_model_values(device=D, region=R, name='node_index@n0')).astype(int)
i1 = np.array(ds.get_edge_model_values(device=D, region=R, name='node_index@n1')).astype(int)
cpl = np.array(ds.get_edge_model_values(device=D, region=R, name='EdgeCouple'))
hinv = np.array(ds.get_edge_model_values(device=D, region=R, name='EdgeInverseLength'))
B = np.array(ds.get_edge_model_values(device=D, region=R, name='Bern01'))
vd = np.array(ds.get_edge_model_values(device=D, region=R, name='vdiff'))
Bm = B + vd  # B(-v)
neff = 0.5 * (n[i0] * B + n[i1] * Bm)
print('max rel asym of n0B(v) vs n1B(-v):', np.max(np.abs(n[i0] * B - n[i1] * Bm) / neff))
g = M.Q * M.MU_N * hinv * cpl * neff
N = len(x)

ext = M.MESH['contact_extent_cm'] * (1 + 1e-12)
src = (np.abs(x - x.min()) < 1e-15) & (y <= ext)
drn = (np.abs(x - x.max()) < 1e-15) & (y <= ext)
bod = np.abs(y - y.max()) < 1e-15
print('contact node counts', src.sum(), drn.sum(), bod.sum())
fixed = src | drn | bod
val = np.where(drn, 1.0, 0.0)

L = sp.coo_matrix((np.concatenate([g, g, -g, -g]),
                   (np.concatenate([i0, i1, i0, i1]), np.concatenate([i0, i1, i1, i0]))), shape=(N, N)).tocsr()
free = np.where(~fixed)[0]; fx = np.where(fixed)[0]
A = L[free][:, free].tocsc()
rhs = -L[free][:, fx] @ val[fx]
# diagonal scaling for conditioning
d = np.sqrt(A.diagonal())
Dm = sp.diags(1 / d)
u = spla.spsolve((Dm @ A @ Dm).tocsc(), rhs / d) / d
phi = val.copy(); phi[free] = u
# iterative refinement
for _ in range(3):
    r = rhs - A @ phi[free]
    du = spla.spsolve((Dm @ A @ Dm).tocsc(), r / d) / d
    phi[free] += du
print('phi range', phi.min(), phi.max())
Je = g * (phi[i0] - phi[i1])  # current n0 -> n1 (phi difference driving)
inj = np.zeros(N); np.add.at(inj, i0, Je); np.add.at(inj, i1, -Je)
G_drain = inj[drn].sum(); G_src = -inj[src].sum(); G_body = -inj[bod].sum()
G_energy = np.sum(g * (phi[i0] - phi[i1]) ** 2)
print(f'network: G_drain={G_drain:.10e} G_src+body={G_src + G_body:.10e} (src {G_src:.6e}, body {G_body:.6e})'
      f' G_energy={G_energy:.10e}', flush=True)
print('max free-node residual / max |inj|', np.max(np.abs(inj[free])), flush=True)

# (2) AC
for f in (1e-6, 1.0):
    ds.solve(type='ac', frequency=f)
    print('AC f=', f, 'V1.I re', ds.get_circuit_node_value(node='V1.I', solution='ssac_real'),
          'im', ds.get_circuit_node_value(node='V1.I', solution='ssac_imag'), flush=True)


def currents():
    out = {c: ds.get_contact_current(device=D, contact=c, equation='ElectronContinuityEquation')
           for c in ('source', 'drain', 'body')}
    out['V1.I'] = ds.get_circuit_node_value(node='V1.I', solution='dcop')
    return out

print('V=0', currents(), flush=True)
res = {}
for dv in (1e-3, -1e-3, 5e-3, -5e-3):
    ds.circuit_alter(name='V1', value=dv)
    M.solve_dc()
    res[dv] = currents()
    print('V=', dv, res[dv], flush=True)
for a in (1e-3, 5e-3):
    print(a, {k: (res[a][k] - res[-a][k]) / (2 * a) for k in res[a]}, flush=True)
