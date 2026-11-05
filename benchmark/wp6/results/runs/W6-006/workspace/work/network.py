# Exact linear response at equilibrium: SG edge current J = g_e (exp(-phi1/Vt)-exp(-phi0/Vt)) with phi=0 at eq.
# => linear resistor network for delta phi_n, conductance per edge G_e = q mu n_i B(v) exp(psi0/Vt) * EdgeCouple*EdgeInverseLength
import sys, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
sys.path.insert(0, '.')
import mosfet_deck as M
import devsim as ds
M.build()
D, R = M.DEV, M.SI
ds.edge_from_node_model(device=D, region=R, node_model='node_index')
g = lambda n, r=R: np.array(ds.get_edge_model_values(device=D, region=r, name=n))
i0 = g('node_index@n0').astype(int); i1 = g('node_index@n1').astype(int)
psi = np.array(ds.get_node_model_values(device=D, region=R, name='Potential'))
nel = np.array(ds.get_node_model_values(device=D, region=R, name='Electrons'))
VT, NI, Q, MU = M.VT, M.NI, M.Q, M.MU_N
couple = g('EdgeCouple'); il = g('EdgeInverseLength')
v = (psi[i0] - psi[i1]) / VT
from scipy.special import exprel
B = 1.0 / exprel(v)       # v/(e^v-1)
# use the solved electron density (eq): B(v)*n0 ; symmetric form check
Ge = Q * MU * VT * il * couple * B * nel[i0] / VT
Ge_alt = Q * MU * VT * il * couple * B * np.exp(v) * nel[i1] / VT
print('symmetry check max rel diff', np.max(np.abs(Ge - Ge_alt) / np.maximum(Ge, 1e-300)))
Ge_th = Q * MU * il * couple * B * NI * np.exp(psi[i0] / VT)
print('n vs n_i exp(psi/VT) max rel', np.max(np.abs(Ge - Ge_th) / np.maximum(Ge, 1e-300)))
N = len(psi)
Lp = sp.coo_matrix((np.r_[Ge, Ge, -Ge, -Ge], (np.r_[i0, i1, i0, i1], np.r_[i0, i1, i1, i0])), shape=(N, N)).tocsr()
def cnodes(c):
    return np.array(sorted(set(int(k) for e in ds.get_element_node_list(device=D, region=R, contact=c) for k in e)))
dn = cnodes('drain'); sn = cnodes('source'); bn = cnodes('body')
print('contact nodes', len(dn), len(sn), len(bn))
fixed = np.zeros(N, bool); fixed[dn] = fixed[sn] = fixed[bn] = True
phi = np.zeros(N); phi[dn] = 1.0
free = ~fixed
A = Lp[free][:, free]; b = -Lp[free][:, fixed] @ phi[fixed]
phi[free] = spl.spsolve(A.tocsc(), b)
I = Lp @ phi
print('residual free', np.max(np.abs(I[free])) / np.max(np.abs(I)))
print('G_DD (drain flux)  =', I[dn].sum())
print('source flux', I[sn].sum(), 'body flux', I[bn].sum(), 'sum', I[dn].sum()+I[sn].sum()+I[bn].sum())
print('energy phi^T L phi =', phi @ (Lp @ phi))
np.savez('work/network.npz', phi=phi, psi=psi)
print('Ge range', Ge.min(), Ge.max(), 'neg', (Ge<0).sum(), 'couple min', couple.min())
r = np.abs(I); r[fixed]=0; k=np.argmax(r); xs,ys=M.coords(R)
print('worst node', k, xs[k], ys[k], I[k], phi[k], 'diag', Lp[k,k])
# symmetric diagonal scaling + refinement
d = np.sqrt(A.diagonal()); Dm = sp.diags(1/d)
As = (Dm @ A @ Dm).tocsc(); bs = b / d
lu = spl.splu(As)
y = lu.solve(bs)
for it in range(5):
    res = bs - As @ y; y += lu.solve(res); print('refine', it, np.max(np.abs(res)))
phi2 = np.zeros(N); phi2[dn]=1; phi2[free] = y / d
I2 = Lp @ phi2
print('scaled: residual', np.max(np.abs(I2[free])), 'drain', I2[dn].sum(), 'source', I2[sn].sum(), 'body', I2[bn].sum(), 'energy', phi2 @ (Lp@phi2))
np.savez('work/network.npz', phi=phi2, psi=psi)
# complementary problem: drain at 0, source and body at 1, so the drain flux is computed where phi is near 0
phc = np.zeros(N); phc[sn] = 1; phc[bn] = 1
bc = -(Lp[free][:, fixed] @ phc[fixed]) / d
y = lu.solve(bc)
for it in range(3): y += lu.solve(bc - As @ y)
phc[free] = y / d
Ic = Lp @ phc
print('complementary: drain flux', Ic[dn].sum())
print('G_DD from problem1 = -(I_S+I_B) =', -(I2[sn].sum() + I2[bn].sum()))
