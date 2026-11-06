"""Geometric factor of the hole path through one n+ well (planning estimate for WP3).
Hole conductivity in the well sigma_p = q mu_p n_i^2/N_D is tiny; the p bulk is ~1e18x more conductive and is
treated as a perfect conductor at 0 (connected to the body or to the other well). Well: 0 <= x' <= l_w (x' measured
from the far-edge contact), 0 <= y <= h_w; contact x' = 0 at 1; junction face x' = l_w and well bottom y = h_w at 0;
oxide top y = 0 Neumann. Returns G/sigma_p (per unit width)."""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla
def g(lw, hw, h=0.0025):
    x = np.arange(0, lw + h / 2, h); y = np.arange(0, hw + h / 2, h); nx, ny = len(x), len(y)
    X, Y = np.meshgrid(x, y, indexing='ij'); idx = np.arange(nx * ny).reshape(nx, ny)
    r = np.r_[idx[:-1, :].ravel(), idx[:, :-1].ravel()]; c = np.r_[idx[1:, :].ravel(), idx[:, 1:].ravel()]
    W = sp.coo_matrix((np.ones(2 * len(r)), (np.r_[r, c], np.r_[c, r])), shape=(nx * ny,) * 2).tocsr()
    Lap = (sp.diags(np.asarray(W.sum(1)).ravel()) - W).tocsr()
    one = (X == 0).ravel(); zero = ((np.abs(X - lw) < 1e-12) | (np.abs(Y - hw) < 1e-12)).ravel() & ~one
    phi = np.zeros(nx * ny); phi[one] = 1; fx = one | zero; fr = np.where(~fx)[0]
    phi[fr] = spla.spsolve(Lap[fr][:, fr].tocsc(), -Lap[fr][:, np.where(fx)[0]] @ phi[fx])
    return float((Lap @ phi)[one].sum())
q, mu_p, ni, ND = 1.602176634e-19, 450.0, 1e10, 1e20
gf = g(0.3, 0.5); sp_ = q * mu_p * ni * ni / ND
print(f"geometric factor {gf:.3f}; sigma_p,well {sp_:.3e} S/cm; hole conductance per well {gf*sp_:.3e} S/cm")
for lab, G0 in (('dev1 psi12', 3.577e-12), ('dev1 psi16', 1.649e-10), ('dev2 psi12', 2.635e-9)):
    print(f"  relative to channel G0 ({lab}, deck-W value {G0:.3e}): {gf*sp_/G0:.2e}")
