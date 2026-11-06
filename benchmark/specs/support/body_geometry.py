"""Geometric factor of the drain-to-body electron path (planning estimate for WP3).
Uniform conductivity 1 in the p region; n+ wells (optionally enlarged by the junction depletion width W_j) are
perfect conductors held at drain = 1, source = 0; Si/SiO2 surface and far edges below the wells Neumann;
bottom y = depth either Dirichlet 0 (ohmic body, B1) or Neumann (blocking body, B0).
Returns per-unit-width currents in units of sigma_bulk (dimensionless geometric factors)."""
import sys, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla

def solve(L, lw, hw, depth, wj, ohmic_body, h=0.005):
    x = np.arange(-lw, L + lw + h / 2, h); y = np.arange(0, depth + h / 2, h)
    nx, ny = len(x), len(y); X, Y = np.meshgrid(x, y, indexing='ij')
    src = (X <= 0 + wj + 1e-12) & (Y <= hw + wj + 1e-12)
    drn = (X >= L - wj - 1e-12) & (Y <= hw + wj + 1e-12)
    body = (np.abs(Y - depth) < 1e-12) & ~src & ~drn if ohmic_body else np.zeros_like(src)
    idx = np.arange(nx * ny).reshape(nx, ny)
    rows, cols = [], []
    rows += list(idx[:-1, :].ravel()); cols += list(idx[1:, :].ravel())
    rows += list(idx[:, :-1].ravel()); cols += list(idx[:, 1:].ravel())
    r, c = np.array(rows), np.array(cols); w = np.ones(len(r))
    W = sp.coo_matrix((np.r_[w, w], (np.r_[r, c], np.r_[c, r])), shape=(nx * ny,) * 2).tocsr()
    Lap = (sp.diags(np.asarray(W.sum(1)).ravel()) - W).tocsr()
    fixed = (src | drn | body).ravel(); phi = np.zeros(nx * ny); phi[drn.ravel()] = 1.0
    fr = np.where(~fixed)[0]; fx = np.where(fixed)[0]
    phi[fr] = spla.spsolve(Lap[fr][:, fr].tocsc(), -Lap[fr][:, fx] @ phi[fx])
    I = Lap @ phi       # unit conductance per edge of a square grid = sigma per unit width
    return dict(I_D=float(I[drn.ravel()].sum()), I_S=float(I[src.ravel()].sum()), I_B=float(I[body.ravel()].sum()),
                energy=float(phi @ (Lap @ phi)))

for name, L, wj in (('dev1 NA1e18 tox2 L1', 1.0, 0.037), ('dev2 NA1e16 tox20 L2', 2.0, 0.351)):
    for lab, w in (('wells only', 0.0), ('wells+W_j', wj)):
        b1 = solve(L, 0.3, 0.5, 1.2, w, True); b0 = solve(L, 0.3, 0.5, 1.2, w, False)
        print(f"{name:22s} {lab:10s} B1: G_DD={b1['I_D']:.4f} G_SD={-b1['I_S']:.4f} G_DB={-b1['I_B']:.4f} (energy {b1['energy']:.4f}) | "
              f"B0: G_DD={b0['I_D']:.4f} G_SD={-b0['I_S']:.4f}   [units of sigma_bulk]")
