"""D-1 mesh-family generator (independent replacement, 2026-10-30).

The 2026-10-29 `phase3/meshmap.py` was not recoverable (filesystem, all Git
history and Project documents were searched). This module is a new
implementation with its own tests (`phase3/tests/test_meshmap.py`).

Construction (R04 v2.1 section 7.3)
-----------------------------------
A coordinate direction is a sequence of SEGMENTS between FIXED breakpoints.
Each segment s has a fixed monotone map  x = X_s(xi), xi in [0, 1]:

  * 'uniform'   : X(xi) = a + (b - a) xi
  * 'geometric' : X(xi) = a + (b - a) (Q^xi - 1)/(Q - 1)   (Q > 0, Q != 1)
                  Q is the total stretch (last/first local spacing ratio of the
                  continuous map, Q = X'(1)/X'(0)); Q < 1 grades the other way.

Level k samples every segment at N_{s,k} + 1 uniform xi points with
    N_{s,k} = round(N_{s,0} * r^k),  r = sqrt(2).
Breakpoints are therefore nodes on every level, and the map itself (not the
cell ratios) is shared by all levels. Consequences, stated explicitly:
  * adjacent-cell ratio inside a geometric segment = Q^(1/N_{s,k}); it CHANGES
    with k (it tends to 1). It is NOT required to be level-independent;
  * effective refinement ratio of a segment is N_{s,k+1}/N_{s,k}, close to r but
    not exactly r because of rounding; it is recorded per segment and level;
  * no convergence order is assumed anywhere in this module.

Quality metrics recorded per level: cell count, min/max spacing, max adjacent
spacing ratio (inside segments and across breakpoints), representative
spacing h = (length / N), and the effective refinement ratio to the next level.

The F0 depth map follows R05 C-2: breakpoints 0, 3 W_dep, x_max; level-0 first
spacing at the interface 0.25 nm; body spacing <= 20 nm; oxide >= 5 cells.
"""
from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, field
from typing import List

import numpy as np

R_REFINE = math.sqrt(2.0)


@dataclass(frozen=True)
class Segment:
    a: float
    b: float
    n0: int
    kind: str = 'uniform'       # 'uniform' | 'geometric'
    Q: float = 1.0              # total stretch for 'geometric'

    def X(self, xi):
        xi = np.asarray(xi, float)
        if self.kind == 'uniform' or abs(self.Q - 1.0) < 1e-15:
            return self.a + (self.b - self.a) * xi
        if self.kind == 'geometric':
            return self.a + (self.b - self.a) * np.expm1(xi * math.log(self.Q)) / (self.Q - 1.0)
        raise ValueError(self.kind)

    def n(self, k, r=R_REFINE):
        return int(round(self.n0 * r ** k))

    def nodes(self, k, r=R_REFINE):
        n = self.n(k, r)
        x = self.X(np.linspace(0.0, 1.0, n + 1))
        x[0], x[-1] = self.a, self.b         # breakpoints exact on every level
        return x


@dataclass(frozen=True)
class AxisMap:
    name: str
    segments: tuple
    units: str = 'cm'

    def nodes(self, k, r=R_REFINE):
        out = [self.segments[0].nodes(k, r)]
        for s in self.segments[1:]:
            out.append(s.nodes(k, r)[1:])
        x = np.concatenate(out)
        return x

    def breakpoints(self):
        return [self.segments[0].a] + [s.b for s in self.segments]

    def quality(self, k, r=R_REFINE):
        x = self.nodes(k, r)
        h = np.diff(x)
        ratios = h[1:] / h[:-1]
        seg = []
        for s in self.segments:
            hs = np.diff(s.nodes(k, r))
            rs = hs[1:] / hs[:-1] if len(hs) > 1 else np.array([1.0])
            seg.append(dict(a=s.a, b=s.b, kind=s.kind, Q=s.Q, n=s.n(k, r),
                            h_first=float(hs[0]), h_last=float(hs[-1]),
                            max_adjacent_ratio=float(np.max(np.maximum(rs, 1 / rs))),
                            h_rep=float((s.b - s.a) / s.n(k, r)),
                            refinement_ratio_to_next=s.n(k + 1, r) / s.n(k, r)))
        joints = []
        idx = 0
        for s in self.segments[:-1]:
            idx += s.n(k, r)
            joints.append(dict(at=float(x[idx]), ratio_right_over_left=float(h[idx] / h[idx - 1])))
        return dict(axis=self.name, level=k, n_nodes=int(len(x)), n_cells=int(len(h)),
                    h_min=float(h.min()), h_max=float(h.max()),
                    max_adjacent_ratio=float(np.max(np.maximum(ratios, 1 / ratios))),
                    monotone=bool(np.all(h > 0)), segments=seg, breakpoint_joints=joints)

    def to_dict(self):
        return dict(name=self.name, units=self.units, segments=[asdict(s) for s in self.segments])


def stretch_for(length, h_first, n):
    """Q such that a geometric segment of `length` with n cells has first cell h_first
    (on the level that uses n cells). Solved by bisection on log Q."""
    f = lambda lq: length * math.expm1(lq / n) / math.expm1(lq) - h_first if abs(lq) > 1e-12 \
        else length / n - h_first
    lo, hi = -50.0, 50.0
    if f(lo) * f(hi) > 0:
        raise ValueError('no stretch reproduces the requested first spacing')
    for _ in range(300):
        mid = .5 * (lo + hi)
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return math.exp(.5 * (lo + hi))


# ----------------------------------------------------------------- F0 family
def wdep_max_cm(na_cm3=1e16, ni_cm3=1e10, eps_r=11.7, vt=0.025851999786435316):
    """Maximum depletion width sqrt(2 eps 2psi_F / (q N_A)) in cm (planning value)."""
    eps = eps_r * 8.8541878128e-14
    psif = vt * math.log(na_cm3 / ni_cm3)
    return math.sqrt(2 * eps * 2 * psif / (1.602176634e-19 * na_cm3))


def f0_depth_map(h0_nm=0.25, hmax_nm=20.0, xmax_um=1.2):
    """Si depth map (x = 0 at the Si/SiO2 interface, increasing into the substrate).
    Segment A [0, 3 W_dep]: geometric, first spacing h0 at level 0 and last spacing
    close to hmax; segment B [3 W_dep, x_max]: uniform with spacing <= hmax."""
    xm = xmax_um * 1e-4
    xb = 3 * wdep_max_cm()
    h0, hm = h0_nm * 1e-7, hmax_nm * 1e-7
    # choose n_A so that the level-0 last spacing of segment A is closest to (<=) hmax
    best = None
    for n in range(20, 2000):
        Q = stretch_for(xb, h0, n)
        last = xb * (Q ** 1.0 - Q ** ((n - 1) / n)) / (Q - 1)
        if last <= hm:
            best = (n, Q); break
    nA, QA = best
    nB = math.ceil((xm - xb) / hm)
    return AxisMap('si_depth', (Segment(0.0, xb, nA, 'geometric', QA), Segment(xb, xm, nB, 'uniform')))


def oxide_map(tox_nm=20.0, n0=5):
    """Oxide from the gate (-t_ox) to the interface (0), uniform (>= 5 cells)."""
    return AxisMap('oxide', (Segment(-tox_nm * 1e-7, 0.0, n0, 'uniform'),))


def lateral_map(L_um=1.0, n0=50):
    """Channel direction, uniform (R05 C-3: y-map uniform, N_y,0 = 50)."""
    return AxisMap('lateral', (Segment(0.0, L_um * 1e-4, n0, 'uniform'),))


def family_record(maps, levels=range(4), r=R_REFINE):
    rec = dict(r=r, maps=[m.to_dict() for m in maps], levels={})
    for k in levels:
        rec['levels'][k] = {m.name: dict(quality=m.quality(k, r), nodes=m.nodes(k, r).tolist())
                            for m in maps}
    return rec


def write_family(path, maps, levels=range(4), r=R_REFINE):
    rec = family_record(maps, levels, r)
    with open(path, 'w') as fh:
        json.dump(rec, fh, indent=1)
    return rec
