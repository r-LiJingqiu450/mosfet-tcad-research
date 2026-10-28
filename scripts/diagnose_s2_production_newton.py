#!/usr/bin/env python3
"""Record the production S2 full-Newton equilibrium trace."""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "semiconductor_sim"))

from common.constants import NA_CM3, ND_CM3  # noqa: E402
from common.grid import make_junction_grid  # noqa: E402
from stage2_pn.solver import newton_solve  # noqa: E402


def main() -> int:
    x = make_junction_grid(2.0e-6)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = newton_solve(x, NA_CM3 * 1e6, ND_CM3 * 1e6, 0.0,
                              tol=1e-9, verbose=True)
    print("RESULT", {"iters": result["iters"],
                     "converged": result["converged"]})
    print("WARNINGS", [str(item.message) for item in caught])
    return 0 if result["converged"] and not caught else 1


if __name__ == "__main__":
    raise SystemExit(main())
