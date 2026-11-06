"""Phase 3 research modules. Import production code read-only; never modify it."""
import os as _os
import sys as _sys

_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_SIM = _os.path.join(_ROOT, "semiconductor_sim")
for _p in (_ROOT, _SIM):
    if _p not in _sys.path:
        _sys.path.insert(0, _p)
