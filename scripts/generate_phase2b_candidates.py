#!/usr/bin/env python3
"""Generate versioned S1/S2 candidate baselines.

These files are explicitly separate from the official production golden
directory.  The script never writes semiconductor_sim/results/golden.
"""

from __future__ import annotations

import hashlib
import json
import platform
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from semiconductor_sim.tests.golden_helpers import stage1_numbers, stage2_numbers  # noqa: E402


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    out = ROOT / "verification/candidates/phase2b_v1_0"
    out.mkdir(parents=True, exist_ok=True)
    official = ROOT / "semiconductor_sim/results/golden/stage3.json"
    with warnings.catch_warnings(record=True) as s1_caught:
        warnings.simplefilter("always")
        s1 = stage1_numbers()
    with warnings.catch_warnings(record=True) as s2_caught:
        warnings.simplefilter("always")
        s2 = stage2_numbers()
    common = {
        "candidate_status": "candidate-not-official-golden",
        "candidate_version": "phase2b_v1.0",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "repository_baseline": "fb7a311303f3b366b556d2a5385def0b9c256399",
        "official_golden_directory": "semiconductor_sim/results/golden",
        "official_s1_s2_status": "missing; not reconstructed",
        "official_stage3_sha256": _sha(official),
        "source_configuration": "semiconductor_sim/tests/golden_helpers.py",
        "candidate_tolerance": {"relative": 1e-6, "absolute": 1e-12},
        "independent_evidence": {
            "s1": "verification/phase2b_s1_verification.json",
            "s2": "verification/phase2b_s2_verification.json",
        },
    }
    s1_doc = {**common, "stage": "stage1", "configuration": {
        "NGRID": 201, "verbose": False, "source_function": "stage1_numbers"},
              "numbers": s1, "warning_summary": sorted({str(item.message) for item in s1_caught})}
    s2_doc = {**common, "stage": "stage2", "configuration": {
        "L_m": 2.0e-6, "v_forward": [0.5], "v_reverse": [-1.0],
        "verbose": False, "source_function": "stage2_numbers"},
              "numbers": s2, "warning_summary": sorted({str(item.message) for item in s2_caught})}
    (out / "stage1.json").write_text(json.dumps(s1_doc, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (out / "stage2.json").write_text(json.dumps(s2_doc, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    manifest = {
        "candidate_status": "candidate-not-official-golden",
        "candidate_version": "phase2b.v1.0",
        "official_golden_untouched": True,
        "official_stage3_sha256": _sha(official),
        "missing_official_files": ["stage1.json", "stage2.json"],
        "files": {
            "stage1.json": _sha(out / "stage1.json"),
            "stage2.json": _sha(out / "stage2.json"),
        },
        "tolerance": common["candidate_tolerance"],
        "provenance": common,
    }
    (out / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(out), "stage1_sha256": manifest["files"]["stage1.json"],
                      "stage2_sha256": manifest["files"]["stage2.json"],
                      "warnings": {"stage1": len(s1_caught), "stage2": len(s2_caught)}}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
