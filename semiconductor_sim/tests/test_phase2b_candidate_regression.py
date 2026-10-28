"""Regression checks for versioned Phase 2B candidate baselines.

These tests intentionally target verification/candidates/ and never the
official semiconductor_sim/results/golden directory.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import warnings

import pytest

from semiconductor_sim.tests.golden_helpers import stage1_numbers, stage2_numbers


ROOT = Path(__file__).resolve().parents[2]
CANDIDATE = ROOT / "verification/candidates/phase2b_v1_0"


def _load(name):
    return json.loads((CANDIDATE / name).read_text(encoding="utf-8"))


def _assert_numbers(actual, expected):
    assert set(actual) == set(expected)
    for key, value in expected.items():
        if isinstance(value, int):
            assert actual[key] == value, key
        else:
            assert actual[key] == pytest.approx(value, rel=1e-6, abs=1e-12), key


def test_candidate_manifest_is_versioned_and_official_golden_is_separate():
    manifest = _load("manifest.json")
    assert manifest["candidate_status"] == "candidate-not-official-golden"
    assert manifest["official_golden_untouched"] is True
    assert manifest["missing_official_files"] == ["stage1.json", "stage2.json"]
    official = ROOT / "semiconductor_sim/results/golden/stage3.json"
    assert hashlib.sha256(official.read_bytes()).hexdigest() == manifest["official_stage3_sha256"]
    for filename, digest in manifest["files"].items():
        assert hashlib.sha256((CANDIDATE / filename).read_bytes()).hexdigest() == digest


def test_stage1_candidate_regression():
    expected = _load("stage1.json")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        actual = stage1_numbers()
    assert not caught
    _assert_numbers(actual, expected["numbers"])


def test_stage2_candidate_regression_and_gummel_warning_is_visible():
    expected = _load("stage2.json")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        actual = stage2_numbers()
    messages = [str(item.message) for item in caught]
    assert any("gummel_solve" in message for message in messages)
    assert any("solve_poisson_pn" in message for message in messages)
    _assert_numbers(actual, expected["numbers"])
