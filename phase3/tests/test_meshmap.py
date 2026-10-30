"""Tests for the D-1 mesh-family generator (phase3/meshmap.py)."""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
import pytest

from phase3 import meshmap as MM

DEPTH = MM.f0_depth_map()
OX = MM.oxide_map()
LAT = MM.lateral_map()


@pytest.mark.parametrize('m', [DEPTH, OX, LAT])
@pytest.mark.parametrize('k', range(5))
def test_breakpoints_are_nodes_on_every_level(m, k):
    x = m.nodes(k)
    for bp in m.breakpoints():
        assert np.min(np.abs(x - bp)) == 0.0


@pytest.mark.parametrize('m', [DEPTH, OX, LAT])
@pytest.mark.parametrize('k', range(5))
def test_monotone_and_counts(m, k):
    x = m.nodes(k)
    assert np.all(np.diff(x) > 0)
    assert len(x) - 1 == sum(int(round(s.n0 * math.sqrt(2) ** k)) for s in m.segments)


def test_level0_spec_values():
    q = DEPTH.quality(0)
    segA, segB = q['segments']
    assert abs(segA['h_first'] - 0.25e-7) < 1e-12          # 0.25 nm at the interface
    assert segA['h_last'] <= 20e-7 + 1e-15 and segB['h_last'] <= 20e-7 + 1e-15
    assert OX.quality(0)['n_cells'] >= 5


def test_same_map_on_all_levels():
    """Nodes of a coarse level lie on the fixed map: X(xi) evaluated at their xi."""
    sA = DEPTH.segments[0]
    for k in range(4):
        n = sA.n(k)
        x = sA.nodes(k)
        assert np.allclose(x, sA.X(np.arange(n + 1) / n), rtol=0, atol=1e-18)


def test_effective_refinement_ratio_close_to_sqrt2():
    for m in (DEPTH, OX, LAT):
        for k in range(3):
            for s in m.quality(k)['segments']:
                assert abs(s['refinement_ratio_to_next'] - math.sqrt(2)) < 0.15


def test_adjacent_ratio_changes_with_level_and_tends_to_one():
    """The intra-segment ratio is Q^(1/N_k): level dependent, by construction."""
    sA = DEPTH.segments[0]
    rr = [DEPTH.quality(k)['segments'][0]['max_adjacent_ratio'] for k in range(4)]
    assert all(abs(r - sA.Q ** (1 / sA.n(k))) < 1e-9 for k, r in enumerate(rr))
    assert rr[0] > rr[1] > rr[2] > rr[3] > 1.0


def test_stretch_solver_roundtrip():
    Q = MM.stretch_for(1.0, 0.001, 100)
    s = MM.Segment(0., 1., 100, 'geometric', Q)
    assert abs(np.diff(s.nodes(0))[0] - 0.001) < 1e-12


def test_family_record_is_reproducible(tmp_path):
    a = MM.write_family(tmp_path / 'a.json', [DEPTH, OX])
    b = MM.write_family(tmp_path / 'b.json', [DEPTH, OX])
    assert (tmp_path / 'a.json').read_bytes() == (tmp_path / 'b.json').read_bytes()
    assert a['levels'][3]['si_depth']['quality']['n_nodes'] > a['levels'][0]['si_depth']['quality']['n_nodes']
