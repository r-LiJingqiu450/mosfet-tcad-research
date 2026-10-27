"""Release-blocking evidence regressions, independent from the new S3 golden."""
import json
from pathlib import Path
import numpy as np
import pytest
from semiconductor_sim.common.jsonio import dumps,write_json


def test_numpy_json_file_and_console_use_same_encoder(tmp_path):
    value=dict(flag=np.bool_(True),integer=np.int64(7),array=np.array([.1,.2]))
    path=tmp_path/'result.json';write_json(path,value)
    assert json.loads(path.read_text())==json.loads(dumps(value))


@pytest.mark.parametrize('bad',[float('nan'),float('inf'),np.float64('-inf')])
def test_nonfinite_json_is_rejected_without_partial_output(tmp_path,bad):
    path=tmp_path/'bad.json'
    with pytest.raises(ValueError):write_json(path,dict(bad=bad))
    assert not path.exists()


def test_background_decomposition_does_not_dominate_reported_subthreshold_points():
    from semiconductor_sim.stage3_mosfet.reference import reference_current
    a=reference_current(.5,.05);b=reference_current(.75,.05)
    for r in [a,b]:assert np.isclose(r['ID'],r['ID_background']+r['ID_excess'],rtol=1e-14)
    assert a['background_fraction']<.001
    assert b['background_fraction']<1e-6
    # A real limitation remains very close to flatband.
    assert reference_current(.05,.05)['background_fraction']>.9


def test_background_is_explicitly_depth_dependent_but_excess_is_not():
    from semiconductor_sim.stage3_mosfet.reference import reference_current
    a=reference_current(.5,.05,x_max=1.2e-6);b=reference_current(.5,.05,x_max=2.4e-6)
    assert np.isclose(b['ID_background'],2*a['ID_background'],rtol=1e-12)
    assert np.isclose(b['ID_excess'],a['ID_excess'],rtol=1e-12)


def test_clipping_sensitivity_is_zero_inside_physical_working_range():
    from semiconductor_sim.common.grid import make_surface_grid
    from semiconductor_sim.stage3_mosfet.solver import solve_poisson_slice
    x=make_surface_grid(1.2e-6,d0=.125e-9,h_max=4e-9)
    a=solve_poisson_slice(.85,0,x);b=solve_poisson_slice(.85,0,x,density_clip_limit=3.)
    assert a['converged'] and b['converged']
    assert b['clip_evaluations']==0
    assert np.array_equal(a['phi'],b['phi'])


def test_density_clipping_is_disabled_and_overflow_is_explicit():
    from semiconductor_sim.stage3_mosfet.solver import _carriers
    p,n=_carriers(np.array([3.1]),np.array([0.]),1e22,1e10)
    pc,nc=_carriers(np.array([3.1]),np.array([0.]),1e22,1e10,3.)
    assert n[0]>nc[0] and p[0]<pc[0]
    with pytest.raises(ValueError,match='floating-point'):_carriers(np.array([100.]),0,1e22,1e10)


def test_production_table_and_golden_metadata_share_one_configuration():
    from semiconductor_sim.stage3_mosfet.config import PRODUCTION,production_metadata
    from semiconductor_sim.stage3_mosfet.golden import read_golden
    from semiconductor_sim.stage3_mosfet.phase1 import production_table
    t,_=production_table()
    assert (t.nps,t.npn,t.sheet_interpolation)==(PRODUCTION.nps,PRODUCTION.npn,PRODUCTION.interpolation)
    assert read_golden()['config']==production_metadata()
    assert t.build_stats['potential_min']>-3 and t.build_stats['potential_max']<3


def test_charge_comparison_is_labelled_algebra_not_physics_validation():
    from semiconductor_sim.stage3_mosfet.reference import charge_budget,reference_slice
    r=reference_slice(2.,0.)
    budget=charge_budget(2.,r['phi_s'],r['ns'],r['qs'])
    assert budget['evidence_type']=='algebraic_reference_comparison'
