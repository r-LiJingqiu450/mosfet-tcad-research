"""Production S3 baseline; unavailable original S1/S2 baselines explicitly skipped."""
import pytest
from semiconductor_sim.stage3_mosfet.golden import (production_numbers,read_golden,compare_production)

@pytest.mark.parametrize('stage',['stage1','stage2','stage3'])
def test_golden_regression(stage):
    if stage!='stage3':pytest.skip('Original '+stage+' golden absent in uploaded source; not reconstructed')
    assert compare_production(production_numbers(),read_golden())

@pytest.mark.parametrize('stage',['stage1','stage2','stage3'])
def test_golden_keys_stable(stage):
    if stage!='stage3':pytest.skip('Original '+stage+' golden absent in uploaded source; not reconstructed')
    assert set(production_numbers())==set(read_golden()['numbers'])


def test_golden_checks_configuration_and_values():
    import copy
    baseline=read_golden();current=dict(baseline['numbers'])
    wrong=copy.deepcopy(baseline);wrong['config']['table']['nps']=240
    with pytest.raises(AssertionError,match='configuration'):compare_production(current,wrong)
    current['id_lin_2']*=1.01
    with pytest.raises(AssertionError,match='regression'):compare_production(current,baseline)
