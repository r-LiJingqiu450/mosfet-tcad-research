"""Production regression baseline. Metadata, freeze and tests use the same config."""
from pathlib import Path
import json
import numpy as np
from .config import production_metadata
from .phase1 import production_table
from .gca import gca_solve,id_vg
from .postprocess import extract_vt,extract_ss

GOLDEN_PATH=Path(__file__).resolve().parents[1]/'results/golden/stage3.json'


def production_numbers(table=None):
    if table is None:table,_=production_table()
    vg=np.linspace(0,2,81);lin=id_vg(table,vg,.05)
    return dict(vt_cc=extract_vt(vg,lin)[0],ss=extract_ss(vg,lin),
        id_lin_2=gca_solve(2,.05,table)['ID'],id_sat_2=gca_solve(2,1.5,table)['ID'],
        ns_source_2=table.solve_vg(2,0)[1],
        id_subthreshold_0p5=gca_solve(.5,.05,table)['ID'],
        id_transition_1p25=gca_solve(1.25,.2,table)['ID'])


def read_golden(path=GOLDEN_PATH):
    if not Path(path).is_file():
        raise FileNotFoundError('Production S3 golden missing; normal reproduction cannot regenerate it')
    return json.loads(Path(path).read_text(encoding='utf-8'))


def compare_production(current,baseline):
    if baseline['config']!=production_metadata():
        raise AssertionError('Golden/production configuration mismatch')
    ref=baseline['numbers']
    if set(current)!=set(ref):raise AssertionError('Golden key set mismatch')
    for key,value in current.items():
        if not np.isfinite(value) or not np.isfinite(ref[key]):
            raise AssertionError(f'Nonfinite golden metric: {key}')
        if not np.isclose(value,ref[key],rtol=1e-6,atol=0):
            raise AssertionError(f'Golden regression failure: {key}: {value} vs {ref[key]}')
    return True
