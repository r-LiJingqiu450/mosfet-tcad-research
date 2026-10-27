"""Fail closed on missing outputs, incomplete pytest, invalid manifests or failed study."""
import argparse,csv,hashlib,json,sys,xml.etree.ElementTree as ET,zipfile
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from semiconductor_sim.common.jsonio import write_json
OUT=Path('semiconductor_sim/results/phase1')
FIGURES=['fig1_continuity','fig2_independent_reference','fig3_clm_ablation','fig4_transfer',
         'fig5_spatial_refinement','fig6_background','fig7_error_attribution']
REQUIRED=[OUT/n for n in ['phase1_results.json','baseline_recomputed.json','independent_reference.csv',
    'continuity.csv','clm_ablation.csv','spatial_refinement.csv','transfer_curves.npz',
    'error_attribution.json','error_attribution.csv','neutral_background.csv']]
REQUIRED += [OUT/(n+'.'+ext) for n in FIGURES for ext in ['png','svg']]
REQUIRED += [Path(p) for p in ['semiconductor_sim/results/golden/stage3.json','tests_phase1.txt','study_log.txt',
    'verification/pytest.xml','verification/pytest_status.json','reports/Phase1_Report.md',
    'reports/S3_Replacement_Sections.md','README.md','requirements.txt','environment_verified.txt','changes.patch']]


def parse_pytest(root):
    root=Path(root);xml=root/'verification/pytest.xml'
    tree=ET.parse(xml);cases=tree.findall('.//testcase')
    if not cases:raise AssertionError('No completed test cases in JUnit')
    failures=sum(c.find('failure') is not None for c in cases)
    errors=sum(c.find('error') is not None for c in cases)
    skipped=[c for c in cases if c.find('skipped') is not None]
    if failures or errors:raise AssertionError('JUnit contains failures/errors')
    expected_skip={('test_golden_regression[stage1]'),('test_golden_regression[stage2]'),
                   ('test_golden_keys_stable[stage1]'),('test_golden_keys_stable[stage2]')}
    if {c.attrib['name'] for c in skipped}!=expected_skip:
        raise AssertionError('Unexpected skipped tests')
    passed=len(cases)-len(skipped)
    log=(root/'tests_phase1.txt').read_text(encoding='utf-8')
    if f'{passed} passed' not in log or f'{len(skipped)} skipped' not in log:
        raise AssertionError('Test log lacks a matching completed summary')
    return dict(collected=len(cases),passed=passed,skipped=len(skipped),failures=failures,errors=errors,
                skipped_names=sorted(c.attrib['name'] for c in skipped))


def validate(root=ROOT,allow_running=False):
    root=Path(root)
    for rel in REQUIRED:
        if not (root/rel).is_file() or (root/rel).stat().st_size==0:
            raise AssertionError(f'Missing/empty deliverable: {rel}')
    results=json.loads((root/OUT/'phase1_results.json').read_text())
    if not all(results['acceptance'].values()):raise AssertionError('Numerical acceptance failed')
    if (root/'semiconductor_sim/results/golden_phase1.json').exists():
        raise AssertionError('Obsolete disconnected golden still present')
    parse_pytest(root)
    state=json.loads((root/'verification/run_status.json').read_text())
    if not allow_running and (state['state']!='completed' or state['exit_code']!=0):
        raise AssertionError('Reproduction did not finish successfully')
    if state.get('golden_sha256_before')!=state.get('golden_sha256_after'):
        raise AssertionError('Golden modified during reproduction')
    rows=list(csv.DictReader((root/OUT/'independent_reference.csv').open(encoding='utf-8-sig')))
    if len(rows)!=42:raise AssertionError('Expected 42 independent-reference rows')
    with np.load(root/OUT/'transfer_curves.npz',allow_pickle=False) as data:
        if set(data.files)!={'vg','id_lin','id_sat'}:raise AssertionError('Invalid transfer keys')
        if any(data[k].shape!=(81,) or not np.all(np.isfinite(data[k])) for k in data.files):
            raise AssertionError('Invalid transfer arrays')
    for name in FIGURES:
        for ext in ['png','svg']:
            if (root/OUT/(name+'.'+ext)).stat().st_size<1000:raise AssertionError('Truncated figure')
    for name in ['study_log.txt','tests_phase1.txt']:
        if 'Traceback (most recent call last)' in (root/name).read_text():
            raise AssertionError('Failed-run traceback in canonical logs')
    files=[]
    for rel in REQUIRED:
        data=(root/rel).read_bytes();files.append(dict(path=rel.as_posix(),size_bytes=len(data),sha256=hashlib.sha256(data).hexdigest()))
    return dict(schema='delivery-manifest-v1',verified=True,files=files,pytest=parse_pytest(root))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--allow-running',action='store_true');a=p.parse_args()
    result=validate(allow_running=a.allow_running)
    write_json(ROOT/'verification/artifact_manifest.json',result)
    print('VERIFIED',len(result['files']),'required files',result['pytest'])
