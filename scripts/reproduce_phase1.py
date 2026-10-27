"""Complete sequential reproduction. Never freeze/regenerate regression baselines.

Exit code and JUnit must agree; reports and delivery checks run only after pytest
finishes. An interrupted run remains running/failed and cannot be released.
"""
import argparse,hashlib,os,subprocess,sys,time,traceback
from pathlib import Path
from datetime import datetime,timezone,timedelta

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from semiconductor_sim.common.jsonio import write_json


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--fresh',action='store_true')
    args=parser.parse_args()
    if args.fresh:
        import shutil
        shutil.rmtree(ROOT/'.cache',ignore_errors=True)
    directory=ROOT/'verification';directory.mkdir(exist_ok=True)
    status=dict(state='running',started_at=datetime.now(timezone(timedelta(hours=8))).isoformat(),
                exit_code=None,steps=[])
    state_path=directory/'run_status.json';write_json(state_path,status)
    golden=ROOT/'semiconductor_sim/results/golden/stage3.json'
    if not golden.is_file():
        status.update(state='failed',exit_code=2,error='Missing production golden')
        write_json(state_path,status)
        raise SystemExit('Production golden is required; reproduction will not create it')
    before=hashlib.sha256(golden.read_bytes()).hexdigest()
    env=os.environ.copy();env.update(PYTHONUNBUFFERED='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',MPLBACKEND='Agg')
    def run(name,command,log):
        print(name,flush=True);t=time.perf_counter()
        step=dict(name=name,state='running',command=command)
        status['steps'].append(step);write_json(state_path,status)
        with (ROOT/log).open('w',encoding='utf-8') as f:
            result=subprocess.run(command,cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT)
        step.update(state='completed' if result.returncode==0 else 'failed',exit_code=result.returncode,
                    elapsed_seconds=time.perf_counter()-t,log=log)
        write_json(state_path,status)
        if result.returncode:raise RuntimeError(f'{name} returned {result.returncode}; see {log}')
    try:
        run('1. Numerical study and attribution',[sys.executable,'-m','semiconductor_sim.stage3_mosfet.phase1'],'study_log.txt')
        run('2. Full pytest suite',[sys.executable,'-m','pytest','semiconductor_sim','-vv','-rA',
            '--junitxml=verification/pytest.xml'],'tests_phase1.txt')
        from verify_delivery import parse_pytest
        pytest_status=parse_pytest(ROOT)
        write_json(directory/'pytest_status.json',pytest_status)
        after=hashlib.sha256(golden.read_bytes()).hexdigest()
        if before!=after:raise RuntimeError('Frozen golden changed during normal reproduction')
        status['golden_sha256_before']=before;status['golden_sha256_after']=after
        status['pytest']=pytest_status
        write_json(state_path,status)
        run('3. Generate reports',[sys.executable,'scripts/write_phase1_reports.py'],'verification/report_generation.txt')
        run('4. Verify all deliverables',[sys.executable,'scripts/verify_delivery.py','--allow-running'],'verification/artifact_check.txt')
        status.update(state='completed',exit_code=0,
            completed_at=datetime.now(timezone(timedelta(hours=8))).isoformat())
        write_json(state_path,status)
        print('COMPLETED: '+str(pytest_status),flush=True)
    except BaseException as exc:
        status.update(state='failed',exit_code=1,error=str(exc));write_json(state_path,status)
        raise

if __name__=='__main__':main()
