"""Explicit, audited production baseline freeze. Never run by reproduction/tests."""
import argparse,sys
from pathlib import Path
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root))
from semiconductor_sim.common.jsonio import write_json
from semiconductor_sim.stage3_mosfet.golden import GOLDEN_PATH,production_numbers
from semiconductor_sim.stage3_mosfet.config import production_metadata
from semiconductor_sim.stage3_mosfet.phase1 import production_table
from semiconductor_sim.stage3_mosfet.reference import reference_current
from semiconductor_sim.stage3_mosfet.gca import gca_solve

def main():
    p=argparse.ArgumentParser();p.add_argument('--reason',required=True)
    p.add_argument('--allow-replace',action='store_true');args=p.parse_args()
    if GOLDEN_PATH.exists() and not args.allow_replace:
        raise SystemExit('Refusing replacement without --allow-replace')
    table,_=production_table()
    verification=[]
    for vg,vd in [(.5,.05),(1,.05),(1.25,.2),(2,1.5)]:
        ref=reference_current(vg,vd)['ID'];err=gca_solve(vg,vd,table)['ID']/ref-1
        if abs(err)>=.005:raise RuntimeError('Cannot freeze: independent-reference target failed')
        verification.append(dict(vg=vg,vd=vd,relative_error=err))
    write_json(GOLDEN_PATH,dict(version='phase1-1.2',reason=args.reason,
        config=production_metadata(),numbers=production_numbers(table),
        numerical_verification_before_freeze=verification,
        evidence_type='new_revision_regression_not_original_baseline'))
    print(GOLDEN_PATH)

if __name__=='__main__':main()
