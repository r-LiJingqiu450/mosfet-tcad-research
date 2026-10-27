"""Run unmodified, recovered v1 code in an isolated temporary interpreter path."""
import json,sys,tempfile,zipfile,warnings
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory() as temp:
    with zipfile.ZipFile(root/'baseline/recovered_v1_source.zip') as z:z.extractall(temp)
    sys.path.insert(0,temp)
    from semiconductor_sim.stage3_mosfet.run import default_table
    from semiconductor_sim.stage3_mosfet.gca import gca_solve,id_vg
    from semiconductor_sim.stage3_mosfet.postprocess import extract_vt,continuity_check
    t=default_table();t.build()
    vg=np.linspace(.5,2,121)
    out={'vg':vg.tolist(),'transfer_legacy':id_vg(t,vg,1.5).tolist(),
         'linear_legacy':id_vg(t,vg,.05).tolist(),'continuity':[],'ablation':[]}
    dense=np.linspace(0,2,81)
    out['vt_cc']=extract_vt(dense,id_vg(t,dense,.05))[0]
    for ny in [81,161,321,641,1281]:
        r=gca_solve(2,1.5,t,y=np.linspace(0,1e-6,ny))
        out['continuity'].append({'ny':ny,'finite_difference':continuity_check(r)[0]})
    for v,d in [(2,.05),(2,1.5),(1.175,2),(1.2,2)]:
        r=gca_solve(v,d,t)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always');r0=gca_solve(v,d,t,clm=False)
        out['ablation'].append(dict(vg=v,vd=d,legacy_ID=r['ID'],
             legacy_no_clm_ID=r0['ID'],delta_L_fraction=1-r['y_p']/1e-6,
             no_clm_profile_finite=bool(np.all(np.isfinite(r0['phi_n']))),
             warning_count=len(w)))
    print(json.dumps(out))
