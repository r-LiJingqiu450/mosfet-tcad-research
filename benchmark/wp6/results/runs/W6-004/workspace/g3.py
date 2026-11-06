import numpy as np, devsim as ds, mosfet_deck as M
D,S=M.DEV,M.SI
M.build_mesh(); M.build_physics()
for e in ['diff(kahan3(Electrons@n1*Bern01, Electrons@n1*vdiff, -Electrons@n0*Bern01), Potential@n0)',
          'diff(Electrons@n1*Bern01 + Electrons@n1*vdiff -Electrons@n0*Bern01, Potential@n0)',
          'diff(kahan3(Electrons@n1*Bern01, Electrons@n1*vdiff, -Electrons@n0*Bern01), Electrons@n0)',
          'diff(kahan3(Holes, -Electrons, NetDoping), Potential)']:
    print(e, '->', ds.symdiff(expr=e) if hasattr(ds,'symdiff') else None)
for nm in ['ElectronCurrent:Potential@n0','ElectronCurrent:Electrons@n0','ElectronCurrent:Electrons@n1','ElectronCurrent:Potential@n1']:
    try: print(nm, ds.get_edge_model_values(device=D, region=S, name=nm)[:3])
    except Exception as ex: print(nm, ex)
