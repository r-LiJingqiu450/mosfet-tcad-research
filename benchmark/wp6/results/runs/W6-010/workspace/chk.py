import devsim as ds
print(ds.symdiff(expr='diff(kahan3(a*x,b,c),x)'))
print(ds.symdiff(expr='simplify(diff(kahan3(a*x,b*x,-c*x),x))'))
import mosfet_deck as M
M.build_mesh(); M.build_physics()
for n in ['ElectronCurrent','ElectronCurrent:Electrons@n0','ElectronCurrent:Electrons@n1','ElectronCurrent:Potential@n0','ElectronCurrent:Potential@n1']:
    print(n, ds.edge_model(device='mos',region='si',name=n+'_x',equation=n))
