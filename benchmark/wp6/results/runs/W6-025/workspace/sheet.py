import numpy as np, json
m=json.load(open('mesh.json')); st=np.load('zero_bias_state.npz'); print(list(st.keys()))
xs=np.array(m['x_cm']); ysi=np.array(m['y_si_cm']); nx=len(xs); ny=len(ysi)
n=st['n_si'].reshape(ny,nx)  # assume row-major j*nx+i
q=1.602176634e-19; mu=1350
# trapezoid in y
Ns=np.trapezoid(n,ysi,axis=0)
R=np.trapezoid(1/(q*mu*Ns),xs)
print('G_1D',1/R)
sel=(xs>=0)&(xs<=1e-4); print('R channel part',np.trapezoid(1/(q*mu*Ns[sel]),xs[sel]), 'R total',R)
print('Ns mid', Ns[nx//2], 'n surface mid', n[0,nx//2])
