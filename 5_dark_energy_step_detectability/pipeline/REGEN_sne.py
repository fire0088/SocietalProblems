from REGEN_lib import *
import numpy as np
from numpy.linalg import inv,cholesky
rng=np.random.default_rng(2)
raw=np.fromfile('realdata/PPcov.cov',sep=' '); N=int(raw[0]); Cf=raw[1:].reshape(N,N); Cf=0.5*(Cf+Cf.T)
with open('realdata/PantheonSH0ES.dat') as f:
    h=f.readline().split(); ix={n:i for i,n in enumerate(h)}; R=[l.split() for l in f]
z=np.array([float(r[ix['zHD']]) for r in R]); cal=np.array([int(r[ix['IS_CALIBRATOR']]) for r in R])
m=(z>0.01)&(cal==0); ZP=z[m]; CPP=Cf[np.ix_(m,m)]
CI=inv(CPP); L=cholesky(CPP); ONE=np.ones(len(ZP)); A1=ONE@CI@ONE
def mu(Om,A): return 5*np.log10((1+ZP)*c_km*np.interp(ZP,ZG,comoving(Om,A)))+25
MU=np.array([[mu(Om,A) for A in A_G] for Om in OM_G]).reshape(-1,len(ZP))
MUI=(CI@MU.T).T; selfI=np.einsum('pi,pi->p',MU,MUI); MU1=MU@(CI@ONE)
def fits(o):
    oI=CI@o; rCr=(o@oI)-2*MU@oI+selfI; B=ONE@oI-MU1
    G=(rCr-B*B/A1).reshape(len(OM_G),len(A_G)); return G[:,J0].min(),G.min()
def truth(A): return mu(0.33,A)
def draw(o): return o+L@rng.standard_normal(len(ZP))
r=validate_and_rate(fits,truth,draw,rng,400,'SNe')
np.savez('/home/claude/r_sne.npz',**{str(k):v for k,v in r.items()})
print("SNe:",{k:round(v) for k,v in r.items()})
