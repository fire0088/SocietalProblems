from REGEN_lib import *
from cc_data import Z as ZC, covariance
import numpy as np
from numpy.linalg import inv,cholesky
C=covariance('ooo'); CI=inv(C); L=cholesky(C); NZ=len(ZC); rng=np.random.default_rng(1)
OMc=np.linspace(0.15,0.5,22); H0c=np.linspace(56,80,22)
BASE=(H0c[:,None,None,None]*np.sqrt(OMc[None,:,None,None]*(1+ZC)**3
      +(1-OMc[None,:,None,None])*(1.0+A_G[None,None,:,None]*Sfun(ZC)[None,None,None,:]))).reshape(-1,len(A_G),NZ)
def fits(d):
    R=d[None,None,:]-BASE; chis=np.einsum('pai,ij,paj->pa',R,CI,R); return chis[:,J0].min(),chis.min()
def truth(A): return 67.0*np.sqrt(0.33*(1+ZC)**3+0.67*(1.0+A*Sfun(ZC)))
def draw(o): return o+L@rng.standard_normal(NZ)
r=validate_and_rate(fits,truth,draw,rng,150,'CC')
np.savez('/home/claude/r_cc.npz',**{str(k):v for k,v in r.items()})
print("CC:",{k:round(v) for k,v in r.items()})
