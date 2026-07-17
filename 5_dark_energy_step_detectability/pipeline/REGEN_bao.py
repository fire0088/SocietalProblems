from REGEN_lib import *
from desi_dr2 import get_arrays
import numpy as np
from numpy.linalg import inv,cholesky
rng=np.random.default_rng(3)
ZB,DMr,DHr,CB=get_arrays(); CI=inv(CB); L=cholesky(CB); H0RD=67.0*147.09
def shape(Om,A):
    chi=comoving(Om,A); DM=np.interp(ZB,ZG,chi); DH=1.0/Efun(ZB,Om,A)
    o=np.empty(2*len(ZB)); o[0::2]=DM; o[1::2]=DH; return o
SH=np.array([[shape(Om,A) for A in A_G] for Om in OM_G]).reshape(-1,2*len(ZB))
def fits(o):
    num=SH@(CI@o); denom=np.einsum('pi,pi->p',SH,(CI@SH.T).T); k=num/denom
    R=o[None,:]-k[:,None]*SH; G=np.einsum('pi,ij,pj->p',R,CI,R).reshape(len(OM_G),len(A_G))
    return G[:,J0].min(),G.min()
def truth(A): return shape(0.33,A)*(c_km/H0RD)
def draw(o): return o+L@rng.standard_normal(2*len(ZB))
r=validate_and_rate(fits,truth,draw,rng,400,'BAO')
np.savez('/home/claude/r_bao.npz',**{str(k):v for k,v in r.items()})
print("BAO:",{k:round(v) for k,v in r.items()})
