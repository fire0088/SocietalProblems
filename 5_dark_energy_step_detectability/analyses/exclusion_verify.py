import numpy as np, sys
from numpy.linalg import inv, cholesky
from scipy.stats import chi2 as chi2d
sys.path.insert(0,'realdata')
from desi_dr2 import get_arrays
c_km=299792.458; ZT=0.7; DZ=0.15
def Sfun(z): return 0.5*(1+np.tanh((z-ZT)/DZ))
def Efun(z,Om,A): return np.sqrt(Om*(1+z)**3+(1-Om)*(1.0+A*Sfun(z)))
ZG=np.linspace(0,3.0,900)
def comoving(Om,A):
    invE=1.0/Efun(ZG,Om,A); return np.concatenate([[0],np.cumsum(0.5*(invE[1:]+invE[:-1])*np.diff(ZG))])
ZB,DMr,DHr,CB=get_arrays(); CBI=inv(CB); LB=cholesky(CB); nb=len(ZB); H0RD=67*147.09
Om_g=np.linspace(0.20,0.42,80); A_g=np.linspace(-1.5,3.0,181)
def shape(Om,A):
    chi=comoving(Om,A); DM=np.interp(ZB,ZG,chi); DH=1.0/Efun(ZB,Om,A)
    o=np.empty(2*nb); o[0::2]=DM; o[1::2]=DH; return o
SH=np.array([[shape(Om,A) for A in A_g] for Om in Om_g]).reshape(-1,2*nb)
SHd=np.einsum('pi,pi->p',SH,(CBI@SH.T).T); nOm,nA=len(Om_g),len(A_g)
def profile_A(o):
    num=SH@(CBI@o); k=num/SHd; R=o[None,:]-k[:,None]*SH
    G=np.einsum('pi,ij,pj->p',R,CBI,R).reshape(nOm,nA)
    return G.min(0)  # profile over Om -> chi2(A)
# COVERAGE CHECK: generate mocks at A=0 truth, does 95% interval contain 0 ~95% of time?
rng=np.random.default_rng(3); truth=shape(0.33,0.0)*(c_km/H0RD)
t95=chi2d.ppf(0.95,1); cover=0; widths=[]
for _ in range(300):
    o=truth+LB@rng.standard_normal(2*nb)
    c=profile_A(o); c-=c.min(); ok=c<t95
    A_in=A_g[ok]
    if A_in.min()<=0<=A_in.max(): cover+=1
    widths.append(A_in.max()-A_in.min())
print(f"Coverage check: 95% interval contains true A=0 in {cover/300*100:.0f}% of mocks (target 95%)")
print(f"median 95% interval half-width: {np.median(widths)/2:.3f}  (real data gave +/-0.22)")
print("-> profile-likelihood interval is well-calibrated" if 90<cover/300*100<99 else "-> check coverage")
