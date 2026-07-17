"""
#4 EXCLUSION: what step amplitudes does REAL DESI DR2 already rule out?
Fit the REAL BAO data with the step model, profile out (Om, scale), get the
profile-likelihood chi2(A). The 95%/99% CL exclusion region is where
chi2(A) - chi2_min > threshold. This turns detectability -> a real measurement:
an upper limit on how big a sudden DE transition can be, given today's data.
Fixed location zt=0.7, dz=0.15 (same as detection study).
"""
import numpy as np, sys
from numpy.linalg import inv
from scipy.stats import chi2 as chi2d, norm
sys.path.insert(0,'realdata')
from desi_dr2 import get_arrays
c_km=299792.458; ZT=0.7; DZ=0.15
def Sfun(z): return 0.5*(1+np.tanh((z-ZT)/DZ))
def Efun(z,Om,A): return np.sqrt(Om*(1+z)**3+(1-Om)*(1.0+A*Sfun(z)))
ZG=np.linspace(0,3.0,900)
def comoving(Om,A):
    invE=1.0/Efun(ZG,Om,A); return np.concatenate([[0],np.cumsum(0.5*(invE[1:]+invE[:-1])*np.diff(ZG))])

ZB,DMr,DHr,CB=get_arrays(); CBI=inv(CB); nb=len(ZB)
obs=np.empty(2*nb); obs[0::2]=DMr; obs[1::2]=DHr

def bao_shape(Om,A):
    chi=comoving(Om,A); DM=np.interp(ZB,ZG,chi); DH=1.0/Efun(ZB,Om,A)
    o=np.empty(2*nb); o[0::2]=DM; o[1::2]=DH; return o

# fine A grid, fine Om grid; profile scale analytically, Om numerically
Om_g=np.linspace(0.20,0.42,120)
A_g=np.linspace(-1.5,3.0,451)
def chi2_of_A(A):
    best=1e18
    for Om in Om_g:
        sh=bao_shape(Om,A)
        k=(sh@CBI@obs)/(sh@CBI@sh)   # profile scale
        r=obs-k*sh
        c=r@CBI@r
        if c<best: best=c
    return best

chi2A=np.array([chi2_of_A(A) for A in A_g])
imin=chi2A.argmin(); A_best=A_g[imin]; chi2min=chi2A[imin]
dchi=chi2A-chi2min
# CL thresholds for 1 param
t68,t95,t99=chi2d.ppf(0.68,1),chi2d.ppf(0.95,1),chi2d.ppf(0.997,1)
# find where dchi crosses each threshold around the minimum
def interval(thr):
    below=dchi<thr
    idx=np.where(below)[0]
    return A_g[idx[0]],A_g[idx[-1]]
lo68,hi68=interval(t68); lo95,hi95=interval(t95); lo99,hi99=interval(t99)
print(f"REAL DESI DR2 step-amplitude fit (fixed zt=0.7, dz=0.15):")
print(f"  best-fit A = {A_best:+.3f}  (A=0 is LCDM/no-step)")
print(f"  chi2_min = {chi2min:.2f}, chi2(A=0) = {chi2A[np.argmin(np.abs(A_g))]:.2f}, Delta = {chi2A[np.argmin(np.abs(A_g))]-chi2min:.2f}")
print(f"  68% CL:  A in [{lo68:+.2f}, {hi68:+.2f}]")
print(f"  95% CL:  A in [{lo95:+.2f}, {hi95:+.2f}]")
print(f"  99.7% CL: A in [{lo99:+.2f}, {hi99:+.2f}]")
# convert to "% density shed/gained" at 95%
def frac(A): return A/(1+A)*100
print(f"\n  => 95% CL EXCLUDES thinning beyond {frac(hi95):.0f}% shed and thickening beyond {frac(lo95):.0f}% (A={lo95:.2f})")
print(f"     i.e. a sudden DE transition at z=0.7 larger than these is ruled out by DESI DR2 alone")
np.savez('/home/claude/exclusion.npz',A_g=A_g,dchi=dchi,A_best=A_best,
    lo95=lo95,hi95=hi95,lo99=lo99,hi99=hi99,lo68=lo68,hi68=hi68)
print("saved")
