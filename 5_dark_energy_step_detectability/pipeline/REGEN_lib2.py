"""Extended pipeline: signed amplitudes (thinning A>0 and thickening A<0).
Same validated grid method. A/(1+A) = fraction shed (A>0) or gained (A<0)."""
import numpy as np, sys
from numpy.linalg import inv, cholesky
from scipy.stats import chi2 as chi2d
sys.path.insert(0,'realdata')
c_km=299792.458; ZT=0.7; DZ=0.15
THR3=chi2d.ppf(0.9973,1); THR2=chi2d.ppf(0.9545,1)
def Sfun(z): return 0.5*(1+np.tanh((z-ZT)/DZ))
def Efun(z,Om,A): return np.sqrt(Om*(1+z)**3+(1-Om)*(1.0+A*Sfun(z)))
ZG=np.linspace(0,3.0,900)
def comoving(Om,A):
    invE=1.0/Efun(ZG,Om,A); return np.concatenate([[0],np.cumsum(0.5*(invE[1:]+invE[:-1])*np.diff(ZG))])
OM_G=np.linspace(0.15,0.45,41)
# amplitude grid must include A=0 AND span negative -> positive
A_G=np.round(np.arange(-0.9,3.0+1e-9,0.05),4); J0=int(np.argmin(np.abs(A_G)))
assert abs(A_G[J0])<1e-9
# signed test amplitudes: thinning (+) and thickening (-)
# For A<0, 1+A must stay >0 for rho_DE>0; A>-1 always. Use symmetric magnitudes.
SIGNED_AMPS=[-0.6,-0.4,-0.2,0.2,0.4,0.6]   # |frac| ~ 0.17,0.29,0.375,0.29,0.17 ... note asymmetric mapping
def frac(A): return A/(1+A)*100  # signed % change in density
def validate_and_rate(fits_fn,truth_fn,draw_fn,rng,nmock,name,amps):
    o0=truth_fn(0.0)
    dn=np.clip([np.subtract(*fits_fn(draw_fn(o0))) for _ in range(nmock)],0,None)
    med=np.median(dn); fpr=np.mean(np.array(dn)>THR3)*100
    if not (med<1.1 and fpr<3.5): raise RuntimeError(f"{name} null FAIL med={med:.3f}")
    print(f"[{name}] null median={med:.3f} FPR@3s={fpr:.1f}% OK")
    out={}
    for A in amps:
        oA=truth_fn(A); dd=np.array([np.subtract(*fits_fn(draw_fn(oA))) for _ in range(nmock)])
        out[A]=np.mean(dd>THR3)*100
    return out
