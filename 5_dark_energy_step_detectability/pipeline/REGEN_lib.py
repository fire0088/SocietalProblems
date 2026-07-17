"""Shared clean pipeline components. Imported by per-probe runners.
No hardcoded rates; each runner validates its own null."""
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
A_G=np.round(np.arange(-0.9,3.0+1e-9,0.05),4); J0=int(np.argmin(np.abs(A_G)))
assert abs(A_G[J0])<1e-9
AMPS=[0.2,0.4,0.6,0.8,1.0]
def validate_and_rate(fits_fn,truth_fn,draw_fn,rng,nmock,name):
    o0=truth_fn(0.0)
    dn=np.clip([np.subtract(*fits_fn(draw_fn(o0))) for _ in range(nmock)],0,None)
    med=np.median(dn); fpr=np.mean(np.array(dn)>THR3)*100
    if not (med<1.0 and fpr<3.5): raise RuntimeError(f"{name} null FAIL med={med:.3f} fpr={fpr:.1f}")
    print(f"[{name}] null median={med:.3f} FPR@3s={fpr:.1f}% OK")
    out={}
    for A in AMPS:
        oA=truth_fn(A); dd=np.array([np.subtract(*fits_fn(draw_fn(oA))) for _ in range(nmock)])
        out[A]=np.mean(dd>THR3)*100
    return out
