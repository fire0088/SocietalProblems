import numpy as np, sys
from numpy.linalg import inv
from scipy.integrate import solve_ivp
from scipy.stats import chi2 as chi2d
sys.path.insert(0,'realdata')
from desi_dr2 import get_arrays
ZB,DMr,DHr,CB=get_arrays(); CBI=inv(CB); nb=len(ZB)
obs=np.empty(2*nb); obs[0::2]=DMr; obs[1::2]=DHr

def bhde_E(Delta,cpar,Om_m0,zmax=3.0,N=1200):
    Om_DE0=1-Om_m0; p=1.0/(2.0-Delta)
    def rhs(x,Y):
        Om=min(max(Y[0],1e-12),1-1e-12)
        return [(2-Delta)*Om*(1-Om)*(1.0/(2-Delta)+(1.0/cpar)*Om**p)]
    xb=np.linspace(0,-np.log(1+zmax),N)
    sb=solve_ivp(rhs,[0,xb[-1]],[Om_DE0],t_eval=xb,rtol=1e-8,atol=1e-10)
    z=np.exp(-sb.t)-1; Om=sb.y[0]; order=np.argsort(z); z=z[order]; Om=Om[order]
    rho_m=Om_m0*(1+z)**3; rho_DE=rho_m*Om/(1-Om); i0=np.argmin(np.abs(z))
    E=np.sqrt((rho_m+rho_DE)/(rho_m[i0]+rho_DE[i0])); return z,E

def shape(zg,E):
    invE=1.0/E; chi=np.concatenate([[0],np.cumsum(0.5*(invE[1:]+invE[:-1])*np.diff(zg))])
    DM=np.interp(ZB,zg,chi); DH=np.interp(ZB,zg,1.0/E)
    o=np.empty(2*nb); o[0::2]=DM; o[1::2]=DH; return o

def chi2_Delta(Delta):
    best=1e18
    for Om in np.linspace(0.24,0.38,29):
        for cpar in np.linspace(0.55,1.45,19):
            zg,E=bhde_E(Delta,cpar,Om); sh=shape(zg,E)
            k=(sh@CBI@obs)/(sh@CBI@sh); r=obs-k*sh; ch=r@CBI@r
            if ch<best: best=ch
    return best

D_grid=np.linspace(0.0,1.0,21)
chi2D=np.array([chi2_Delta(D) for D in D_grid])
imin=chi2D.argmin(); D_best=D_grid[imin]; chi2min=chi2D[imin]; dchi=chi2D-chi2min
t95=chi2d.ppf(0.95,1); t99=chi2d.ppf(0.997,1)
b95=D_grid[dchi<t95]; b99=D_grid[dchi<t99]
print("Barrow HDE (future event horizon) vs REAL DESI DR2 BAO")
print(f"  chi2_min = {chi2min:.2f}  at Delta = {D_best:.3f}")
print(f"  chi2(Delta=0, standard HDE) = {chi2D[0]:.2f}")
print(f"  Delta=0 vs best Delta-chi2 = {chi2D[0]-chi2min:.2f}")
print(f"  95% CL: Delta in [{b95.min():.2f},{b95.max():.2f}]")
print(f"  99.7% CL: Delta in [{b99.min():.2f},{b99.max():.2f}]")
if b95.max()<1.0: print(f"  => DESI DR2 EXCLUDES Delta > {b95.max():.2f} at 95%")
else: print(f"  => Delta not excluded up to 1.0")
# also report chi2 vs Delta curve
print("  Delta:",np.round(D_grid[::4],2))
print("  dchi2:",np.round(dchi[::4],2))
np.savez('/home/claude/barrow_exclusion.npz',D=D_grid,dchi=dchi,D_best=D_best,chi2D=chi2D,chi2min=chi2min)
print("saved")
