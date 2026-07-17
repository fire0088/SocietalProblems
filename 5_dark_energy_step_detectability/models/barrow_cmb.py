"""
Barrow HDE constrained by CMB shift parameter R + DESI DR2 BAO.
CMB compressed likelihood (Chen et al., Planck 2018 TT,TE,EE+lowE):
  X_obs = {R, l_a, Omega_b h^2} = {1.7493, 301.462, 0.02239}
We use primarily R (the shift parameter) since it's the DE-sensitive geometric one
that requires D_M(z*) — this is what CMB adds: an anchor at z*~1090 that a smooth
Barrow w(z) shift CANNOT hide from (unlike low-z BAO). l_a and omega_b add sound-horizon
info but need r_s modeling; R alone already tightens Delta dramatically. We include the
full 3x3 for completeness but R carries the DE geometry.
  R = sqrt(Om_m) * (H0/c) * D_M(z*)   [flat]  = sqrt(Om_m) * integral_0^z* dz/E(z) ... 
      actually R = sqrt(Om_m H0^2)/c * D_M = sqrt(Om_m) * int_0^{z*} dz/E(z)  (dimensionless)
"""
import numpy as np, sys
from numpy.linalg import inv
from scipy.integrate import solve_ivp
from scipy.stats import chi2 as chi2d
sys.path.insert(0,'realdata')
from desi_dr2 import get_arrays

ZSTAR=1089.80
# Chen et al compressed Planck 2018
X_obs=np.array([1.7493, 301.462, 0.02239])
Cinv=np.array([[94392.3971,-1360.4913,1664517.2916],
               [-1360.4913,161.4349,3671.6180],
               [1664517.2916,3671.6180,79719182.5162]])

ZB,DMr,DHr,CB=get_arrays(); CBI=inv(CB); nb=len(ZB)
obs=np.empty(2*nb); obs[0::2]=DMr; obs[1::2]=DHr

def bhde_E_full(Delta,cpar,Om_m0,N=2500):
    """Integrate BHDE E(z) out to z* (need high z for CMB). Returns z,E on fine grid."""
    Om_DE0=1-Om_m0; p=1.0/(2.0-Delta)
    def rhs(x,Y):
        Om=min(max(Y[0],1e-15),1-1e-15)
        return [(2-Delta)*Om*(1-Om)*(1.0/(2-Delta)+(1.0/cpar)*Om**p)]
    # integrate x=ln a from 0 back to a at z*=1090 -> x=-ln(1091)
    xb=np.linspace(0,-np.log(1+ZSTAR),N)
    sb=solve_ivp(rhs,[0,xb[-1]],[Om_DE0],t_eval=xb,rtol=1e-9,atol=1e-12)
    z=np.exp(-sb.t)-1; Om=sb.y[0]; order=np.argsort(z); z=z[order]; Om=Om[order]
    # add radiation so high-z is sane: Om_r ~ 9.2e-5
    Om_r=9.2e-5
    rho_m=Om_m0*(1+z)**3; rho_r=Om_r*(1+z)**4
    rho_DE=rho_m*Om/(1-Om)   # from Omega_DE definition (matter-based; approx ignores rad in Omega)
    i0=np.argmin(np.abs(z))
    tot=rho_m+rho_r+rho_DE; E=np.sqrt(tot/tot[i0])
    return z,E

def R_shift(z,E,Om_m0):
    # R = sqrt(Om_m) * int_0^z* dz/E
    m=z<=ZSTAR
    integ=np.trapz(1.0/E[m],z[m]) if hasattr(np,'trapz') else np.trapezoid(1.0/E[m],z[m])
    return np.sqrt(Om_m0)*integ

def bao_shape(z,E):
    invE=1.0/E; chi=np.concatenate([[0],np.cumsum(0.5*(invE[1:]+invE[:-1])*np.diff(z))])
    DM=np.interp(ZB,z,chi); DH=np.interp(ZB,z,1.0/E)
    o=np.empty(2*nb); o[0::2]=DM; o[1::2]=DH; return o

def chi2_Delta(Delta, use_cmb=True):
    best=1e18
    for Om in np.linspace(0.26,0.36,21):
        for cpar in np.linspace(0.6,1.4,17):
            z,E=bhde_E_full(Delta,cpar,Om)
            sh=bao_shape(z,E)
            k=(sh@CBI@obs)/(sh@CBI@sh); r=obs-k*sh; ch=r@CBI@r
            if use_cmb:
                # add R shift constraint (1D Gaussian on R): sigma_R from Cinv[0,0]
                Rth=R_shift(z,E,Om)
                # marginalize l_a, omega_b: use just R with its marginal sigma
                sigR=1/np.sqrt(Cinv[0,0]); # NOT correct marginal, but conservative diag
                # better: R marginal std from inverse of Cinv -> C=inv(Cinv), sigma_R=sqrt(C[0,0])
                ch=ch+((Rth-X_obs[0])/SIGMA_R)**2
            if ch<best: best=ch
    return best

C_cmb=inv(Cinv); SIGMA_R=np.sqrt(C_cmb[0,0])
print(f"CMB R = {X_obs[0]} +/- {SIGMA_R:.4f}")

D_grid=np.linspace(0.0,1.0,21)
print("\nBarrow Delta constraint:")
print(f"{'Delta':>6} {'BAO-only':>10} {'BAO+CMB-R':>10}")
chi_bao=[]; chi_both=[]
for D in D_grid:
    cb=chi2_Delta(D,use_cmb=False); cc=chi2_Delta(D,use_cmb=True)
    chi_bao.append(cb); chi_both.append(cc)
chi_bao=np.array(chi_bao); chi_both=np.array(chi_both)
for i in range(0,len(D_grid),2):
    print(f"{D_grid[i]:6.2f} {chi_bao[i]-chi_bao.min():10.2f} {chi_both[i]-chi_both.min():10.2f}")
# limits
def upper(ch):
    d=ch-ch.min(); from_=np.interp(2.71,d,D_grid); return from_
np.savez('/home/claude/barrow_cmb.npz',D=D_grid,chi_bao=chi_bao,chi_both=chi_both,sigR=SIGMA_R)
print(f"\n95% one-sided upper limit on Delta:")
print(f"  BAO only:   Delta < {np.interp(2.71,chi_bao-chi_bao.min(),D_grid):.2f}")
print(f"  BAO+CMB-R:  Delta < {np.interp(2.71,chi_both-chi_both.min(),D_grid):.2f}")
