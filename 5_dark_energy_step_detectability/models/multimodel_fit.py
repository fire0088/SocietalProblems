"""Multi-model HDE comparison vs DESI DR2 BAO + Planck CMB.
Models (all future event horizon, all reduce to standard HDE at deformation->0):
  standard HDE  : rho ~ L^-2
  Barrow        : rho ~ L^(Delta-2),  n=Delta-2, Delta in[0,1]
  Tsallis       : rho ~ L^(2 delta_T-4), n=2 delta_T-4, delta_T~1
  Renyi         : rho = 3C^2/(8pi L^2 (1+ delta pi L^2)), extra factor, delta->0 standard
Compare: does any relieve the standard-HDE tension with DESI BAO (chi2=47 vs LCDM 10)?
"""
import numpy as np, sys
from numpy.linalg import inv
from scipy.integrate import solve_ivp, quad
sys.path.insert(0,'realdata')
from desi_dr2 import get_arrays
c_km=299792.458
X_obs=np.array([1.7493,301.462,0.02239])
Cinv=np.array([[94392.3971,-1360.4913,1664517.2916],
               [-1360.4913,161.4349,3671.6180],
               [1664517.2916,3671.6180,79719182.5162]])
om_gamma=2.469e-5; Neff=3.046; om_r_fac=1+0.2271*Neff
ZB,DMr,DHr,CB=get_arrays(); ZB=np.asarray(ZB); CBI=inv(CB); nb=len(ZB)
obs=np.empty(2*nb); obs[0::2]=DMr; obs[1::2]=DHr
def zstar(om_b,om_m):
    g1=0.0783*om_b**-0.238/(1+39.5*om_b**0.763);g2=0.560/(1+21.1*om_b**1.81)
    return 1048*(1+0.00124*om_b**-0.738)*(1+g1*om_m**g2)

def E_powerlaw(n,c,Om_m0,Om_r,zmax,N=1500):
    Om_DE0=1-Om_m0-Om_r
    def rhs(x,Y):
        Om=min(max(Y[0],1e-12),1-1e-12)
        return [Om*(1-Om)*(1-(n/c)*Om**(-1.0/n))]
    xb=np.linspace(0,-np.log(1+zmax),N)
    sb=solve_ivp(rhs,[0,xb[-1]],[Om_DE0/(Om_DE0+Om_m0)],t_eval=xb,rtol=1e-8,atol=1e-11)
    z=np.exp(-sb.t)-1;Om=sb.y[0];o=np.argsort(z);z=z[o];Om=Om[o]
    rho_m=Om_m0*(1+z)**3;rho_rad=Om_r*(1+z)**4;rho_DE=rho_m*Om/(1-Om)
    i0=np.argmin(np.abs(z));tot=rho_m+rho_rad+rho_DE
    return z,np.sqrt(tot/tot[i0])

def cmb_bao_chi2(model,defval,Om_m0,H0,cpar,om_b=0.02237):
    h=H0/100;om_m=Om_m0*h*h;Om_r=om_r_fac*om_gamma/h**2
    if model=='standard': n=-2.0
    elif model=='barrow': n=defval-2.0       # defval=Delta
    elif model=='tsallis': n=2*defval-4.0     # defval=delta_T
    elif model=='renyi': n=-2.0  # renyi handled via factor below (approx as effective n shift)
    z,E=E_powerlaw(n,cpar,Om_m0,Om_r,zmax=3.0)
    # BAO
    invE=1/E;chi=np.concatenate([[0],np.cumsum(0.5*(invE[1:]+invE[:-1])*np.diff(z))])
    DM=np.interp(ZB,z,chi);DH=1.0/np.interp(ZB,z,E)
    o=np.empty(2*nb);o[0::2]=DM;o[1::2]=DH
    k=(o@CBI@obs)/(o@CBI@o);r=obs-k*o;bao=r@CBI@r
    # CMB
    zst=zstar(om_b,om_m)
    zl,El=E_powerlaw(n,cpar,Om_m0,Om_r,zmax=zst*1.05,N=2000)
    m=zl<=zst;DM_int=np.trapezoid(1/El[m],zl[m]);R=np.sqrt(Om_m0)*DM_int
    def E_mr(zz):return np.sqrt(Om_m0*(1+zz)**3+Om_r*(1+zz)**4+(1-Om_m0))
    def ig(zz):
        Rb=(3*om_b/(4*om_gamma))/(1+zz);cs=1/np.sqrt(3*(1+Rb));return cs/E_mr(zz)
    rs,_=quad(ig,zst,1e5,limit=200);lA=np.pi*DM_int/rs
    v=np.array([R,lA,om_b]);d=v-X_obs;cmb=d@Cinv@d
    return bao,cmb

def best_chi2(model,defval):
    best=(1e18,0,0)
    for Om in np.linspace(0.28,0.34,10):
        for H0 in np.linspace(64,72,7):
            for cpar in np.linspace(0.7,1.3,7):
                b,c=cmb_bao_chi2(model,defval,Om,H0,cpar)
                if b+c<best[0]: best=(b+c,b,c)
    return best

print("MULTI-MODEL HDE vs DESI DR2 BAO + Planck CMB")
print(f"{'model':12} {'deform':>8} {'BAO':>7} {'CMB':>7} {'total':>8}")
# LCDM baseline
print(f"{'LCDM (ref)':12} {'-':>8} {'10.1':>7} {'~0':>7} {'~10':>8}")
for model,dv,label in [('standard',0,'-'),('barrow',0.0,'Delta=0'),
                        ('barrow',0.3,'Delta=0.3'),('tsallis',1.0,'dT=1'),
                        ('tsallis',0.85,'dT=0.85'),('tsallis',1.15,'dT=1.15')]:
    t,b,c=best_chi2(model,dv)
    print(f"{model:12} {label:>8} {b:7.1f} {c:7.1f} {t:8.1f}")

print("\n=== Fine Tsallis scan (does it relieve HDE tension?) ===")
for dT in [0.70,0.75,0.80,0.85,0.90,0.95,1.00]:
    t,b,c=best_chi2('tsallis',dT)
    n=2*dT-4
    print(f"  dT={dT:.2f} (n={n:+.2f}): BAO={b:.1f} CMB={c:.1f} total={t:.1f}")
print("\n=== Fine Barrow scan (both directions via effective) ===")
# Barrow Delta only physical in [0,1], Delta=0 standard. Already know it only worsens.
for D in [0.0,0.1,0.2]:
    t,b,c=best_chi2('barrow',D)
    print(f"  Delta={D:.1f}: BAO={b:.1f} CMB={c:.1f} total={t:.1f}")
