"""
Full 3-param CMB distance priors for Barrow HDE — corrected sound horizon.
Key fix: above z*, dark energy is utterly negligible, so r_s uses the standard
matter+radiation E(z) analytically (integrate z* to ~1e5). Below z*, D_M uses
the full BHDE E(z). This is standard practice and avoids integrating the DE ODE
to absurd redshift.
"""
import numpy as np, sys
from numpy.linalg import inv
from scipy.integrate import solve_ivp, quad
from scipy.stats import chi2 as chi2d
sys.path.insert(0,'realdata')
from desi_dr2 import get_arrays
c_km=299792.458
X_obs=np.array([1.7493, 301.462, 0.02239])
Cinv=np.array([[94392.3971,-1360.4913,1664517.2916],
               [-1360.4913,161.4349,3671.6180],
               [1664517.2916,3671.6180,79719182.5162]])
om_gamma=2.469e-5; Neff=3.046; om_r_fac=1+0.2271*Neff

def zstar(om_b,om_m):
    g1=0.0783*om_b**-0.238/(1+39.5*om_b**0.763)
    g2=0.560/(1+21.1*om_b**1.81)
    return 1048*(1+0.00124*om_b**-0.738)*(1+g1*om_m**g2)

def bhde_E_lowz(Delta,cpar,Om_m0,Om_r,zmax,N=2500):
    Om_DE0=1-Om_m0; p=1.0/(2.0-Delta)
    def rhs(x,Y):
        Om=min(max(Y[0],1e-15),1-1e-15)
        return [(2-Delta)*Om*(1-Om)*(1.0/(2-Delta)+(1.0/cpar)*Om**p)]
    xb=np.linspace(0,-np.log(1+zmax),N)
    sb=solve_ivp(rhs,[0,xb[-1]],[Om_DE0],t_eval=xb,rtol=1e-9,atol=1e-12)
    z=np.exp(-sb.t)-1; Om=sb.y[0]; o=np.argsort(z); z=z[o]; Om=Om[o]
    rho_m=Om_m0*(1+z)**3; rho_rad=Om_r*(1+z)**4; rho_DE=rho_m*Om/(1-Om)
    i0=np.argmin(np.abs(z)); tot=rho_m+rho_rad+rho_DE
    return z, np.sqrt(tot/tot[i0])

def cmb_vector(Delta,cpar,Om_m0,H0,om_b):
    h=H0/100.0; om_m=Om_m0*h*h
    Om_r=om_r_fac*om_gamma/(h*h)
    zst=zstar(om_b,om_m)
    # D_M(z*): integrate BHDE E from 0 to z*
    z,E=bhde_E_lowz(Delta,cpar,Om_m0,Om_r,zmax=zst*1.05)
    m=z<=zst; DM_int=np.trapezoid(1.0/E[m],z[m])
    R=np.sqrt(Om_m0)*DM_int
    # r_s(z*): DE negligible above z*, use matter+radiation E analytically
    def E_mr(zz): return np.sqrt(Om_m0*(1+zz)**3+Om_r*(1+zz)**4+(1-Om_m0))  # DE~const, tiny at high z
    def integrand(zz):
        Rb=(3*om_b/(4*om_gamma))/(1+zz)
        cs=1.0/np.sqrt(3*(1+Rb))
        return cs/E_mr(zz)
    rs_dimless,_=quad(integrand, zst, 1e5, limit=200)
    lA=np.pi*DM_int/rs_dimless
    return np.array([R,lA,om_b]), zst

for (Om,H0,omb) in [(0.31,67.5,0.02237),(0.30,68.0,0.02240)]:
    v,zst=cmb_vector(0.0,1.0,Om,H0,omb)
    d=v-X_obs; chi2=d@Cinv@d
    print(f"Om={Om} H0={H0} omb={omb}: R={v[0]:.4f} l_A={v[1]:.2f} z*={zst:.1f} chi2={chi2:.1f}")
print("(targets: R=1.7493, l_A=301.46)")

print("\n=== FULL CMB+BAO Barrow constraint (profile Om,H0,om_b,c,scale) ===")
ZB,DMr,DHr,CB=get_arrays(); CBI=inv(CB); nb=len(ZB)
obs=np.empty(2*nb); obs[0::2]=DMr; obs[1::2]=DHr

def bao_chi2(Delta,cpar,Om_m0,H0,om_b):
    h=H0/100; Om_r=om_r_fac*om_gamma/h**2
    z,E=bhde_E_lowz(Delta,cpar,Om_m0,Om_r,zmax=3.0)
    invE=1/E; chi=np.concatenate([[0],np.cumsum(0.5*(invE[1:]+invE[:-1])*np.diff(z))])
    DM=np.interp(ZB,z,chi); DH=np.interp(ZB,z,1/E)
    o=np.empty(2*nb); o[0::2]=DM; o[1::2]=DH
    k=(o@CBI@obs)/(o@CBI@o); r=obs-k*o
    return r@CBI@r

def total_chi2(Delta):
    best=1e18
    for Om in np.linspace(0.28,0.34,13):
        for H0 in np.linspace(64,72,9):
            for cpar in np.linspace(0.7,1.3,7):
                om_b=0.02237  # tightly constrained, fix near Planck
                v,_=cmb_vector(Delta,cpar,Om,H0,om_b)
                d=v-X_obs; cmb=d@Cinv@d
                bao=bao_chi2(Delta,cpar,Om,H0,om_b)
                t=cmb+bao
                if t<best: best=t
    return best

D_grid=np.array([0.0,0.05,0.1,0.2,0.3,0.5])
ch=np.array([total_chi2(D) for D in D_grid])
dch=ch-ch.min()
print("Delta:", D_grid)
print("dchi2:", np.round(dch,2))
from numpy import interp
print(f"95% one-sided upper limit: Delta < {interp(2.71,dch,D_grid):.3f}")
np.savez('/home/claude/barrow_cmb_full.npz',D=D_grid,dchi=dch,chimin=ch.min())
