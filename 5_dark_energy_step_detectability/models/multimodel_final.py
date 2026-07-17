"""Final multi-model: standard, Barrow, Tsallis (power-law) + Renyi (coupled), 
vs DESI DR2 BAO + full Planck CMB. Also settle the H0 question with finer grid."""
import numpy as np, sys
from numpy.linalg import inv
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
sys.path.insert(0,'realdata')
from desi_dr2 import get_arrays
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
        Om=min(max(Y[0],1e-12),1-1e-12); return [Om*(1-Om)*(1-(n/c)*Om**(-1.0/n))]
    xb=np.linspace(0,-np.log(1+zmax),N)
    sb=solve_ivp(rhs,[0,xb[-1]],[Om_DE0/(Om_DE0+Om_m0)],t_eval=xb,rtol=1e-8,atol=1e-11)
    z=np.exp(-sb.t)-1;Om=sb.y[0];o=np.argsort(z);z=z[o];Om=Om[o]
    rm=Om_m0*(1+z)**3;rr=Om_r*(1+z)**4;rd=rm*Om/(1-Om);i0=np.argmin(np.abs(z))
    tot=rm+rr+rd;return z,np.sqrt(tot/tot[i0])

def E_renyi(d_hat,Om_m0,Om_r,zmax):
    Om_DE0=1-Om_m0-Om_r
    def integ(K):
        E_inf=np.sqrt(d_hat/(K-1))
        Nf=np.log(300.0);Nlo=np.log(1/(1+zmax))
        def rhs(N,Y):
            u=Y[0];a=np.exp(N);rd=K/(u**2*(1+d_hat*u**2))
            E=np.sqrt(max(Om_m0*a**-3+Om_r*a**-4+rd,1e-12));return [u-1/E]
        return solve_ivp(rhs,[Nf,Nlo],[1/E_inf],dense_output=True,rtol=1e-8,atol=1e-11,max_step=0.05)
    def resid(K):
        s=integ(K);u1=s.sol(0.0)[0];rd=K/(u1**2*(1+d_hat*u1**2))
        E1=np.sqrt(Om_m0+Om_r+rd);return rd/E1**2-Om_DE0
    try: K=brentq(resid,1.001,30,xtol=1e-6)
    except: return None
    s=integ(K);Ng=np.linspace(0,np.log(1/(1+zmax)),400)
    us=np.array([s.sol(N)[0] for N in Ng]);a=np.exp(Ng);z=1/a-1
    rd=K/(us**2*(1+d_hat*us**2));E=np.sqrt(Om_m0*a**-3+Om_r*a**-4+rd)
    o=np.argsort(z);return z[o],E[o]

def chi2(model,dv,Om_m0,H0,cpar,om_b=0.02237):
    h=H0/100;om_m=Om_m0*h*h;Om_r=om_r_fac*om_gamma/h**2
    if model=='standard': z,E=E_powerlaw(-2,cpar,Om_m0,Om_r,3.0)
    elif model=='barrow': z,E=E_powerlaw(dv-2,cpar,Om_m0,Om_r,3.0)
    elif model=='tsallis': z,E=E_powerlaw(2*dv-4,cpar,Om_m0,Om_r,3.0)
    elif model=='renyi':
        r=E_renyi(dv,Om_m0,Om_r,3.0)
        if r is None: return 1e9
        z,E=r
    invE=1/E;chi=np.concatenate([[0],np.cumsum(0.5*(invE[1:]+invE[:-1])*np.diff(z))])
    DM=np.interp(ZB,z,chi);DH=1/np.interp(ZB,z,E)
    o=np.empty(2*nb);o[0::2]=DM;o[1::2]=DH
    k=(o@CBI@obs)/(o@CBI@o);rr=obs-k*o;bao=rr@CBI@rr
    zst=zstar(om_b,om_m)
    if model in ('standard','barrow','tsallis'):
        nn={'standard':-2,'barrow':dv-2,'tsallis':2*dv-4}[model]
        zl,El=E_powerlaw(nn,cpar,Om_m0,Om_r,zst*1.05,2000)
    else:
        r=E_renyi(dv,Om_m0,Om_r,zst*1.05)
        if r is None: return 1e9
        zl,El=r
    m=zl<=zst;DMi=np.trapezoid(1/El[m],zl[m]);R=np.sqrt(Om_m0)*DMi
    def Emr(zz):return np.sqrt(Om_m0*(1+zz)**3+Om_r*(1+zz)**4+(1-Om_m0))
    def ig(zz):Rb=(3*om_b/(4*om_gamma))/(1+zz);return (1/np.sqrt(3*(1+Rb)))/Emr(zz)
    rs,_=quad(ig,zst,1e5,limit=200);lA=np.pi*DMi/rs
    v=np.array([R,lA,om_b]);dd=v-X_obs;return bao+dd@Cinv@dd

def best(model,dv,H0range=np.linspace(63,75,13)):
    b=1e18;bp=None
    for Om in np.linspace(0.27,0.35,13):
        for H0 in H0range:
            for c in np.linspace(0.6,1.4,9):
                v=chi2(model,dv,Om,H0,c)
                if v<b: b=v;bp=(Om,H0,c)
    return b,bp

print("FINAL MULTI-MODEL (DESI DR2 BAO + Planck CMB), best-fit total chi2:")
print(f"{'model':16}{'total':>8}{'Om':>7}{'H0':>7}{'c':>6}")
for model,dv,lab in [('standard',0,'standard HDE'),('barrow',0.1,'Barrow D=0.1'),
                     ('tsallis',0.7,'Tsallis dT=0.7'),('tsallis',0.64,'Tsallis dT=0.64'),
                     ('renyi',0.5,'Renyi d=0.5'),('renyi',1.0,'Renyi d=1.0'),('renyi',2.0,'Renyi d=2.0')]:
    b,bp=best(model,dv)
    print(f"{lab:16}{b:8.1f}{bp[0]:7.3f}{bp[1]:7.1f}{bp[2]:6.2f}")
print(f"{'LCDM ref':16}{'~10':>8}")
