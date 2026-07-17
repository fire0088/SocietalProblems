"""
#3 LOCATION SCAN: detection rate vs transition redshift z_t.
Everything so far fixed z_t=0.7. Slide it across 0.3-1.5 and map where each
probe is sensitive. Key question: is BAO's advantage uniform, or concentrated
near its tracer redshifts (0.51-2.33)? Fixed amplitude A=0.6 (37.5% shed),
dz=0.15. Validated grid method, real covariances, null-checked per z_t.
"""
import numpy as np, sys
from numpy.linalg import inv, cholesky
from scipy.stats import chi2 as chi2d
sys.path.insert(0,'realdata')
from desi_dr2 import get_arrays
c_km=299792.458; DZ=0.15; THR3=chi2d.ppf(0.9973,1); A_TEST=0.6
def Sfun(z,zt): return 0.5*(1+np.tanh((z-zt)/DZ))
def Efun(z,Om,A,zt): return np.sqrt(Om*(1+z)**3+(1-Om)*(1.0+A*Sfun(z,zt)))
ZG=np.linspace(0,3.0,900)
def comoving(Om,A,zt):
    invE=1.0/Efun(ZG,Om,A,zt); return np.concatenate([[0],np.cumsum(0.5*(invE[1:]+invE[:-1])*np.diff(ZG))])
OM_G=np.linspace(0.15,0.45,41); A_G=np.round(np.arange(-0.9,3.0+1e-9,0.05),4); J0=int(np.argmin(np.abs(A_G)))
rng=np.random.default_rng(11); NM=300; H0RD=67*147.09

# BAO
ZB,DMr,DHr,CB=get_arrays(); CBI=inv(CB); LB=cholesky(CB); nb=len(ZB)
def bao_shape(Om,A,zt):
    chi=comoving(Om,A,zt); DM=np.interp(ZB,ZG,chi); DH=1.0/Efun(ZB,Om,A,zt)
    o=np.empty(2*nb); o[0::2]=DM; o[1::2]=DH; return o
def bao_rate(zt):
    SH=np.array([[bao_shape(Om,A,zt) for A in A_G] for Om in OM_G]).reshape(-1,2*nb)
    SHd=np.einsum('pi,pi->p',SH,(CBI@SH.T).T)
    def fits(o):
        k=(SH@(CBI@o))/SHd; R=o[None,:]-k[:,None]*SH
        G=np.einsum('pi,ij,pj->p',R,CBI,R).reshape(len(OM_G),len(A_G)); return G[:,J0].min(),G.min()
    o0=bao_shape(0.33,0.0,zt)*(c_km/H0RD)
    dn=np.clip([np.subtract(*fits(o0+LB@rng.standard_normal(2*nb))) for _ in range(NM)],0,None)
    if np.median(dn)>1.3: return np.nan
    oa=bao_shape(0.33,A_TEST,zt)*(c_km/H0RD)
    dd=np.array([np.subtract(*fits(oa+LB@rng.standard_normal(2*nb))) for _ in range(NM)])
    return np.mean(dd>THR3)*100

# SNe (real Pantheon+)
raw=np.fromfile('realdata/PPcov.cov',sep=' '); N=int(raw[0]); Cf=raw[1:].reshape(N,N); Cf=0.5*(Cf+Cf.T)
with open('realdata/PantheonSH0ES.dat') as f:
    h=f.readline().split(); ix={n:i for i,n in enumerate(h)}; R=[l.split() for l in f]
z=np.array([float(r[ix['zHD']]) for r in R]); cal=np.array([int(r[ix['IS_CALIBRATOR']]) for r in R])
m=(z>0.01)&(cal==0); ZP=z[m]; CPP=Cf[np.ix_(m,m)]; CPPI=inv(CPP); LP=cholesky(CPP)
ONE=np.ones(len(ZP)); A1P=ONE@CPPI@ONE; CppI1=CPPI@ONE
def sn_mu(Om,A,zt): return 5*np.log10((1+ZP)*c_km*np.interp(ZP,ZG,comoving(Om,A,zt)))+25
def sn_rate(zt):
    MU=np.array([[sn_mu(Om,A,zt) for A in A_G] for Om in OM_G]).reshape(-1,len(ZP))
    MUI=(CPPI@MU.T).T; selfI=np.einsum('pi,pi->p',MU,MUI); MU1=MU@CppI1
    def fits(o):
        oI=CPPI@o; rCr=(o@oI)-2*MU@oI+selfI; B=(ONE@oI)-MU1
        G=(rCr-B*B/A1P).reshape(len(OM_G),len(A_G)); return G[:,J0].min(),G.min()
    o0=sn_mu(0.33,0.0,zt)
    dn=np.clip([np.subtract(*fits(o0+LP@rng.standard_normal(len(ZP)))) for _ in range(NM)],0,None)
    if np.median(dn)>1.3: return np.nan
    oa=sn_mu(0.33,A_TEST,zt)
    dd=np.array([np.subtract(*fits(oa+LP@rng.standard_normal(len(ZP)))) for _ in range(NM)])
    return np.mean(dd>THR3)*100

zt_scan=np.array([0.3,0.5,0.7,0.9,1.1,1.3,1.5])
print("LOCATION SCAN (A=0.6/37.5% shed, vary transition redshift zt)")
print(f"{'zt':>5} {'BAO':>6} {'SNe':>6}")
baoZ=[]; sneZ=[]
for zt in zt_scan:
    b=bao_rate(zt); baoZ.append(b)
    print(f"{zt:5.1f} {b:6.0f}", end="", flush=True); 
    s=sn_rate(zt); sneZ.append(s)
    print(f" {s:6.0f}")
np.savez('/home/claude/location.npz',zt=zt_scan,bao=baoZ,sne=sneZ,
    tracers=ZB)
print("\nDESI tracer redshifts:",np.round(ZB,2))
