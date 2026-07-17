import numpy as np, sys
from numpy.linalg import inv, cholesky
from scipy.stats import chi2 as chi2d
sys.path.insert(0,'realdata')
from desi_dr2 import get_arrays
c_km=299792.458; ZT=0.7; DZ=0.15; THR3=chi2d.ppf(0.9973,1)
def Sfun(z): return 0.5*(1+np.tanh((z-ZT)/DZ))
def Efun(z,Om,A): return np.sqrt(Om*(1+z)**3+(1-Om)*(1.0+A*Sfun(z)))
ZG=np.linspace(0,3.0,900)
def comoving(Om,A):
    invE=1.0/Efun(ZG,Om,A); return np.concatenate([[0],np.cumsum(0.5*(invE[1:]+invE[:-1])*np.diff(ZG))])
OM_G=np.linspace(0.15,0.45,41); A_G=np.round(np.arange(-0.9,3.0+1e-9,0.05),4); J0=int(np.argmin(np.abs(A_G)))
rng=np.random.default_rng(7); NM=300; H0RD=67.0*147.09
nOm,nA=len(OM_G),len(A_G)

# SNe setup - precompute CppI @ MU for all templates ONCE
raw=np.fromfile('realdata/PPcov.cov',sep=' '); N=int(raw[0]); Cf=raw[1:].reshape(N,N); Cf=0.5*(Cf+Cf.T)
with open('realdata/PantheonSH0ES.dat') as f:
    h=f.readline().split(); ix={n:i for i,n in enumerate(h)}; R=[l.split() for l in f]
z=np.array([float(r[ix['zHD']]) for r in R]); cal=np.array([int(r[ix['IS_CALIBRATOR']]) for r in R])
m=(z>0.01)&(cal==0); ZP=z[m]; CPP=Cf[np.ix_(m,m)]; CPPI=inv(CPP); LP=cholesky(CPP)
ONE=np.ones(len(ZP)); A1P=ONE@CPPI@ONE; CppI_ONE=CPPI@ONE
def sn_mu(Om,A): return 5*np.log10((1+ZP)*c_km*np.interp(ZP,ZG,comoving(Om,A)))+25
MU=np.array([[sn_mu(Om,A) for A in A_G] for Om in OM_G]).reshape(-1,len(ZP))  # (P,Nsn)
MUI=(CPPI@MU.T).T; selfI=np.einsum('pi,pi->p',MU,MUI); MU_ONE=MU@CppI_ONE
def sn_chi2_grid(o):
    oI=CPPI@o
    rCr=(o@oI)-2*MU@oI+selfI; B=(ONE@oI)-MU_ONE
    return (rCr-B*B/A1P).reshape(nOm,nA)

# BAO setup
ZB,DMr,DHr,CB=get_arrays(); CBI=inv(CB); LB=cholesky(CB); nb=len(ZB)
def bao_shape(Om,A):
    chi=comoving(Om,A); DM=np.interp(ZB,ZG,chi); DH=1.0/Efun(ZB,Om,A)
    o=np.empty(2*nb); o[0::2]=DM; o[1::2]=DH; return o
SH=np.array([[bao_shape(Om,A) for A in A_G] for Om in OM_G]).reshape(-1,2*nb)
SHI=(CBI@SH.T).T; SH_selfden=np.einsum('pi,pi->p',SH,SHI)
def bao_chi2_grid(o):
    num=SH@(CBI@o); k=num/SH_selfden
    Rr=o[None,:]-k[:,None]*SH
    return np.einsum('pi,ij,pj->p',Rr,CBI,Rr).reshape(nOm,nA)

def truth_sn(A): return sn_mu(0.33,A)
def truth_bao(A): return bao_shape(0.33,A)*(c_km/H0RD)
def draw_sn(o): return o+LP@rng.standard_normal(len(ZP))
def draw_bao(o): return o+LB@rng.standard_normal(2*nb)

def rate(kind,amps=[0.2,0.4,0.6]):
    def fits(osn,obao):
        if kind=='bao':   G=bao_chi2_grid(obao)
        elif kind=='joint': G=sn_chi2_grid(osn)+bao_chi2_grid(obao)
        return G[:,J0].min(),G.min()
    o0s,o0b=truth_sn(0.0),truth_bao(0.0)
    dn=np.clip([np.subtract(*fits(draw_sn(o0s),draw_bao(o0b))) for _ in range(NM)],0,None)
    med=np.median(dn)
    out={}
    for A in amps:
        os,ob=truth_sn(A),truth_bao(A)
        dd=np.array([np.subtract(*fits(draw_sn(os),draw_bao(ob))) for _ in range(NM)])
        out[A]=np.mean(dd>THR3)*100
    print(f"{kind:8} null={med:.3f} | "+"  ".join(f"A={A}:{out[A]:3.0f}%" for A in amps))
    return out

print("#1 SNe+BAO vs BAO alone (real covariances)\n")
rb=rate('bao'); rj=rate('joint')
np.savez('/home/claude/combine.npz',bao=[rb[A] for A in [0.2,0.4,0.6]],
    joint=[rj[A] for A in [0.2,0.4,0.6]],amps=[0.2,0.4,0.6])
print("\nGain from adding SNe to BAO:")
for A in [0.2,0.4,0.6]:
    print(f"  A={A} ({A/(1+A)*100:.0f}% shed): BAO={rb[A]:.0f}%  SNe+BAO={rj[A]:.0f}%  {rj[A]-rb[A]:+.0f}pp")
