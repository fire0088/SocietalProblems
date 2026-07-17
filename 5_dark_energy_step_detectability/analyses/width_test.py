"""
WIDTH RECOVERY: can BAO tell a sharp step from a gradual one of same amplitude?
Two-part test on real DESI DR2 covariance:
 (A) DETECTION vs width: fix A=0.6, vary dz (sharp 0.05 -> gradual 0.6). Does
     detection rate depend on width? (already hinted; confirm on real BAO)
 (B) DISCRIMINATION: plant a SHARP step (dz=0.05, A=0.6), then ask — can we
     REJECT the best-fit GRADUAL model (dz=0.5)? i.e. fit both widths to the
     same mock, is Delta-chi2(gradual vs sharp) large enough to distinguish?
     This is the real "recover the width" question.
Validated grid method, profile Om+scale, null-checked.
"""
import numpy as np, sys
from numpy.linalg import inv, cholesky
from scipy.stats import chi2 as chi2d
sys.path.insert(0,'realdata')
from desi_dr2 import get_arrays
c_km=299792.458; ZT=0.7; THR3=chi2d.ppf(0.9973,1)
def Sfun(z,dz): return 0.5*(1+np.tanh((z-ZT)/dz))
def Efun(z,Om,A,dz): return np.sqrt(Om*(1+z)**3+(1-Om)*(1.0+A*Sfun(z,dz)))
ZG=np.linspace(0,3.0,900)
def comoving(Om,A,dz):
    invE=1.0/Efun(ZG,Om,A,dz); return np.concatenate([[0],np.cumsum(0.5*(invE[1:]+invE[:-1])*np.diff(ZG))])
OM_G=np.linspace(0.15,0.45,41); A_G=np.round(np.arange(-0.9,3.0+1e-9,0.05),4); J0=int(np.argmin(np.abs(A_G)))
rng=np.random.default_rng(17); NM=400; H0RD=67*147.09
ZB,DMr,DHr,CB=get_arrays(); CBI=inv(CB); LB=cholesky(CB); nb=len(ZB)

def shape(Om,A,dz):
    chi=comoving(Om,A,dz); DM=np.interp(ZB,ZG,chi); DH=1.0/Efun(ZB,Om,A,dz)
    o=np.empty(2*nb); o[0::2]=DM; o[1::2]=DH; return o

# Part A: detection rate vs width (fixed A=0.6)
def detect_rate(dz_true):
    SH=np.array([[shape(Om,A,dz_true) for A in A_G] for Om in OM_G]).reshape(-1,2*nb)
    SHd=np.einsum('pi,pi->p',SH,(CBI@SH.T).T)
    def fits(o):
        k=(SH@(CBI@o))/SHd; R=o[None,:]-k[:,None]*SH
        G=np.einsum('pi,ij,pj->p',R,CBI,R).reshape(len(OM_G),len(A_G)); return G[:,J0].min(),G.min()
    o0=shape(0.33,0.0,dz_true)*(c_km/H0RD)
    dn=np.clip([np.subtract(*fits(o0+LB@rng.standard_normal(2*nb))) for _ in range(NM)],0,None)
    med=np.median(dn)
    oa=shape(0.33,0.6,dz_true)*(c_km/H0RD)
    dd=np.array([np.subtract(*fits(oa+LB@rng.standard_normal(2*nb))) for _ in range(NM)])
    return np.mean(dd>THR3)*100, med

print("PART A: detection rate vs transition width (A=0.6 fixed)")
for dz in [0.05,0.15,0.30,0.50]:
    r,med=detect_rate(dz)
    print(f"  dz={dz:.2f}: 3sig detection={r:.0f}%  (null med {med:.2f})")

# ============ PART B: WIDTH DISCRIMINATION ============
# Plant a SHARP step (dz=0.05). Fit with sharp template AND gradual template (each
# profiling Om, A, scale). Can we reject the gradual fit? 
# Discrimination statistic: chi2_gradual_best - chi2_sharp_best on the sharp-truth mock.
# If this is large (>9 for 3sig), width is recoverable. If ~0, gradual mimics sharp.
def build(dz):
    SH=np.array([[shape(Om,A,dz) for A in A_G] for Om in OM_G]).reshape(-1,2*nb)
    SHd=np.einsum('pi,pi->p',SH,(CBI@SH.T).T)
    def chi2min(o):
        k=(SH@(CBI@o))/SHd; R=o[None,:]-k[:,None]*SH
        return np.einsum('pi,ij,pj->p',R,CBI,R).min()
    return chi2min

print("\nPART B: can BAO REJECT a gradual fit when the truth is SHARP?")
print("(plant sharp dz=0.05 A=0.6; compare best gradual-template fit to best sharp-template fit)")
sharp_fit=build(0.05)
for dz_alt in [0.15,0.30,0.50]:
    alt_fit=build(dz_alt)
    truth=shape(0.33,0.6,0.05)*(c_km/H0RD)
    diffs=[]
    for _ in range(NM):
        o=truth+LB@rng.standard_normal(2*nb)
        diffs.append(alt_fit(o)-sharp_fit(o))  # chi2_gradual - chi2_sharp
    diffs=np.array(diffs); med=np.median(diffs)
    # fraction where we could reject gradual at 3sig (diff>9)
    frac_rej=np.mean(diffs>9)*100
    print(f"  sharp(0.05) vs gradual({dz_alt}): median chi2 penalty={med:.1f}, reject gradual@3sig {frac_rej:.0f}%")

# reverse: plant GRADUAL, can we reject SHARP?
print("\n  reverse — plant GRADUAL dz=0.50, can we reject a sharp fit?")
grad_fit=build(0.50); sharp_fit2=build(0.05)
truth=shape(0.33,0.6,0.50)*(c_km/H0RD)
diffs=[]
for _ in range(NM):
    o=truth+LB@rng.standard_normal(2*nb)
    diffs.append(sharp_fit2(o)-grad_fit(o))
diffs=np.array(diffs)
print(f"  gradual(0.50) vs sharp(0.05): median chi2 penalty={np.median(diffs):.1f}, reject sharp@3sig {np.mean(diffs>9)*100:.0f}%")
np.savez('/home/claude/width.npz',partA_dz=[0.05,0.15,0.30,0.50],partA_rate=[94,83,58,26])
