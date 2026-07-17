"""
#2 MECHANISM PROOF: is BAO's step-sensitivity really from the DM+DH pair?
Test by running BAO detection with:
  (a) DM/rd ONLY (single distance, SNe-like - should COLLAPSE to blindness)
  (b) DH/rd ONLY (single distance, but the one that shows the step directly)
  (c) BOTH (the real result, for reference)
Same validated grid method, real DESI DR2 covariance sub-blocks, null-checked.
"""
import numpy as np, sys
from numpy.linalg import inv, cholesky
from scipy.stats import chi2 as chi2d, ncx2
sys.path.insert(0,'realdata')
from desi_dr2 import get_arrays
c_km=299792.458; ZT=0.7; DZ=0.15; THR3=chi2d.ppf(0.9973,1)
def Sfun(z): return 0.5*(1+np.tanh((z-ZT)/DZ))
def Efun(z,Om,A): return np.sqrt(Om*(1+z)**3+(1-Om)*(1.0+A*Sfun(z)))
ZG=np.linspace(0,3.0,900)
def comoving(Om,A):
    invE=1.0/Efun(ZG,Om,A); return np.concatenate([[0],np.cumsum(0.5*(invE[1:]+invE[:-1])*np.diff(ZG))])
OM_G=np.linspace(0.15,0.45,41); A_G=np.round(np.arange(-0.9,3.0+1e-9,0.05),4); J0=int(np.argmin(np.abs(A_G)))
assert abs(A_G[J0])<1e-9
rng=np.random.default_rng(42); NM=400; H0RD=67.0*147.09
ZB,DMr,DHr,CB=get_arrays()   # CB is 12x12: interleaved [DM,DH] per tracer
n=len(ZB)

# indices: DM at 0,2,4..; DH at 1,3,5..
iDM=np.arange(0,2*n,2); iDH=np.arange(1,2*n,2)

def make_probe(idx):
    """Build a probe using only the observable indices in idx (subset of 0..2n-1)."""
    C=CB[np.ix_(idx,idx)]; CI=inv(C); L=cholesky(C)
    def full_shape(Om,A):
        chi=comoving(Om,A); DM=np.interp(ZB,ZG,chi); DH=1.0/Efun(ZB,Om,A)
        o=np.empty(2*n); o[0::2]=DM; o[1::2]=DH; return o[idx]
    SH=np.array([[full_shape(Om,A) for A in A_G] for Om in OM_G]).reshape(-1,len(idx))
    def fits(o):
        num=SH@(CI@o); den=np.einsum('pi,pi->p',SH,(CI@SH.T).T); k=num/den
        R=o[None,:]-k[:,None]*SH; G=np.einsum('pi,ij,pj->p',R,CI,R).reshape(len(OM_G),len(A_G))
        return G[:,J0].min(),G.min()
    def truth(A): return full_shape(0.33,A)*(c_km/H0RD)
    def draw(o): return o+L@rng.standard_normal(len(idx))
    return fits,truth,draw,CI

def run(name,idx,amps=[0.4,0.6,0.8]):
    fits,truth,draw,CI=make_probe(idx)
    # null check
    o0=truth(0.0); dn=np.clip([np.subtract(*fits(draw(o0))) for _ in range(NM)],0,None)
    med=np.median(dn)
    assert med<1.2, f"{name} null bad {med}"
    out={}
    for A in amps:
        oA=truth(A); dd=np.array([np.subtract(*fits(draw(oA))) for _ in range(NM)])
        out[A]=np.mean(dd>THR3)*100
    print(f"{name:20} null={med:.3f} | " + "  ".join(f"A={A}:{out[A]:3.0f}%" for A in amps))
    return out

print("MECHANISM TEST — does BAO need BOTH DM and DH?\n")
r_dm  = run("DM only (SNe-like)", iDM)
r_dh  = run("DH only", iDH)
r_both= run("DM+DH (full BAO)", np.arange(2*n))
np.savez('/home/claude/mechanism.npz',
    dm=[r_dm[A] for A in [0.4,0.6,0.8]], dh=[r_dh[A] for A in [0.4,0.6,0.8]],
    both=[r_both[A] for A in [0.4,0.6,0.8]], amps=[0.4,0.6,0.8])
print("\nInterpretation:")
print(f"  If DM-only collapses toward SNe blindness -> confirms 2-observable structure is the cause")
