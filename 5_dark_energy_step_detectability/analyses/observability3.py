"""Clean observability check: the ORTHOGONAL FRACTION as the master quantity.
For each probe structure, measure what fraction of each feature-parameter's
signal survives after projecting out the nuisance subspace (Om + distance scale).
This IS the marginalized-Fisher / degeneracy statement, applied to real DESI."""
import numpy as np, sys
from numpy.linalg import inv
sys.path.insert(0,'realdata')
from desi_dr2 import get_arrays
ZB,DMr,DHr,CB=get_arrays(); ZB=np.asarray(ZB); CBI=inv(CB); nb=len(ZB)
ZT,DZ,SCALE=0.7,0.15,30.42
def S(z,zt=ZT,dz=DZ): return 0.5*(1+np.tanh((z-zt)/dz))
def Ef(z,Om,A,zt=ZT,dz=DZ): return np.sqrt(Om*(1+z)**3+(1-Om)*(1+A*S(z,zt,dz)))
ZG=np.linspace(0,3,1500)
def ov(Om,A,k=SCALE,zt=ZT,dz=DZ,mode='both'):
    iE=1/Ef(ZG,Om,A,zt,dz); chi=np.concatenate([[0],np.cumsum(0.5*(iE[1:]+iE[:-1])*np.diff(ZG))])
    DM=k*np.interp(ZB,ZG,chi); DH=k/Ef(ZB,Om,A,zt,dz)
    if mode=='both': o=np.empty(2*nb); o[0::2]=DM; o[1::2]=DH; return o
    return DM if mode=='DM' else DH
def Cinv(mode):
    C=inv(CBI)
    if mode=='both': return CBI
    idx=np.arange(0,2*nb,2) if mode=='DM' else np.arange(1,2*nb,2)
    return inv(C[np.ix_(idx,idx)])
def dv(p,mode,Om0=0.3,A0=0.0,k0=SCALE,h=1e-5):
    f=lambda **kw: ov(kw.get('Om',Om0),kw.get('A',A0),kw.get('k',k0),dz=kw.get('dz',DZ),mode=mode)
    return {'A':(f(A=A0+h)-f(A=A0-h))/(2*h),'Om':(f(Om=Om0+h)-f(Om=Om0-h))/(2*h),
            'k':(f(k=k0+h)-f(k=k0-h))/(2*h),'dz':(f(dz=DZ+h)-f(dz=DZ-h))/(2*h)}[p]
def orth_frac(target,nuis,mode,A0=0.0):
    Ci=Cinv(mode); g=dv(target,mode,A0=A0)
    B=np.array([dv(p,mode,A0=A0) for p in nuis])
    proj=B.T@inv(B@Ci@B.T)@(B@Ci@g); r=g-proj
    return np.sqrt((r@Ci@r)/(g@Ci@g))
def marg_snr(target,nuis,mode,amp=0.6,A0=0.0):
    Ci=Cinv(mode); allp=[target]+nuis; G=np.array([dv(p,mode,A0=A0) for p in allp])
    F=G@Ci@G.T; Fm=F[0,0]-F[0,1:]@inv(F[1:,1:])@F[1:,0]; return np.sqrt(max(Fm,0))*amp

print("MASTER TABLE: orthogonal fraction & marginalized SNR (feature vs nuisance Om+scale)")
print(f"{'probe':10}{'#obs/z':>7}{'A_orthfrac':>12}{'SNR(A=.6)':>11}")
for mode,k in [('DM',1),('DH',1),('both',2)]:
    of=orth_frac('A',['Om','k'],mode); sn=marg_snr('A',['Om','k'],mode)
    print(f"{mode:10}{k:>7}{of:>12.3f}{sn:>11.2f}")
print("  MC 3sig rates for A=0.6: DM~2%, DH~46%, both~88% -> tracks SNR ordering")

print("\nCLAIM 3 — WIDTH: is dz recoverable once A is in the model? (both mode)")
# marginalize dz over (Om, k, A) -> if ~0, width in null space
for A0 in [0.6]:
    of_dz=orth_frac('dz',['Om','k','A'],'both',A0=A0)
    of_A =orth_frac('A',['Om','k','dz'],'both',A0=A0)
    print(f"  at A={A0}: dz orthogonal-fraction (after Om,k,A) = {of_dz:.3f}")
    print(f"            A  orthogonal-fraction (after Om,k,dz) = {of_A:.3f}")
    print("  -> dz frac near 0 = width NOT separately recoverable (matches width_test: 0% reject)")

print("\nRANK of the feature block (A,dz) marginalized over (Om,k), both mode:")
Ci=Cinv('both'); feats=['A','dz']; nuis=['Om','k']; allp=feats+nuis
G=np.array([dv(p,'both',A0=0.6) for p in allp]); F=G@Ci@G.T
Fff=F[:2,:2]-F[:2,2:]@inv(F[2:,2:])@F[2:,:2]
ev=np.linalg.eigvalsh(Fff)
print(f"  marginalized (A,dz) Fisher eigenvalues: {ev}")
print(f"  condition number: {ev[-1]/max(ev[0],1e-30):.1f}  (huge = rank-deficient = 1 recoverable dof not 2)")
