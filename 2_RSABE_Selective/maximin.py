import numpy as np
from scipy import stats
from scipy.optimize import linprog
ln125=np.log(1.25); Theta=(ln125/0.25)**2; alpha=0.05
n=8; nu=n-2; tau=1.0/n     # the worst case for F-test (where LP beat it)
Dg=np.linspace(-1.3,1.3,200); sg=np.linspace(0.04,1.7,170)
area=(Dg[1]-Dg[0])*(sg[1]-sg[0]); DD,SS=np.meshgrid(Dg,sg,indexing='ij')
def dens(D,s): return (stats.norm.pdf(DD,D,np.sqrt(tau)*s)*stats.chi2.pdf(nu*SS**2/s**2,nu)*(nu/s**2)*2*SS).ravel()
W=(DD**2/(tau*SS**2)).ravel()
wq=stats.ncf.ppf(alpha,1,nu,nc=Theta*n); phiF=(W<wq).astype(float)

# rebuild the LP test that BEAT the F-test at alternative (rho1=0.5*Theta, sigma=0.30)
sig_c=np.linspace(0.10,0.62,24)
A_eq=np.array([dens(np.sqrt(Theta)*s,s)*area for s in sig_c]); b_eq=A_eq@phiF
rho1=Theta*0.5; sig_star=0.30
cobj=-dens(np.sqrt(rho1)*sig_star,sig_star)*area
res=linprog(cobj,A_eq=A_eq,b_eq=b_eq,bounds=[(0,1)]*len(phiF),method='highs')
phiLP=res.x

# POWER ACROSS THE WHOLE ALTERNATIVE SLICE {rho = rho1}, i.e. Delta = sqrt(rho1)*sigma, varying sigma
print(f"n={n} (nu={nu}). Alternative slice rho={rho1:.3f} (=0.5*Theta). Delta=sqrt(rho)*sigma.")
print("F-test power should be CONSTANT in sigma (invariance). LP test was tuned at sigma=0.30.\n")
print(f"{'sigma':>7} {'F-test power':>13} {'LP-test power':>14} {'LP - F':>9}")
sigmas=[0.15,0.20,0.25,0.30,0.35,0.45,0.60,0.80]
pf_list=[];pl_list=[]
for s in sigmas:
    f=dens(np.sqrt(rho1)*s,s)*area
    pf=f@phiF; pl=f@phiLP; pf_list.append(pf); pl_list.append(pl)
    mark="  <-- LP tuned here" if abs(s-sig_star)<1e-9 else ""
    print(f"{s:>7.2f} {pf:>13.4f} {pl:>14.4f} {pl-pf:>+9.4f}{mark}")
print(f"\nWORST-CASE (minimum) power over the slice:")
print(f"   F-test : {min(pf_list):.4f}")
print(f"   LP-test: {min(pl_list):.4f}")
print(f"   => F-test maximin advantage: {min(pf_list)-min(pl_list):+.4f}")
print("\n(If positive, the LP's local win is bought by sacrificing power elsewhere on the slice,")
print(" and the F-test is the maximin choice -- exactly as Hunt-Stein predicts.)")
