"""Final figure — reads ONLY the freshly-computed r_*.npz. No hardcoded rates."""
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import sys; sys.path.insert(0,'realdata'); from desi_dr2 import get_arrays

AMPS=[0.2,0.4,0.6,0.8,1.0]; shed=np.array([A/(1+A)*100 for A in AMPS])
def load(f): d=np.load(f); return np.array([float(d[str(A)]) for A in AMPS])
cc=load('r_cc.npz'); sne=load('r_sne.npz'); bao=load('r_bao.npz')
print("Loaded fresh rates:")
print(f"  shed%: {np.round(shed,1)}")
print(f"  CC:  {cc}\n  SNe: {sne}\n  BAO: {bao}")

CC='#2980b9'; SNE='#c0392b'; BAO='#27ae60'
fig,ax=plt.subplots(1,2,figsize=(13.5,5.4))

ax[0].plot(shed,cc,'o-',color=CC,lw=2.2,ms=7,label='Cosmic chronometers (Moresco)')
ax[0].plot(shed,sne,'s-',color=SNE,lw=2.2,ms=7,label='SNe (real Pantheon+, full cov)')
ax[0].plot(shed,bao,'^-',color=BAO,lw=2.2,ms=7,label='BAO (real DESI DR2)')
ax[0].axhline(99.73,ls=':',c='gray',lw=1); ax[0].axhline(50,ls='--',c='gray',lw=.7,alpha=.6)
ax[0].set_xlabel('dark energy density shed at transition  [%]')
ax[0].set_ylabel('3σ detection rate  [%]')
ax[0].set_title('(a)  Detecting a dark-energy step — real data\n(fixed location z=0.7, width Δz=0.15)',fontsize=11)
ax[0].legend(fontsize=9,loc='center left'); ax[0].grid(alpha=.2); ax[0].set_ylim(-4,105)
ax[0].annotate('BAO sees it',xy=(shed[3],bao[3]),xytext=(30,72),fontsize=9,color=BAO,
    arrowprops=dict(arrowstyle='->',color=BAO,alpha=.7))
ax[0].annotate('SNe & CC blind',xy=(shed[3],sne[3]),xytext=(38,22),fontsize=9,color=SNE,
    arrowprops=dict(arrowstyle='->',color=SNE,alpha=.7))

# Panel B: mechanism — step imprint on DM vs DH
ZB,_,_,_=get_arrays(); c_km=299792.458
def Sfun(z): return 0.5*(1+np.tanh((z-0.7)/0.15))
def Efun(z,A): return np.sqrt(0.33*(1+z)**3+0.67*(1.0+A*Sfun(z)))
zg=np.linspace(0.02,2.5,300)
def DMDH(A):
    invE=1.0/Efun(zg,A); chi=np.concatenate([[0],np.cumsum(0.5*(invE[1:]+invE[:-1])*np.diff(zg))])
    return chi,1.0/Efun(zg,A)
chi0,DH0=DMDH(0.0); chiS,DHS=DMDH(0.8)
with np.errstate(invalid='ignore',divide='ignore'):
    dM=100*(chiS-chi0)/chi0
ax[1].plot(zg,dM,color='#8e44ad',lw=2.2,label='D$_M$/r$_d$  (transverse, integrated)')
ax[1].plot(zg,100*(DHS-DH0)/DH0,color='#e67e22',lw=2.2,label='D$_H$/r$_d$  (radial, local H)')
ax[1].axhline(0,color='gray',lw=.6); ax[1].axvline(0.7,ls=':',c='k',alpha=.4)
ax[1].scatter(ZB,[100*(np.interp(z,zg,DHS)-np.interp(z,zg,DH0))/np.interp(z,zg,DH0) for z in ZB],
    color='#e67e22',s=45,zorder=5,edgecolor='white',label='DESI DR2 tracers')
ax[1].set_xlabel('redshift z'); ax[1].set_ylabel('% change from ΛCDM (44% step)')
ax[1].set_title('(b)  Why BAO wins: step shows in D$_H$, hides in D$_M$',fontsize=11)
ax[1].legend(fontsize=8.5); ax[1].grid(alpha=.2)
plt.tight_layout(); plt.savefig('/mnt/user-data/outputs/real_detectability.png',dpi=140,bbox_inches='tight')
print("\nsaved figure from fresh data only")
