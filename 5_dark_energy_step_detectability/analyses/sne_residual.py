"""
Address the Adversarial Reviewer directly: inject a 44% step (A=0.8) into the SNe
Hubble diagram, refit LCDM (floating Om + offset), and show the residuals across
z=0.7. Demonstrate that the Om shift absorbs the step WITHOUT breaking the low-z fit
= the mechanism of SNe blindness, made visual.
"""
import numpy as np, sys, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from numpy.linalg import inv
from scipy.integrate import quad
sys.path.insert(0,'realdata')
dat=np.genfromtxt('realdata/PantheonSH0ES.dat',names=True,dtype=None,encoding=None)
# Hubble-flow, non-calibrator
m=(dat['IS_CALIBRATOR']==0)&(dat['zHD']>0.01)
z=dat['zHD'][m]; mu_obs=dat['MU_SH0ES'][m]; mu_err=dat['MU_SH0ES_ERR_DIAG'][m]
order=np.argsort(z); z=z[order]; mu_obs=mu_obs[order]; mu_err=mu_err[order]
print(f'{len(z)} Hubble-flow SNe, z={z.min():.3f}-{z.max():.2f}')
c=299792.458; ZT,DZ=0.7,0.15
def S(zz): return 0.5*(1+np.tanh((zz-ZT)/DZ))
def Ez(zz,Om,A): return np.sqrt(Om*(1+zz)**3+(1-Om)*(1+A*S(zz)))
def mu_model(zz,Om,A,M):
    # luminosity distance modulus (H0 absorbed in M offset)
    dl=np.array([(1+zi)*quad(lambda x:1/Ez(x,Om,A),0,zi)[0] for zi in zz])
    return 5*np.log10(dl)+M
# 1. inject a 44% step (A=0.8) into a fiducial truth (Om=0.3), make mock 'observed' mu
Om_true=0.30; A_true=0.8
mu_true=mu_model(z,Om_true,A_true,0.0)
M_off=np.median(mu_obs-mu_true)  # align to real data zeropoint
mu_injected=mu_true+M_off        # noiseless injected-step data (show mechanism cleanly)
# 2. fit LCDM (A=0) floating Om + M to the injected-step data
from scipy.optimize import minimize
def chi2(p):
    Om,M=p; r=mu_injected-mu_model(z,Om,0.0,M); return np.sum((r/mu_err)**2)
res=minimize(chi2,[0.3,M_off],method='Nelder-Mead')
Om_fit,M_fit=res.x
print(f'Injected A=0.8 step at Om={Om_true}. Best LCDM (no-step) fit: Om={Om_fit:.3f}')
print(f'  Om shifted by {Om_fit-Om_true:+.3f} to absorb the step')
# 3. residuals of the no-step fit to the stepped data
resid=mu_injected-mu_model(z,Om_fit,0.0,M_fit)
chi2_val=np.sum((resid/mu_err)**2); dof=len(z)-2
print(f'  chi2/dof of no-step fit to stepped data = {chi2_val:.1f}/{dof} = {chi2_val/dof:.3f}')
print(f'  RMS residual = {np.std(resid)*1000:.1f} mmag; typical SNe error {np.median(mu_err)*1000:.0f} mmag')
# binned residuals for clarity
bins=np.linspace(0.01,z.max(),16); bc=0.5*(bins[1:]+bins[:-1])
idx=np.digitize(z,bins); rb=[np.mean(resid[idx==i]) for i in range(1,len(bins))]
rbe=[np.std(resid[idx==i])/np.sqrt(max((idx==i).sum(),1)) for i in range(1,len(bins))]

fig,ax=plt.subplots(2,1,figsize=(9,7),height_ratios=[2,1],sharex=True)
# top: the step in mu vs the absorbing fit
zsm=np.linspace(0.01,z.max(),200)
ax[0].plot(zsm,mu_model(zsm,Om_true,A_true,M_off)-mu_model(zsm,0.3,0,M_off),color='#c0392b',lw=2,label='injected 44% step (A=0.8)')
ax[0].plot(zsm,mu_model(zsm,Om_fit,0,M_fit)-mu_model(zsm,0.3,0,M_off),color='#2980b9',lw=2,ls='--',label=f'best no-step fit (Ωm={Om_fit:.3f})')
ax[0].axvline(0.7,color='gray',ls=':',alpha=0.6); ax[0].set_ylabel('Δμ vs fiducial [mag]')
ax[0].legend(fontsize=9); ax[0].set_title('How a shifted Ωm absorbs a 44% dark-energy step in SNe')
# bottom: residuals of the no-step fit
ax[1].errorbar(bc,np.array(rb)*1000,yerr=np.array(rbe)*1000,fmt='o',color='#27ae60',ms=5)
ax[1].axhline(0,color='k',lw=0.8); ax[1].axvline(0.7,color='gray',ls=':',alpha=0.6)
ax[1].set_ylabel('binned residual [mmag]'); ax[1].set_xlabel('redshift z')
ax[1].set_title(f'Residuals of no-step fit: RMS {np.std(resid)*1000:.0f} mmag, χ²/dof={chi2_val/dof:.2f} (no visible kink at z=0.7)')
plt.tight_layout(); plt.savefig('/mnt/user-data/outputs/sne_residual.png',dpi=140,bbox_inches='tight')
print('saved figure')
np.savez('/home/claude/sne_residual.npz',Om_fit=Om_fit,chi2=chi2_val,dof=dof,rms=np.std(resid))
