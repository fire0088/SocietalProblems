import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
d=np.load('exclusion.npz')
A=d['A_g']; dchi=d['dchi']; lo95,hi95=float(d['lo95']),float(d['hi95']); lo99,hi99=float(d['lo99']),float(d['hi99'])
from scipy.stats import chi2 as chi2d
t68,t95,t99=chi2d.ppf(0.68,1),chi2d.ppf(0.95,1),chi2d.ppf(0.997,1)
fig,ax=plt.subplots(figsize=(9,5.2))
m=(A>-1.0)&(A<1.5)
ax.plot(A[m],dchi[m],color='#27ae60',lw=2.5)
for t,lab,c in [(t68,'68%','#999'),(t95,'95%','#e67e22'),(t99,'99.7%','#c0392b')]:
    ax.axhline(t,ls='--',lw=1,color=c,label=f'{lab} CL')
# shade allowed region at 95%
ax.axvspan(lo95,hi95,color='#27ae60',alpha=0.08)
ax.axvline(0,ls=':',color='k',alpha=.5)
ax.text(0.02,ax.get_ylim()[1]*0.92,'ΛCDM\n(no step)',fontsize=8,color='k',alpha=.6)
# annotate excluded regions
ax.annotate('excluded:\nthinning >19% shed',xy=(hi95,t95),xytext=(0.55,7),fontsize=8.5,color='#c0392b',
    arrowprops=dict(arrowstyle='->',color='#c0392b',alpha=.6))
ax.annotate('excluded:\nthickening',xy=(lo95,t95),xytext=(-0.85,7),fontsize=8.5,color='#c0392b',
    arrowprops=dict(arrowstyle='->',color='#c0392b',alpha=.6))
ax.set_xlabel('step amplitude  A   (negative = thickening, positive = thinning)')
ax.set_ylabel(r'$\Delta\chi^2$  (profile likelihood)')
ax.set_title('What DESI DR2 already excludes: a sudden dark-energy step at z=0.7')
ax.legend(fontsize=9,loc='upper center'); ax.grid(alpha=.2); ax.set_ylim(0,12)
plt.tight_layout(); plt.savefig('/mnt/user-data/outputs/exclusion.png',dpi=140,bbox_inches='tight')
print(f"95% allowed: A in [{lo95:.2f},{hi95:.2f}] = shed/gain within [{lo95/(1+lo95)*100:.0f}%,{hi95/(1+hi95)*100:.0f}%]")
print("saved")
