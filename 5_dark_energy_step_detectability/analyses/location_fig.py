import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
d=np.load('location.npz')
zt=d['zt']; bao=d['bao']; sne=d['sne']; tracers=d['tracers']
fig,ax=plt.subplots(figsize=(9.5,5.2))
# mark DESI tracer redshifts
for i,zt_tr in enumerate(tracers):
    ax.axvline(zt_tr,ls=':',color='#27ae60',alpha=0.35,lw=1)
ax.text(tracers[0],103,'DESI tracers',fontsize=8,color='#27ae60',alpha=.8)
ax.plot(zt,bao,'^-',color='#27ae60',lw=2.3,ms=8,label='BAO (real DESI DR2)')
ax.plot(zt,sne,'s-',color='#c0392b',lw=2.3,ms=8,label='SNe (real Pantheon+)')
ax.axhline(50,ls='--',color='gray',lw=.7,alpha=.6)
ax.set_xlabel('transition redshift  $z_t$   (when the step happens)')
ax.set_ylabel('3σ detection rate  [%]')
ax.set_title('Where in cosmic history can each probe catch a dark-energy step?\n(fixed amplitude: 37.5% shed)')
ax.legend(fontsize=10,loc='upper right'); ax.grid(alpha=.2); ax.set_ylim(-3,108)
ax.annotate('SNe only competitive\nat very low z',xy=(0.3,22),xytext=(0.55,35),fontsize=8.5,color='#c0392b',
    arrowprops=dict(arrowstyle='->',color='#c0392b',alpha=.6))
ax.annotate('BAO sensitivity peaks\nwhere its tracers are dense',xy=(0.6,87),xytext=(0.85,60),fontsize=8.5,color='#27ae60',
    arrowprops=dict(arrowstyle='->',color='#27ae60',alpha=.6))
plt.tight_layout(); plt.savefig('/mnt/user-data/outputs/location_scan.png',dpi=140,bbox_inches='tight')
print("saved")
print("BAO peak:",zt[np.argmax(bao)],"=",max(bao),"%")
print("BAO falls to <10% beyond zt~1.1; SNe only reach >20% at zt<=0.3")
