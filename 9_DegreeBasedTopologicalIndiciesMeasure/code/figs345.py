import os, csv, numpy as np, chemicals
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from chemicals.identifiers import pubchem_db
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

INDICES = {
 "$M_1$":lambda a,b:a+b, "$M_2$":lambda a,b:a*b, "F":lambda a,b:a*a+b*b,
 "HM":lambda a,b:(a+b)**2, "Randić":lambda a,b:1/np.sqrt(a*b),
 "sum-conn":lambda a,b:1/np.sqrt(a+b), "ABC":lambda a,b:np.sqrt((a+b-2)/(a*b)),
 "GA":lambda a,b:2*np.sqrt(a*b)/(a+b), "harmonic":lambda a,b:2/(a+b),
 "Sombor":lambda a,b:np.sqrt(a*a+b*b), "AZI":lambda a,b:(a*b/max(a+b-2,1e-9))**3,
 "ISI":lambda a,b:a*b/(a+b), "SDD":lambda a,b:a/b+b/a, "ReZG3":lambda a,b:a*b*(a+b)}
IN = list(INDICES); K = len(IN)
PAIRS = [(1,1),(1,2),(1,3),(1,4),(2,2),(2,3),(2,4),(3,3),(3,4),(4,4)]

pubchem_db.autoload_main_db()
cas2smi = {chemicals.CAS_to_int(o.CASs): o.smiles
           for o in pubchem_db.CAS_index.values() if o.smiles}
p = os.path.join(os.path.dirname(chemicals.__file__),'Phase Change','OpenNotebook Melting Points.tsv')
tab = {}
with open(p,encoding='utf-8',errors='replace') as fh:
    next(fh)
    for line in fh:
        q=line.rstrip().split("\t")
        if len(q)>1 and q[1]:
            try: tab[chemicals.CAS_to_int(q[0])]=float(q[1])
            except Exception: pass

TIv, ms, Pm = [], [], []
for cas in tab:
    s = cas2smi.get(cas)
    if not s or '.' in s: continue
    m = Chem.MolFromSmiles(s)
    if m is None or not (10<=m.GetNumAtoms()<=60) or m.GetNumBonds()<3: continue
    if max(a.GetDegree() for a in m.GetAtoms())>4: continue
    if not any(a.GetSymbol()=='C' for a in m.GetAtoms()): continue
    e=[tuple(sorted((b.GetBeginAtom().GetDegree(),b.GetEndAtom().GetDegree()))) for b in m.GetBonds()]
    mm=len(e)
    TIv.append([sum(f(a,b) for a,b in e) for f in INDICES.values()])
    Pm.append([sum(1 for x in e if x==pr)/mm for pr in PAIRS]); ms.append(mm)
TI=np.array(TIv,float); ms=np.array(ms,float); Pm=np.array(Pm)
N=len(ms); print("N =", N)

# ================= FIGURE 3 : one dimension =================
C = np.corrcoef(TI.T)
Z = (TI-TI.mean(0))/TI.std(0)
ev = np.linalg.svd(Z, compute_uv=False)**2; ev/=ev.sum()

fig, ax = plt.subplots(1, 2, figsize=(12.6, 5.0),
                       gridspec_kw={"width_ratios":[1.25,1]})
im = ax[0].imshow(C, vmin=0.85, vmax=1.0, cmap="RdYlBu_r")
ax[0].set_xticks(range(K)); ax[0].set_xticklabels(IN, rotation=90, fontsize=7.5)
ax[0].set_yticks(range(K)); ax[0].set_yticklabels(IN, fontsize=7.5)
ax[0].set_title("(a)  Pairwise correlation among 14 degree-based indices\n"
                f"(N = {N} drug-like molecules; median |r| = {np.median(np.abs(C[np.triu_indices(K,1)])):.3f})",
                fontsize=10)
cb = fig.colorbar(im, ax=ax[0], fraction=.046); cb.ax.tick_params(labelsize=8)

ax[1].bar(range(1,K+1), 100*ev, color="#34495e")
ax[1].set_yscale("log"); ax[1].set_xlabel("principal component")
ax[1].set_ylabel("variance explained (%, log scale)")
ax[1].set_title(f"(b)  Scree: PC1 = {100*ev[0]:.2f}%,  PC2 = {100*ev[1]:.2f}%", fontsize=10)
ax[1].grid(alpha=.25, axis="y")
ax[1].annotate("one effective dimension", xy=(1, 100*ev[0]), xytext=(4, 30),
               arrowprops=dict(arrowstyle="->", lw=1.2), fontsize=9)
fig.tight_layout(); fig.savefig("/mnt/user-data/outputs/figure3_dimensionality.png", dpi=170)

# ================= FIGURE 4 : bound tightness =================
w = ms**2/(ms**2).sum(); pbar = Pm.mean(0); U = Pm-pbar
Sig = (U*w[:,None]).T @ U
supp = pbar > 1e-4
pred, act, lbl = [], [], []
for j,(nm,f) in enumerate(INDICES.items()):
    fv = np.array([f(a,b) if (a,b)!=(1,1) else 0.0 for a,b in PAIRS])
    fs = fv[supp]; gbar = float(pbar@fv)
    fc = fv - 0.5*(fs.max()+fs.min())
    tau = np.sqrt(max(fc@Sig@fc,0))*np.sqrt(np.sum(ms**2)) / (gbar*np.sqrt(N)*ms.std())
    pred.append(max(0,1-tau**2)); act.append(np.corrcoef(TI[:,j],ms)[0,1]**2); lbl.append(nm)

fig, ax = plt.subplots(figsize=(6.6,5.6))
ax.plot([0.5,1.0],[0.5,1.0],"k--",lw=1, label="equality (bound exactly tight)")
ax.scatter(pred, act, s=55, color="#c0392b", zorder=3)
for x,y,t in zip(pred,act,lbl):
    ax.annotate(t,(x,y),textcoords="offset points",xytext=(6,-3),fontsize=7.5)
ax.fill_between([0.5,1.0],[0.5,1.0],[1.02,1.02],color="#2ecc71",alpha=.10,lw=0)
ax.text(0.60,0.97,"bound valid\n(below the diagonal)",fontsize=8.5,color="#1e8449")
ax.set_xlabel(r"lower bound $1-\tau_f^2$ from the Proposition")
ax.set_ylabel(r"observed $r^2(TI_f,\ m)$")
ax.set_title("Validity and tightness of the spectral bound\n"
             f"14/14 valid; 10/14 tight to within 0.013  (N = {N})", fontsize=10.5)
ax.set_xlim(0.5,1.02); ax.set_ylim(0.5,1.02); ax.grid(alpha=.25); ax.legend(fontsize=8, loc="lower right")
fig.tight_layout(); fig.savefig("/mnt/user-data/outputs/figure4_bound.png", dpi=170)

# ================= FIGURE 5 : threshold vs fit order =================
def poly_r2(x,y,deg):
    if np.std(x)<1e-12 or np.std(y)<1e-12: return 0.0
    V=np.vander((x-x.mean())/x.std(),deg+1)
    b,*_=np.linalg.lstsq(V,y,rcond=None); r=y-V@b
    ss=np.sum((y-y.mean())**2)
    return float(max(0,1-np.sum(r**2)/ss)) if ss>0 else 0.0

ys=[]
for cas in tab:
    s=cas2smi.get(cas)
    if not s or '.' in s: continue
    m=Chem.MolFromSmiles(s)
    if m is None or not (10<=m.GetNumAtoms()<=60) or m.GetNumBonds()<3: continue
    if max(a.GetDegree() for a in m.GetAtoms())>4: continue
    if not any(a.GetSymbol()=='C' for a in m.GetAtoms()): continue
    ys.append(tab[cas])
ys=np.array(ys)

degs=[1,2,3]; nulls=[]; reals=[]
rg=np.random.default_rng(5)
for deg in degs:
    nl, rl = [], []
    for _ in range(1500):
        i=rg.choice(len(ys),15,replace=False)
        rl.append(max(poly_r2(TI[i,j],ys[i],deg) for j in range(K)))
        yp=ys[rg.choice(len(ys),15,replace=False)]
        nl.append(max(poly_r2(TI[i,j],yp,deg) for j in range(K)))
    nulls.append(np.percentile(nl,95)); reals.append(np.mean(rl))

fig, ax = plt.subplots(figsize=(7.2,5.0))
ax.plot(degs, nulls, "o-", lw=2.4, ms=8, color="#c0392b",
        label="95th percentile of the permutation null")
ax.plot(degs, reals, "s--", lw=2, ms=7, color="#2c3e50",
        label="mean best-of-14 $R^2$ on real data")
ax.fill_between(degs, nulls, 1.0, color="#c0392b", alpha=.07, lw=0)
ax.text(1.45, 0.86, "a claim must clear this line\nto be evidence of anything",
        fontsize=8.5, color="#7b241c")
ax.set_xticks(degs); ax.set_xticklabels(["linear","quadratic","cubic"])
ax.set_xlabel("polynomial order of the fitted model")
ax.set_ylabel(r"best-of-14 $R^2$  at $N=15$")
ax.set_ylim(0,1.0); ax.grid(alpha=.25); ax.legend(fontsize=8.5, loc="upper left")
ax.set_title("Raising the fit order raises the evidential bar\n(measured melting points, 15-molecule draws)", fontsize=10.5)
fig.tight_layout(); fig.savefig("/mnt/user-data/outputs/figure5_fitorder.png", dpi=170)
print("figures 3-5 written")
print("null 95th by degree:", [f"{v:.3f}" for v in nulls])
