import numpy as np, csv, os, chemicals
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from chemicals.identifiers import pubchem_db
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen
RDLogger.DisableLog('rdApp.*')

INDICES = {
 "M1":lambda a,b:a+b, "M2":lambda a,b:a*b, "F":lambda a,b:a*a+b*b,
 "HM":lambda a,b:(a+b)**2, "Randic":lambda a,b:1/np.sqrt(a*b),
 "sum-conn":lambda a,b:1/np.sqrt(a+b), "ABC":lambda a,b:np.sqrt((a+b-2)/(a*b)),
 "GA":lambda a,b:2*np.sqrt(a*b)/(a+b), "harmonic":lambda a,b:2/(a+b),
 "Sombor":lambda a,b:np.sqrt(a*a+b*b), "AZI":lambda a,b:(a*b/max(a+b-2,1e-9))**3,
 "ISI":lambda a,b:a*b/(a+b), "SDD":lambda a,b:a/b+b/a, "ReZG3":lambda a,b:a*b*(a+b)}
INAMES = list(INDICES); K = len(INAMES)

def feats(mol):
    e = [(b.GetBeginAtom().GetDegree(), b.GetEndAtom().GetDegree()) for b in mol.GetBonds()]
    return [mol.GetNumAtoms()] + [sum(f(a,b) for a,b in e) for f in INDICES.values()]
def r2(x,y):
    if np.std(x)<1e-12 or np.std(y)<1e-12: return 0.0
    return float(np.corrcoef(x,y)[0,1]**2)

# ---- ESOL for panel (a) ----
X, W = [], []
for row in csv.DictReader(open("data/delaney-processed.csv")):
    m = Chem.MolFromSmiles(row["smiles"].strip())
    if m is None or m.GetNumBonds() < 3: continue
    X.append(feats(m)); W.append(Crippen.MolMR(m))
X = np.array(X,float); W = np.array(W); nat, TI = X[:,0], X[:,1:]

rng = np.random.default_rng(11)
cvs, r2b, dr2 = [], [], []
for w in [0,1,2,3,4,6,8,11,15,20,30]:
    for _ in range(400):
        c = rng.choice(nat); pool = np.where(np.abs(nat-c) <= w)[0]
        if len(pool) < 15: continue
        i = rng.choice(pool, 15, replace=False)
        if np.std(W[i]) < 1e-9: continue
        best = max(r2(TI[i,j], W[i]) for j in range(K))
        cvs.append(nat[i].std()/nat[i].mean()); r2b.append(best)
        dr2.append(best - r2(nat[i], W[i]))
cvs, r2b, dr2 = map(np.array, (cvs, r2b, dr2))
bins = np.linspace(0, .55, 12); ctr = .5*(bins[:-1]+bins[1:])
d = np.digitize(cvs, bins)-1
mr, md = [], []
for i in range(len(ctr)):
    m = d==i
    mr.append(r2b[m].mean() if m.sum()>15 else np.nan)
    md.append(dr2[m].mean() if m.sum()>15 else np.nan)

# ---- Tm / Tb for panel (b) ----
pubchem_db.autoload_main_db()
cas2smi = {chemicals.CAS_to_int(o.CASs): o.smiles
           for o in pubchem_db.CAS_index.values() if o.smiles}
def load(fn):
    p = os.path.join(os.path.dirname(chemicals.__file__), 'Phase Change', fn)
    out = {}
    with open(p, encoding='utf-8', errors='replace') as fh:
        next(fh)
        for line in fh:
            q = line.rstrip().split("\t")
            if len(q) > 1 and q[1]:
                try: out[chemicals.CAS_to_int(q[0])] = float(q[1])
                except Exception: pass
    return out
def build(tab):
    F, y = [], []
    for cas, v in tab.items():
        s = cas2smi.get(cas)
        if not s or '.' in s: continue
        m = Chem.MolFromSmiles(s)
        if m is None or not (10 <= m.GetNumAtoms() <= 60) or m.GetNumBonds() < 3: continue
        if max(a.GetDegree() for a in m.GetAtoms()) > 4: continue
        if not any(a.GetSymbol()=='C' for a in m.GetAtoms()): continue
        F.append(feats(m)); y.append(v)
    return np.array(F,float), np.array(y)

fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.6))
ax[0].plot(ctr, mr, "o-", color="#c0392b", lw=2, label=r"reported $R^2$ (best of 14)")
ax[0].plot(ctr, md, "s-", color="#16a085", lw=2, label=r"$\Delta R^2$ over atom-count null")
ax[0].axhline(0, color="k", lw=.7)
ax[0].set_xlabel("size dispersion of the sample,  CV(n)")
ax[0].set_ylabel(r"$R^2$")
ax[0].set_title("(a)  The inversion: as the headline $R^2$ rises,\nthe index's contribution falls  (molar refractivity)", fontsize=10)
ax[0].legend(fontsize=8.5); ax[0].grid(alpha=.25)

col = {"melting point":"#2980b9", "boiling point":"#d35400"}
for lab, fn in [("melting point","OpenNotebook Melting Points.tsv"),
                ("boiling point","Yaws Boiling Points.tsv")]:
    F, y = build(load(fn)); T = F[:,1:]
    rg = np.random.default_rng(5)
    real, perm = [], []
    for _ in range(4000):
        i = rg.choice(len(y),15,replace=False)
        real.append(max(r2(T[i,j], y[i]) for j in range(K)))
        yp = y[rg.choice(len(y),15,replace=False)]
        perm.append(max(r2(T[i,j], yp) for j in range(K)))
    ax[1].hist(real, bins=45, alpha=.50, color=col[lab], density=True,
               label=f"{lab}: real draws (N={len(y)})")
    ax[1].hist(perm, bins=45, histtype="step", lw=1.8, color=col[lab], density=True,
               label=f"{lab}: permutation null")
ax[1].axvline(0.35, color="k", ls=":", lw=1.5)
ax[1].text(0.36, ax[1].get_ylim()[1]*.80,
           "95th pct of the\npermutation null\n"+r"($R^2<0.35$ at $N{=}15$"+"\nis not evidence)", fontsize=7.4)
ax[1].set_xlabel(r"best-of-14 $R^2$ from a draw of 15 molecules")
ax[1].set_ylabel("density")
ax[1].set_title("(b)  Sampling variance of the reported statistic,\nagainst its true permutation null", fontsize=10)
ax[1].legend(fontsize=7.2); ax[1].grid(alpha=.25)
fig.tight_layout()
fig.savefig("/mnt/user-data/outputs/figure2.png", dpi=170)
print("saved")
