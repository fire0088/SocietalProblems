import os, numpy as np, chemicals, csv
from chemicals.identifiers import pubchem_db
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
RDLogger.DisableLog('rdApp.*')

INDICES = {
 "M1":lambda a,b:a+b, "M2":lambda a,b:a*b, "F":lambda a,b:a*a+b*b,
 "HM":lambda a,b:(a+b)**2, "Randic":lambda a,b:1/np.sqrt(a*b),
 "sum-conn":lambda a,b:1/np.sqrt(a+b), "ABC":lambda a,b:np.sqrt((a+b-2)/(a*b)),
 "GA":lambda a,b:2*np.sqrt(a*b)/(a+b), "harmonic":lambda a,b:2/(a+b),
 "Sombor":lambda a,b:np.sqrt(a*a+b*b), "AZI":lambda a,b:(a*b/max(a+b-2,1e-9))**3,
 "ISI":lambda a,b:a*b/(a+b), "SDD":lambda a,b:a/b+b/a, "ReZG3":lambda a,b:a*b*(a+b)}
INAMES = list(INDICES); K = len(INAMES)

def r2(x,y):
    if np.std(x)<1e-12 or np.std(y)<1e-12: return 0.0
    return float(np.corrcoef(x,y)[0,1]**2)

pubchem_db.autoload_main_db()
cas2smi = {chemicals.CAS_to_int(o.CASs): o.smiles
           for o in pubchem_db.CAS_index.values() if o.smiles}
def load(fn):
    p = os.path.join(os.path.dirname(chemicals.__file__),'Phase Change',fn); out={}
    with open(p,encoding='utf-8',errors='replace') as fh:
        next(fh)
        for line in fh:
            q=line.rstrip().split("\t")
            if len(q)>1 and q[1]:
                try: out[chemicals.CAS_to_int(q[0])]=float(q[1])
                except Exception: pass
    return out
def build(tab, maxdeg=4):
    F,y,mm = [],[],[]
    for cas,v in tab.items():
        s=cas2smi.get(cas)
        if not s or '.' in s: continue
        m=Chem.MolFromSmiles(s)
        if m is None or not (10<=m.GetNumAtoms()<=60) or m.GetNumBonds()<3: continue
        if max(a.GetDegree() for a in m.GetAtoms())>maxdeg: continue
        if not any(a.GetSymbol()=='C' for a in m.GetAtoms()): continue
        e=[(b.GetBeginAtom().GetDegree(),b.GetEndAtom().GetDegree()) for b in m.GetBonds()]
        F.append([m.GetNumAtoms(), len(e)]+[sum(f(a,b) for a,b in e) for f in INDICES.values()])
        y.append(v); mm.append(m)
    return np.array(F,float), np.array(y), mm

print("="*72)
print("CHECK 1 -- is Figure 2b actually a NULL distribution?")
print("="*72)
for lab,fn in [("melting point","OpenNotebook Melting Points.tsv"),
               ("boiling point","Yaws Boiling Points.tsv")]:
    F,y,_ = build(load(fn)); TI = F[:,2:]
    rg = np.random.default_rng(5)
    real, perm = [], []
    for _ in range(4000):
        i = rg.choice(len(y),15,replace=False)
        real.append(max(r2(TI[i,j],y[i]) for j in range(K)))
        yp = y[rg.choice(len(y),15,replace=False)]      # response from unrelated molecules
        perm.append(max(r2(TI[i,j],yp) for j in range(K)))
    real,perm = np.array(real),np.array(perm)
    print(f"\n{lab}  (full-sample R2 best-of-14 = {max(r2(TI[:,j],y) for j in range(K)):.4f})")
    print(f"   REAL data,  N=15 draws : mean {real.mean():.3f}  95th {np.percentile(real,95):.3f}")
    print(f"   PERMUTED null, N=15    : mean {perm.mean():.3f}  95th {np.percentile(perm,95):.3f}")
    print(f"   -> published R2=0.698 is at percentile {100*(real<0.698).mean():.1f} of REAL,"
          f" {100*(perm<0.698).mean():.1f} of NULL")

print("\n"+"="*72)
print("CHECK 2 -- correlation-transfer: does r(TI,m) bound |r(TI,P) - r(m,P)| for ANY P?")
print("="*72)
print("Lemma: for centred unit vectors, |r(X,Z)-r(Y,Z)| <= ||x-y|| = sqrt(2(1-r(X,Y)))")
F,y,mols = build(load('OpenNotebook Melting Points.tsv'))
n, m, TI = F[:,0], F[:,1], F[:,2:]
props = {"Tm (measured)": y,
         "MolWt": np.array([Descriptors.MolWt(q) for q in mols]),
         "MolMR": np.array([Crippen.MolMR(q) for q in mols]),
         "TPSA":  np.array([rdMolDescriptors.CalcTPSA(q) for q in mols]),
         "cLogP": np.array([Crippen.MolLogP(q) for q in mols])}
print(f"\n{'index':<10}{'r(TI,m)':>9}{'bound':>9}   " + "".join(f"{k:>16}" for k in props))
worst = 0
for j,nm in enumerate(INAMES):
    rho = np.corrcoef(TI[:,j], m)[0,1]
    bnd = np.sqrt(max(2*(1-rho),0))
    cells = []
    for k,P in props.items():
        d = abs(np.corrcoef(TI[:,j],P)[0,1] - np.corrcoef(m,P)[0,1])
        worst = max(worst, d-bnd); cells.append(f"{d:16.4f}")
    print(f"{nm:<10}{rho:9.4f}{bnd:9.4f}   " + "".join(cells))
print(f"\nlargest violation of the bound across all index x property pairs: {worst:+.5f}")

print("\n"+"="*72)
print("CHECK 3 -- is heavy-atom count a fair null, or a strawman? try stronger nulls")
print("="*72)
rings = np.array([rdMolDescriptors.CalcNumRings(q) for q in mols], float)
MW = props["MolWt"]
def cv_r2(F_, yy, reps=100, frac=.7, seed=0):
    rg=np.random.default_rng(seed); N=len(yy); Fm=np.column_stack([F_,np.ones(N)]); sc=[]
    for _ in range(reps):
        p=rg.permutation(N); tr,te=p[:int(frac*N)],p[int(frac*N):]
        b,*_=np.linalg.lstsq(Fm[tr],yy[tr],rcond=None); pr=Fm[te]@b
        sc.append(1-np.sum((yy[te]-pr)**2)/np.sum((yy[te]-yy[te].mean())**2))
    return np.mean(sc)
best_j = int(np.argmax([r2(TI[:,j],y) for j in range(K)]))
print(f"predicting measured Tm (N={len(y)}), cross-validated R2:")
for lab,Fm in [("heavy-atom count n", n[:,None]),
               ("molecular weight",   MW[:,None]),
               ("n + ring count",     np.column_stack([n,rings])),
               (f"best index ({INAMES[best_j]})", TI[:,[best_j]]),
               ("all 14 indices",     TI),
               ("all 14 + n + rings", np.column_stack([TI,n,rings]))]:
    print(f"   {lab:<28} {cv_r2(Fm,y):.4f}")

print("\n"+"="*72)
print("CHECK 4 -- does the maxdeg<=4 filter bias the sample?")
print("="*72)
F4,y4,_ = build(load('OpenNotebook Melting Points.tsv'), maxdeg=4)
F6,y6,_ = build(load('OpenNotebook Melting Points.tsv'), maxdeg=6)
print(f"   maxdeg<=4: N={len(y4)}   maxdeg<=6: N={len(y6)}   excluded {len(y6)-len(y4)}"
      f" ({100*(len(y6)-len(y4))/len(y6):.1f}%)")
for lab,(Fx,yx) in [("maxdeg<=4",(F4,y4)),("maxdeg<=6",(F6,y6))]:
    T=Fx[:,2:]
    print(f"   {lab}: R2(n)={r2(Fx[:,0],yx):.4f}  R2(best)={max(r2(T[:,j],yx) for j in range(K)):.4f}")
