import numpy as np, itertools
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
from drugs import DRUGS

DRUGS = {k: v for k, v in DRUGS.items() if k not in ("ritonavir", "lopinavir")}

# ---- degree-based topological indices: TI = sum over edges of f(du,dv) ----
INDICES = {
 "M1 (1st Zagreb)":      lambda a,b: a+b,
 "M2 (2nd Zagreb)":      lambda a,b: a*b,
 "F (forgotten)":        lambda a,b: a*a+b*b,
 "HM (hyper-Zagreb)":    lambda a,b: (a+b)**2,
 "Randic":               lambda a,b: 1/np.sqrt(a*b),
 "sum-connectivity":     lambda a,b: 1/np.sqrt(a+b),
 "ABC":                  lambda a,b: np.sqrt((a+b-2)/(a*b)),
 "GA (geom-arith)":      lambda a,b: 2*np.sqrt(a*b)/(a+b),
 "harmonic":             lambda a,b: 2/(a+b),
 "Sombor":               lambda a,b: np.sqrt(a*a+b*b),
 "AZI (augmented)":      lambda a,b: (a*b/(a+b-2))**3,
 "ISI (inv sum indeg)":  lambda a,b: a*b/(a+b),
 "SDD (symm div deg)":   lambda a,b: a/b + b/a,
 "redefined ReZG3":      lambda a,b: a*b*(a+b),
}

names, edgeprofiles, rows, props = [], [], [], []
for name, smi in DRUGS.items():
    mol = Chem.MolFromSmiles(smi)              # hydrogen-suppressed graph
    degs = [at.GetDegree() for at in mol.GetAtoms()]
    edges = [(b.GetBeginAtom().GetDegree(), b.GetEndAtom().GetDegree())
             for b in mol.GetBonds()]
    n, m = mol.GetNumAtoms(), mol.GetNumBonds()
    vals = [sum(f(a,b) for a,b in edges) for f in INDICES.values()]
    names.append(name)
    rows.append([n, m] + vals)
    edgeprofiles.append(edges)
    props.append([Descriptors.MolWt(mol), Crippen.MolMR(mol),
                  rdMolDescriptors.CalcLabuteASA(mol), rdMolDescriptors.CalcTPSA(mol)])

X = np.array(rows, float)          # cols: n, m, then the 14 indices
P = np.array(props, float)
cols = ["n (heavy atoms)", "m (bonds)"] + list(INDICES)
propnames = ["MolWt", "MolMR (molar refractivity)", "LabuteASA (surface area)", "TPSA"]
N = len(names)
print(f"N = {N} drugs;  heavy atoms {X[:,0].min():.0f}-{X[:,0].max():.0f}, "
      f"CV(m) = {X[:,1].std()/X[:,1].mean():.3f}\n")

# ---- 1. realized edge-degree-pair profile ----
from collections import Counter
tot = Counter()
for edges in edgeprofiles:
    for a,b in edges: tot[tuple(sorted((a,b)))] += 1
T = sum(tot.values())
print("Realized edge-degree-pair profile across all molecules:")
for pair, c in sorted(tot.items(), key=lambda kv: -kv[1]):
    print(f"   {pair}: {100*c/T:5.2f}%")

# ---- 2. worst-case vs realized eta for each index ----
allpairs = [(a,b) for a in range(1,5) for b in range(a,5) if (a,b)!=(1,1)]
realized = [p for p in tot if tot[p]/T > 0.001]
print("\n{:<22} {:>10} {:>10}   {:>10}".format("index","eta_worst","eta_real","r(TI,m)"))
for j,(nm,f) in enumerate(INDICES.items()):
    fw = [f(a,b) for a,b in allpairs]
    fr = [f(a,b) for a,b in realized]
    ew = (max(fw)-min(fw))/(max(fw)+min(fw))
    er = (max(fr)-min(fr))/(max(fr)+min(fr))
    r  = np.corrcoef(X[:,2+j], X[:,1])[0,1]
    print(f"{nm:<22} {ew:10.3f} {er:10.3f}   {r:10.5f}")

# ---- 3. collinearity among the indices ----
TI = X[:,2:]
C = np.corrcoef(TI.T)
off = C[np.triu_indices(len(INDICES),1)]
print(f"\nPairwise |r| among the {len(INDICES)} indices: "
      f"min {np.abs(off).min():.4f}, median {np.median(np.abs(off)):.4f}, max {np.abs(off).max():.4f}")
print(f"  fraction of pairs with |r| > 0.99: {np.mean(np.abs(off)>0.99):.2%}")

Z = (TI - TI.mean(0)) / TI.std(0)
sv = np.linalg.svd(Z, compute_uv=False)
ev = sv**2 / (sv**2).sum()
print(f"  PCA variance explained: PC1 {ev[0]:.4%}, PC2 {ev[1]:.4%}, PC3 {ev[2]:.4%}")

# ---- 4. null baseline: does heavy-atom count beat the 'winning' index? ----
def r2(x, y):
    r = np.corrcoef(x, y)[0,1]; return r*r
print("\n{:<28} {:>10} {:>10} {:>12} {:>10}".format(
      "property","R2(n only)","best TI R2","best TI","gain"))
for k, pn in enumerate(propnames):
    y = P[:,k]
    base = r2(X[:,0], y)
    scores = [(r2(TI[:,j], y), nm) for j,nm in enumerate(INDICES)]
    best, bestnm = max(scores)
    print(f"{pn:<28} {base:10.4f} {best:10.4f} {bestnm:>12} {best-base:+10.4f}")
    lo = min(s for s,_ in scores)
    print(f"{'':<28} {'':>10} (worst index R2 = {lo:.4f})")
