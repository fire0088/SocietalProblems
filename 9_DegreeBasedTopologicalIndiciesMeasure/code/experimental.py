import numpy as np, csv, sys
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
RDLogger.DisableLog('rdApp.*')

INDICES = {
 "M1":lambda a,b:a+b, "M2":lambda a,b:a*b, "F":lambda a,b:a*a+b*b,
 "HM":lambda a,b:(a+b)**2, "Randic":lambda a,b:1/np.sqrt(a*b),
 "sum-conn":lambda a,b:1/np.sqrt(a+b), "ABC":lambda a,b:np.sqrt((a+b-2)/(a*b)),
 "GA":lambda a,b:2*np.sqrt(a*b)/(a+b), "harmonic":lambda a,b:2/(a+b),
 "Sombor":lambda a,b:np.sqrt(a*a+b*b), "AZI":lambda a,b:(a*b/max(a+b-2,1e-9))**3,
 "ISI":lambda a,b:a*b/(a+b), "SDD":lambda a,b:a/b+b/a, "ReZG3":lambda a,b:a*b*(a+b),
}
INAMES = list(INDICES)

def featurize(smiles):
    keep, feats = [], []
    for i, smi in enumerate(smiles):
        mol = Chem.MolFromSmiles(smi)
        if mol is None or mol.GetNumBonds() < 3: continue
        edges = [(b.GetBeginAtom().GetDegree(), b.GetEndAtom().GetDegree())
                 for b in mol.GetBonds()]
        row = [mol.GetNumAtoms(), len(edges)] + \
              [sum(f(a,b) for a,b in edges) for f in INDICES.values()]
        feats.append(row); keep.append(i)
    return np.array(keep), np.array(feats, float)

def r2(x, y):
    if np.std(x) == 0: return 0.0
    return np.corrcoef(x, y)[0,1]**2

def heldout_rmse(x, y, seed=0, reps=200, frac=0.7):
    rng = np.random.default_rng(seed); n = len(y); errs = []
    for _ in range(reps):
        p = rng.permutation(n); tr, te = p[:int(frac*n)], p[int(frac*n):]
        if np.std(x[tr]) == 0: continue
        b, a = np.polyfit(x[tr], y[tr], 1)
        errs.append(np.mean((y[te] - (b*x[te]+a))**2))
    return np.sqrt(np.mean(errs))

# ================= PART 1: ZINC drug-like, large N =================
smis = []
with open("data/250k_rndm_zinc_drugs_clean_3.csv") as fh:
    rd = csv.DictReader(fh)
    for i, row in enumerate(rd):
        if i >= 20000: break
        smis.append(row["smiles"].strip())
rng = np.random.default_rng(1)
smis = [smis[i] for i in rng.choice(len(smis), 6000, replace=False)]
_, Z = featurize(smis)
print(f"PART 1 - ZINC drug-like sample: N = {len(Z)}, "
      f"heavy atoms {Z[:,0].min():.0f}-{Z[:,0].max():.0f}, CV(m) = {Z[:,1].std()/Z[:,1].mean():.3f}")
TIz = Z[:,2:]
S = (TIz - TIz.mean(0))/TIz.std(0)
ev = np.linalg.svd(S, compute_uv=False)**2; ev /= ev.sum()
C = np.corrcoef(TIz.T); off = np.abs(C[np.triu_indices(len(INAMES),1)])
print(f"  PC1 {ev[0]:.3%}  PC2 {ev[1]:.3%}  PC3 {ev[2]:.3%}")
print(f"  pairwise |r|: median {np.median(off):.4f}, min {off.min():.4f}")
print(f"  r(TI, bond count m): min {min(np.corrcoef(TIz[:,j],Z[:,1])[0,1] for j in range(len(INAMES))):.4f}")

# ================= PART 2: ESOL, computed vs experimental =================
ids, smiles, exp = [], [], []
with open("data/delaney-processed.csv") as fh:
    for row in csv.DictReader(fh):
        smiles.append(row["smiles"].strip())
        exp.append(float(row["measured log solubility in mols per litre"]))
keep, X = featurize(smiles)
exp = np.array(exp)[keep]
mols = [Chem.MolFromSmiles(smiles[i]) for i in keep]

PROPS = {
 "MolWt  [additive]":        np.array([Descriptors.MolWt(m) for m in mols]),
 "MolMR  [additive/Crippen]":np.array([Crippen.MolMR(m) for m in mols]),
 "LabuteASA [additive]":     np.array([rdMolDescriptors.CalcLabuteASA(m) for m in mols]),
 "logS  [EXPERIMENTAL]":     exp,
}
print(f"\nPART 2 - ESOL: N = {len(exp)}, heavy atoms {X[:,0].min():.0f}-{X[:,0].max():.0f}")
print("\n{:<28} {:>9} {:>9} {:>10} {:>9} {:>18}".format(
      "property","R2(n)","R2(best)","dR2","best idx","held-out RMSE n/best"))
for pn, y in PROPS.items():
    base = r2(X[:,0], y)
    sc = [(r2(X[:,2+j], y), INAMES[j]) for j in range(len(INAMES))]
    best, bnm = max(sc)
    rb = heldout_rmse(X[:,0], y); rt = heldout_rmse(X[:,2+INAMES.index(bnm)], y)
    print(f"{pn:<28} {base:9.4f} {best:9.4f} {best-base:+10.4f} {bnm:>9} "
          f"{rb:8.3f} /{rt:8.3f}  ({100*(rb-rt)/rb:+.1f}%)")

# ---- how much of the winner's edge is selection over 14 collinear predictors? ----
print("\nSelection bias check on logS: shuffle y, take max R2 over the 14 indices")
null_max, null_n = [], []
rr = np.random.default_rng(7)
for _ in range(2000):
    ys = exp[rr.permutation(len(exp))]
    null_max.append(max(r2(X[:,2+j], ys) for j in range(len(INAMES))))
    null_n.append(r2(X[:,0], ys))
print(f"  null max-over-14 R2: mean {np.mean(null_max):.5f}, 95th pct {np.percentile(null_max,95):.5f}")
print(f"  null single-predictor R2: mean {np.mean(null_n):.5f}")
print(f"  observed dR2 (best index - n) = {max(r2(X[:,2+j],exp) for j in range(len(INAMES))) - r2(X[:,0],exp):+.5f}")

# ---- multivariate: does adding ALL indices to n help on experimental logS? ----
def cv_r2(F, y, seed=0, reps=200, frac=0.7):
    rng = np.random.default_rng(seed); n = len(y); sc = []
    F = np.column_stack([F, np.ones(len(y))])
    for _ in range(reps):
        p = rng.permutation(n); tr, te = p[:int(frac*n)], p[int(frac*n):]
        beta, *_ = np.linalg.lstsq(F[tr], y[tr], rcond=None)
        pred = F[te] @ beta
        sc.append(1 - np.sum((y[te]-pred)**2)/np.sum((y[te]-y[te].mean())**2))
    return np.mean(sc)
print("\nCross-validated R2 on EXPERIMENTAL logS (70/30, 200 splits):")
print(f"  n alone                        {cv_r2(X[:,[0]], exp):.4f}")
print(f"  best single index              {cv_r2(X[:,[2+INAMES.index(max((r2(X[:,2+j],exp),INAMES[j]) for j in range(len(INAMES)))[1])]], exp):.4f}")
print(f"  all 14 indices                 {cv_r2(X[:,2:], exp):.4f}")
print(f"  n + all 14 indices             {cv_r2(X[:,[0]+list(range(2,2+len(INAMES)))], exp):.4f}")
