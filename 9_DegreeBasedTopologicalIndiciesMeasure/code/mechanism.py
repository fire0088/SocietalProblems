import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
from drugs import DRUGS
from analyze import INDICES  # reuse definitions

DRUGS = {k: v for k, v in DRUGS.items() if k not in ("ritonavir", "lopinavir")}
PAIRS = [(1,2),(1,3),(1,4),(2,2),(2,3),(2,4),(3,3),(3,4),(4,4)]

profiles, ms, TIs, props = [], [], [], []
for name, smi in DRUGS.items():
    mol = Chem.MolFromSmiles(smi)
    edges = [tuple(sorted((b.GetBeginAtom().GetDegree(), b.GetEndAtom().GetDegree())))
             for b in mol.GetBonds()]
    m = len(edges)
    p = np.array([sum(1 for e in edges if e == pr) for pr in PAIRS], float) / m
    profiles.append(p); ms.append(m)
    TIs.append([sum(f(a,b) for a,b in edges) for f in INDICES.values()])
    props.append([Descriptors.MolWt(mol), Crippen.MolMR(mol),
                  rdMolDescriptors.CalcLabuteASA(mol), rdMolDescriptors.CalcTPSA(mol)])

Pf = np.array(profiles); ms = np.array(ms, float)
TI = np.array(TIs); PR = np.array(props); N = len(ms)

# ---- mechanism: per-molecule profile concentration ----
pbar = Pf.mean(0)
U = Pf - pbar
l1 = np.abs(U).sum(1)
print("Mean edge-degree-pair profile p-bar:")
for pr, v in zip(PAIRS, pbar): print(f"   {pr}: {v:6.3f}")
print(f"\nPer-molecule deviation ||p_i - p_bar||_1: "
      f"mean {l1.mean():.3f}, median {np.median(l1):.3f}, max {l1.max():.3f}")

print("\n{:<22} {:>9} {:>9} {:>9} {:>9}".format(
      "index","eta_real","tau_bnd","1-tau^2","r2 actual"))
for j,(nm,f) in enumerate(INDICES.items()):
    fv = np.array([f(a,b) for a,b in PAIRS])
    fmin, fmax = fv[pbar>0].min(), fv[pbar>0].max()
    eta = (fmax-fmin)/(fmax+fmin)
    g = Pf @ fv                      # intensive part TI/m
    tau = np.std(ms*(g-g.mean()))/(g.mean()*np.std(ms))
    r2a = np.corrcoef(TI[:,j], ms)[0,1]**2
    print(f"{nm:<22} {eta:9.3f} {tau:9.3f} {max(0,1-tau**2):9.4f} {r2a:9.4f}")

# ---- winner stability under bootstrap ----
def r2(x,y):
    r = np.corrcoef(x,y)[0,1]; return r*r
rng = np.random.default_rng(0)
names = list(INDICES)
print("\nBootstrap stability of the 'best predictor' claim (2000 resamples):")
for k, pn in enumerate(["MolWt","MolMR","LabuteASA","TPSA"]):
    wins = {}
    for _ in range(2000):
        idx = rng.integers(0, N, N)
        s = [r2(TI[idx,j], PR[idx,k]) for j in range(len(names))]
        wins[names[int(np.argmax(s))]] = wins.get(names[int(np.argmax(s))],0)+1
    top = sorted(wins.items(), key=lambda kv:-kv[1])[:3]
    share = ", ".join(f"{n} {100*c/2000:.0f}%" for n,c in top)
    print(f"   {pn:<12} distinct winners: {len(wins):2d}   top-3: {share}")

# ---- leave-one-out: does the best index beat heavy-atom count out of sample? ----
n_atoms = np.array([Chem.MolFromSmiles(s).GetNumAtoms() for s in DRUGS.values()], float)
def loo_rmse(x, y):
    e = []
    for i in range(N):
        msk = np.arange(N) != i
        b, a = np.polyfit(x[msk], y[msk], 1)
        e.append(y[i] - (b*x[i]+a))
    return np.sqrt(np.mean(np.square(e)))
print("\nLeave-one-out RMSE (lower is better):")
for k, pn in enumerate(["MolWt","MolMR","LabuteASA","TPSA"]):
    base = loo_rmse(n_atoms, PR[:,k])
    best = min((loo_rmse(TI[:,j], PR[:,k]), names[j]) for j in range(len(names)))
    print(f"   {pn:<12} n only {base:8.3f} | best index {best[0]:8.3f} ({best[1]})"
          f"  -> {100*(base-best[0])/base:+.1f}%")
