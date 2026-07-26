import os, numpy as np, chemicals
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
 "ISI":lambda a,b:a*b/(a+b), "SDD":lambda a,b:a/b+b/a, "ReZG3":lambda a,b:a*b*(a+b),
}
INAMES = list(INDICES); K = len(INAMES)

pubchem_db.autoload_main_db()
cas2smi = {}
for obj in pubchem_db.CAS_index.values():
    if obj.smiles: cas2smi[chemicals.CAS_to_int(obj.CASs)] = obj.smiles

def load(fname, col):
    d = os.path.join(os.path.dirname(chemicals.__file__), 'Phase Change', fname)
    out = {}
    with open(d, encoding='utf-8', errors='replace') as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2 or not p[1]: continue
            try: out[chemicals.CAS_to_int(p[0])] = float(p[1])
            except Exception: pass
    return out

Tm = load('OpenNotebook Melting Points.tsv', 'Tm')
Tb = load('Yaws Boiling Points.tsv', 'Tb')
print(f"melting points {len(Tm)}, boiling points {len(Tb)}, CAS->SMILES {len(cas2smi)}")

def build(table, lo=10, hi=60):
    """drug-like organic molecules only: heavy atoms in [lo,hi], must contain carbon"""
    rows, ys, mols = [], [], []
    for cas, val in table.items():
        smi = cas2smi.get(cas)
        if not smi: continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None: continue
        na = mol.GetNumAtoms()
        if na < lo or na > hi or mol.GetNumBonds() < 3: continue
        if not any(a.GetSymbol() == 'C' for a in mol.GetAtoms()): continue
        if '.' in smi: continue                      # drop salts / mixtures
        e = [(b.GetBeginAtom().GetDegree(), b.GetEndAtom().GetDegree()) for b in mol.GetBonds()]
        if max(max(p) for p in e) > 4: continue      # keep degrees in {1..4}
        rows.append([na, len(e)] + [sum(f(a,b) for a,b in e) for f in INDICES.values()])
        ys.append(val); mols.append(mol)
    return np.array(rows,float), np.array(ys,float), mols

def r2(x,y):
    if np.std(x)<1e-12 or np.std(y)<1e-12: return 0.0
    return float(np.corrcoef(x,y)[0,1]**2)

for label, tab in [("MELTING POINT Tm (experimental)", Tm), ("BOILING POINT Tb (experimental)", Tb)]:
    X, y, mols = build(tab)
    n, TI = X[:,0], X[:,2:]
    add = {"MolWt (additive)": np.array([Descriptors.MolWt(m) for m in mols]),
           "MolMR (additive)": np.array([Crippen.MolMR(m) for m in mols])}
    print(f"\n=== {label} ===  N={len(y)}, heavy atoms {n.min():.0f}-{n.max():.0f}, CV(n)={n.std()/n.mean():.3f}")
    for pn, yy in list(add.items()) + [(label.split()[0]+" [NON-ADDITIVE]", y)]:
        base = r2(n, yy)
        sc = [r2(TI[:,j], yy) for j in range(K)]
        b = int(np.argmax(sc))
        print(f"   {pn:<28} R2(n)={base:.4f}  R2(best)={sc[b]:.4f} [{INAMES[b]}]  dR2={sc[b]-base:+.4f}")

    # small-sample replication of a published-style study
    rng = np.random.default_rng(3); wins = {}; r2s = []
    for _ in range(3000):
        idx = rng.choice(len(y), 15, replace=False)
        sc = [r2(TI[idx,j], y[idx]) for j in range(K)]
        b = int(np.argmax(sc)); wins[INAMES[b]] = wins.get(INAMES[b],0)+1; r2s.append(sc[b])
    top = sorted(wins.items(), key=lambda kv:-kv[1])[:3]
    print(f"   N=15 published-style draws: {len(wins)} distinct winners; "
          + ", ".join(f"{a} {100*c/3000:.0f}%" for a,c in top))
    print(f"   reported R2 would be {np.mean(r2s):.3f} on average "
          f"(5th-95th pct: {np.percentile(r2s,5):.3f}-{np.percentile(r2s,95):.3f})")
