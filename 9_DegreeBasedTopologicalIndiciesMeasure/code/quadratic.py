import os, csv, numpy as np, chemicals
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

def poly_r2(x, y, deg):
    """in-sample R^2 of a degree-`deg` polynomial fit, numerically stabilised"""
    if np.std(x) < 1e-12 or np.std(y) < 1e-12: return 0.0
    xs = (x - x.mean())/x.std()
    V = np.vander(xs, deg+1)
    try: b, *_ = np.linalg.lstsq(V, y, rcond=None)
    except np.linalg.LinAlgError: return 0.0
    res = y - V @ b
    ss = np.sum((y - y.mean())**2)
    return float(max(0.0, 1 - np.sum(res**2)/ss)) if ss > 0 else 0.0

def poly_cv_r2(x, y, deg, reps=200, frac=.7, seed=0):
    rg = np.random.default_rng(seed); N = len(y); sc = []
    xs = (x - x.mean())/(x.std()+1e-12)
    for _ in range(reps):
        p = rg.permutation(N); tr, te = p[:int(frac*N)], p[int(frac*N):]
        V = np.vander(xs[tr], deg+1)
        b, *_ = np.linalg.lstsq(V, y[tr], rcond=None)
        pr = np.vander(xs[te], deg+1) @ b
        ss = np.sum((y[te]-y[te].mean())**2)
        if ss > 0: sc.append(1 - np.sum((y[te]-pr)**2)/ss)
    return float(np.mean(sc))

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
def featurize(mol):
    e=[(b.GetBeginAtom().GetDegree(),b.GetEndAtom().GetDegree()) for b in mol.GetBonds()]
    return [mol.GetNumAtoms()]+[sum(f(a,b) for a,b in e) for f in INDICES.values()]
def build(tab):
    F,y,ms=[],[],[]
    for cas,v in tab.items():
        s=cas2smi.get(cas)
        if not s or '.' in s: continue
        m=Chem.MolFromSmiles(s)
        if m is None or not (10<=m.GetNumAtoms()<=60) or m.GetNumBonds()<3: continue
        if not any(a.GetSymbol()=='C' for a in m.GetAtoms()): continue
        F.append(featurize(m)); y.append(v); ms.append(m)
    return np.array(F,float), np.array(y), ms

DATA = {}
Ftm, ytm, mtm = build(load('OpenNotebook Melting Points.tsv'))
DATA["melting point (measured)"] = (Ftm, ytm)
DATA["molecular weight (additive)"] = (Ftm, np.array([Descriptors.MolWt(q) for q in mtm]))
Ftb, ytb, _ = build(load('Yaws Boiling Points.tsv'))
DATA["boiling point (measured)"] = (Ftb, ytb)

print("="*78)
print("A.  FULL SAMPLE -- does a higher-order fit change the null-baseline verdict?")
print("="*78)
print(f"{'property':<30}{'deg':>4}{'R2(n)':>9}{'R2(best)':>10}{'dR2':>9}{'CV R2(n)':>10}{'CV best':>9}")
for name,(F,y) in DATA.items():
    n, TI = F[:,0], F[:,1:]
    for deg in (1,2,3):
        rn = poly_r2(n,y,deg)
        sc = [poly_r2(TI[:,j],y,deg) for j in range(K)]
        b = int(np.argmax(sc))
        cn = poly_cv_r2(n,y,deg); cb = poly_cv_r2(TI[:,b],y,deg)
        print(f"{name if deg==1 else '':<30}{deg:>4}{rn:9.4f}{sc[b]:10.4f}{sc[b]-rn:+9.4f}{cn:10.4f}{cb:9.4f}")

print("\n"+"="*78)
print("B.  N=15 -- the published regime.  In-sample R2, as papers report it.")
print("="*78)
print(f"{'property':<26}{'deg':>4}{'real mean':>11}{'real 95th':>11}{'NULL mean':>11}{'NULL 95th':>11}")
thresholds = {}
for name in ["melting point (measured)","boiling point (measured)"]:
    F,y = DATA[name]; TI = F[:,1:]
    for deg in (1,2,3):
        rg = np.random.default_rng(5); real, perm = [], []
        for _ in range(3000):
            i = rg.choice(len(y),15,replace=False)
            real.append(max(poly_r2(TI[i,j],y[i],deg) for j in range(K)))
            yp = y[rg.choice(len(y),15,replace=False)]
            perm.append(max(poly_r2(TI[i,j],yp,deg) for j in range(K)))
        real,perm = np.array(real),np.array(perm)
        thresholds[(name,deg)] = np.percentile(perm,95)
        print(f"{name if deg==1 else '':<26}{deg:>4}{real.mean():11.3f}"
              f"{np.percentile(real,95):11.3f}{perm.mean():11.3f}{np.percentile(perm,95):11.3f}")

print("\n"+"="*78)
print("C.  Where the published r=0.836 (R2=0.698) boiling-point claim now sits")
print("="*78)
for deg in (1,2,3):
    t = thresholds[("boiling point (measured)",deg)]
    verdict = "INSIDE the null -- not evidence" if 0.698 <= t else "above the null"
    print(f"   degree {deg}: permutation-null 95th pct = {t:.3f}   ->  {verdict}")
