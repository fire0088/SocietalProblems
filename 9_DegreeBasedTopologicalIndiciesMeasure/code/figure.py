import numpy as np, csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
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

# ---------- load ESOL ----------
smiles, logS = [], []
with open("data/delaney-processed.csv") as fh:
    for row in csv.DictReader(fh):
        smiles.append(row["smiles"].strip())
        logS.append(float(row["measured log solubility in mols per litre"]))

X, Y = [], []
for smi, ls in zip(smiles, logS):
    mol = Chem.MolFromSmiles(smi)
    if mol is None or mol.GetNumBonds() < 3: continue
    e = [(b.GetBeginAtom().GetDegree(), b.GetEndAtom().GetDegree()) for b in mol.GetBonds()]
    X.append([mol.GetNumAtoms()] + [sum(f(a,b) for a,b in e) for f in INDICES.values()])
    Y.append([Descriptors.MolWt(mol), Crippen.MolMR(mol), ls])
X = np.array(X, float); Y = np.array(Y, float)
natoms = X[:,0]; TI = X[:,1:]
PROPS = ["MolWt  (additive)", "MolMR  (additive)", "logS  (experimental)"]
print(f"ESOL: N={len(X)}, heavy atoms {natoms.min():.0f}-{natoms.max():.0f}")

def r2(x, y):
    if np.std(x) < 1e-12 or np.std(y) < 1e-12: return 0.0
    return float(np.corrcoef(x, y)[0,1]**2)

# ---------- sweep: subsamples of n=15 (typical paper size) at controlled size dispersion ----------
rng = np.random.default_rng(42)
SAMPLE_N = 15
widths = [0, 1, 2, 3, 4, 6, 8, 11, 15, 20, 30]
res = {p: {"cv": [], "best": [], "null": [], "delta": [], "winner": []} for p in PROPS}

for w in widths:
    for _ in range(400):
        c = rng.choice(natoms)
        pool = np.where(np.abs(natoms - c) <= w)[0]
        if len(pool) < SAMPLE_N: continue
        idx = rng.choice(pool, SAMPLE_N, replace=False)
        cv = natoms[idx].std() / natoms[idx].mean()
        for k, p in enumerate(PROPS):
            y = Y[idx, k]
            if np.std(y) < 1e-9: continue
            sc = [r2(TI[idx, j], y) for j in range(K)]
            b = int(np.argmax(sc)); nul = r2(natoms[idx], y)
            res[p]["cv"].append(cv); res[p]["best"].append(sc[b])
            res[p]["null"].append(nul); res[p]["delta"].append(sc[b]-nul)
            res[p]["winner"].append(INAMES[b])

# ---------- bin and plot ----------
bins = np.array([0, .02, .05, .08, .12, .16, .20, .26, .32, .40, .55])
ctr = 0.5*(bins[:-1]+bins[1:])
fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
COL = {"best":"#c0392b", "null":"#2c3e50"}

summary = {}
for ax, p in zip(axes, PROPS):
    cv = np.array(res[p]["cv"])
    d = np.digitize(cv, bins) - 1
    mb, mn, md, lo, hi = [], [], [], [], []
    for i in range(len(ctr)):
        m = d == i
        if m.sum() < 15: mb.append(np.nan); mn.append(np.nan); md.append(np.nan); lo.append(np.nan); hi.append(np.nan); continue
        b = np.array(res[p]["best"])[m]; nn = np.array(res[p]["null"])[m]
        mb.append(b.mean()); mn.append(nn.mean()); md.append((b-nn).mean())
        lo.append(np.percentile(b,25)); hi.append(np.percentile(b,75))
    summary[p] = (ctr, np.array(mb), np.array(mn), np.array(md))
    ax.fill_between(ctr, lo, hi, color=COL["best"], alpha=.15, lw=0)
    ax.plot(ctr, mb, "o-", color=COL["best"], lw=2, ms=5, label="best of 14 indices")
    ax.plot(ctr, mn, "s--", color=COL["null"], lw=2, ms=5, label="heavy-atom count (null)")
    ax.set_xlabel("size dispersion of the sample,  CV(n)")
    ax.set_ylabel(r"$R^2$"); ax.set_title(p, fontsize=11)
    ax.set_ylim(-.03, 1.03); ax.grid(alpha=.25); ax.legend(fontsize=8, loc="lower right")
    ax.axvspan(0.28, 0.55, color="orange", alpha=.10, lw=0)

axes[0].text(0.40, .12, "range typical of\npublished drug sets",
             fontsize=7.5, ha="center", color="#8a5a00")
fig.suptitle("Reported $R^2$ in degree-based-index QSPR is a function of how much size variation the sample contains\n"
             "(ESOL, subsamples of 15 molecules — the typical published study size)", fontsize=11.5)
fig.tight_layout(rect=[0,0,1,0.90])
fig.savefig("/mnt/user-data/outputs/size_dispersion.png", dpi=170)

# ---------- numbers for the text ----------
print("\n{:<22} {:>8} {:>10} {:>10} {:>10}".format("property","CV(n)","R2 best","R2 null","dR2"))
for p in PROPS:
    c, mb, mn, md = summary[p]
    for i in (1, 5, 9):
        if not np.isnan(mb[i]):
            print(f"{p:<22} {c[i]:8.3f} {mb[i]:10.4f} {mn[i]:10.4f} {md[i]:+10.4f}")

print("\nWinner instability at low vs high size dispersion (logS):")
cv = np.array(res[PROPS[2]]["cv"]); win = np.array(res[PROPS[2]]["winner"])
for lab, m in [("CV < 0.05", cv < .05), ("CV > 0.30", cv > .30)]:
    if m.sum() > 20:
        u, ct = np.unique(win[m], return_counts=True)
        top = sorted(zip(ct,u), reverse=True)[:3]
        print(f"  {lab:<10} {len(u):2d} distinct winners; top: " +
              ", ".join(f"{n} {100*c/m.sum():.0f}%" for c,n in top))

# ---------- the isomer limit, stated exactly ----------
print("\nIsomer limit (CV(n) = 0 exactly): within a constitutional-isomer class")
print("  MolWt, molar refraction, molar volume, polarizability are CONSTANT ->")
print("  R^2 is undefined; these properties carry zero information at fixed n.")
