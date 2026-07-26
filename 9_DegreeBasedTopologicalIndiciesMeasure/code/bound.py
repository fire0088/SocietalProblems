import os, numpy as np, chemicals
from chemicals.identifiers import pubchem_db
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
from nonadditive import INDICES, INAMES, K, cas2smi, Tm

PAIRS = [(1,1),(1,2),(1,3),(1,4),(2,2),(2,3),(2,4),(3,3),(3,4),(4,4)]

def profiles(mols):
    Pm, ms = [], []
    for mol in mols:
        e = [tuple(sorted((b.GetBeginAtom().GetDegree(), b.GetEndAtom().GetDegree())))
             for b in mol.GetBonds()]
        m = len(e)
        Pm.append([sum(1 for x in e if x == pr)/m for pr in PAIRS]); ms.append(m)
    return np.array(Pm), np.array(ms, float)

mols = []
for cas in Tm:
    smi = cas2smi.get(cas)
    if not smi or '.' in smi: continue
    mol = Chem.MolFromSmiles(smi)
    if mol is None or not (10 <= mol.GetNumAtoms() <= 60) or mol.GetNumBonds() < 3: continue
    if max(a.GetDegree() for a in mol.GetAtoms()) > 4: continue
    if not any(a.GetSymbol()=='C' for a in mol.GetAtoms()): continue
    mols.append(mol)
Pm, ms = profiles(mols)
N = len(ms); CV = ms.std()/ms.mean()
print(f"N = {N} drug-like molecules; CV(m) = {CV:.3f}")

w = ms**2 / (ms**2).sum()                      # m^2 weights (E_i = m_i<u_i,f>)
pbar = Pm.mean(0)
U = Pm - pbar
Ubar = np.sqrt(np.sum(w * (np.abs(U).sum(1)**2)))     # m-weighted rms of ||u||_1
Sig = (U * w[:,None]).T @ U                          # m^2-weighted covariance
supp = pbar > 1e-4
print(f"profile support: {[PAIRS[i] for i in range(len(PAIRS)) if supp[i]]}")
print(f"m-weighted rms ||p_i - p_bar||_1 = {Ubar:.3f}")

def centered(v): return v - v.mean()
print("\n{:<10} {:>7} {:>8} {:>9} {:>9} {:>9} {:>9}".format(
      "index","omega","g_bar","tau_exact","tau_Hold","tau_spec","r2 actual"))
rows = []
for j,(nm,f) in enumerate(INDICES.items()):
    fv = np.array([f(a,b) if (a,b)!=(1,1) else 0.0 for a,b in PAIRS])
    fs = fv[supp]
    omega = 0.5*(fs.max()-fs.min())
    gbar = float(pbar @ fv)
    TI = ms * (Pm @ fv)
    Ev  = ms * (U @ fv)
    tau_ex   = np.linalg.norm(centered(Ev)) / (gbar*np.linalg.norm(centered(ms)))
    tau_hold = (omega/gbar) * Ubar * np.sqrt(1+CV**2)/CV
    fc = fv - 0.5*(fs.max()+fs.min())
    tau_spec = np.sqrt(max(fc @ Sig @ fc, 0))*np.sqrt(N) / (gbar*np.linalg.norm(centered(ms))) \
               * np.sqrt(np.sum(ms**2))/np.sqrt(N)
    r2a = np.corrcoef(TI, ms)[0,1]**2
    rows.append((nm, omega, gbar, tau_ex, tau_hold, tau_spec, r2a))
    print(f"{nm:<10} {omega:7.2f} {gbar:8.3f} {tau_ex:9.4f} {tau_hold:9.3f} {tau_spec:9.4f} {r2a:9.5f}")

print("\nBound quality (1 - tau^2 vs actual r^2; a valid bound must be <= actual):")
ok_h = ok_s = 0
for nm,om,g,te,th,ts,r2a in rows:
    bh, bs = max(0,1-th**2), max(0,1-ts**2)
    ok_h += bh <= r2a + 1e-9; ok_s += bs <= r2a + 1e-9
    print(f"  {nm:<10} Holder {bh:7.4f} {'ok' if bh<=r2a else 'VIOLATED'} | "
          f"spectral {bs:7.4f} {'ok' if bs<=r2a else 'VIOLATED'} | actual {r2a:.5f}"
          f" | spectral slack {r2a-bs:+.4f}")
print(f"\nvalid: Holder {ok_h}/{K}, spectral {ok_s}/{K}")
print(f"Holder non-vacuous (tau<1) for {sum(1 for r in rows if r[4]<1)}/{K} indices")
print(f"spectral non-vacuous for {sum(1 for r in rows if r[5]<1)}/{K} indices")
