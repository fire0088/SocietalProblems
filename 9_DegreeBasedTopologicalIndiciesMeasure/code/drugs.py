from rdkit import Chem
from rdkit.Chem import Descriptors

# Drugs that recur across the QSPR topological-index literature:
# anticancer, antiviral/COVID, and common small-molecule therapeutics.
DRUGS = {
 "aspirin":            "CC(=O)Oc1ccccc1C(=O)O",
 "paracetamol":        "CC(=O)Nc1ccc(O)cc1",
 "ibuprofen":          "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
 "naproxen":           "COc1ccc2cc(ccc2c1)C(C)C(=O)O",
 "caffeine":           "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
 "metformin":          "CN(C)C(=N)NC(N)=N",
 "warfarin":           "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O",
 "diclofenac":         "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl",
 "indomethacin":       "COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1",
 "celecoxib":          "Cc1ccc(cc1)-c1cc(nn1-c1ccc(cc1)S(N)(=O)=O)C(F)(F)F",
 "atorvastatin":       "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O",
 "simvastatin":        "CCC(C)(C)C(=O)O[C@H]1C[C@H](C)C=C2C=C[C@H](C)[C@H](CC[C@@H]3C[C@@H](O)CC(=O)O3)[C@@H]12",
 "omeprazole":         "COc1ccc2[nH]c(nc2c1)S(=O)Cc1ncc(C)c(OC)c1C",
 "fluoxetine":         "CNCCC(Oc1ccc(cc1)C(F)(F)F)c1ccccc1",
 "sertraline":         "CN[C@H]1CC[C@@H](c2ccc(Cl)c(Cl)c2)c2ccccc21",
 "diazepam":           "CN1c2ccc(Cl)cc2C(=NCC1=O)c1ccccc1",
 "morphine":           "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5",
 "codeine":            "COc1ccc2C[C@H]3N(C)CC[C@@]45[C@@H](Oc1c24)[C@H](O)C=C[C@H]35",
 # antivirals / COVID
 "favipiravir":        "NC(=O)c1nc(F)cnc1O",
 "molnupiravir":       "CC(C)C(=O)O[C@@H]1[C@H](O)[C@@H](CO)O[C@H]1N1C=CC(=NO)NC1=O",
 "remdesivir":         "CCC(CC)COC(=O)[C@H](C)N[P@](=O)(OC[C@H]1O[C@](C#N)([C@H](O)[C@@H]1O)c1ccc2c(N)ncnn12)Oc1ccccc1",
 "chloroquine":        "CCN(CC)CCCC(C)Nc1ccnc2cc(Cl)ccc12",
 "hydroxychloroquine": "CCN(CCO)CCCC(C)Nc1ccnc2cc(Cl)ccc12",
 "ritonavir":          "CC(C)c1nc(CN(C)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](Cc2ccccc2)C[C@H](O)[C@H](Cc2ccccc2)NC(=O)OCc2cncs2)cs1",
 "lopinavir":          "CC(C)[C@H](NC(=O)N1CCCNC1=O)C(=O)N[C@H](C[C@H](O)[C@H](Cc1ccccc1)NC(=O)COc1c(C)cccc1C)Cc1ccccc1",
 "oseltamivir":        "CCOC(=O)C1=C[C@@H](OC(CC)CC)[C@H](NC(C)=O)[C@@H](N)C1",
 "acyclovir":          "Nc1nc2n(COCCO)cnc2c(=O)[nH]1",
 "zidovudine":         "Cc1cn([C@H]2C[C@H](N=[N+]=[N-])[C@@H](CO)O2)c(=O)[nH]c1=O",
 # anticancer
 "gemcitabine":        "NC1=NC(=O)N(C=C1)[C@H]1O[C@H](CO)[C@@H](O)C1(F)F",
 "5-fluorouracil":     "O=c1[nH]cc(F)c(=O)[nH]1",
 "capecitabine":       "CCCCCOC(=O)Nc1nc(=O)n([C@@H]2O[C@H](C)[C@@H](O)[C@H]2O)cc1F",
 "methotrexate":       "CN(Cc1cnc2nc(N)nc(N)c2n1)c1ccc(cc1)C(=O)N[C@@H](CCC(=O)O)C(=O)O",
 "imatinib":           "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(n1)-c1cccnc1",
 "erlotinib":          "COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC",
 "gefitinib":          "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1",
 "tamoxifen":          "CC/C(=C(\\c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1",
 "topotecan":          "CC[C@@]1(O)C(=O)OCc2c1cc1n(c2=O)Cc2cc3c(CN(C)C)c(O)ccc3nc2-1",
 "irinotecan":         "CC[C@@]1(O)C(=O)OCc2c1cc1n(c2=O)Cc2cc3c(CC)c(OC(=O)N4CCC(CC4)N4CCCCC4)ccc3nc2-1",
 "doxorubicin":        "COc1cccc2c1C(=O)c1c(O)c3c(c(O)c1C2=O)C[C@@](O)(C(=O)CO)C[C@@H]3O[C@H]1C[C@H](N)[C@H](O)[C@H](C)O1",
 "cyclophosphamide":   "ClCCN(CCCl)P1(=O)NCCCO1",
}

known_mw = {"aspirin":180.2,"paracetamol":151.2,"ibuprofen":206.3,"caffeine":194.2,
            "metformin":129.2,"warfarin":308.3,"diclofenac":296.1,"celecoxib":381.4,
            "atorvastatin":558.6,"omeprazole":345.4,"fluoxetine":309.3,"diazepam":284.7,
            "morphine":285.3,"favipiravir":157.1,"remdesivir":602.6,"chloroquine":319.9,
            "ritonavir":720.9,"lopinavir":628.8,"oseltamivir":312.4,"acyclovir":225.2,
            "gemcitabine":263.2,"5-fluorouracil":130.1,"methotrexate":454.4,
            "imatinib":493.6,"erlotinib":393.4,"gefitinib":446.9,"tamoxifen":371.5,
            "topotecan":421.4,"irinotecan":586.7,"doxorubicin":543.5,"zidovudine":267.2}

bad = []
for name, smi in DRUGS.items():
    m = Chem.MolFromSmiles(smi)
    if m is None:
        bad.append((name, "UNPARSEABLE")); continue
    mw = Descriptors.MolWt(m)
    if name in known_mw and abs(mw - known_mw[name]) > 1.0:
        bad.append((name, f"MW {mw:.1f} vs expected {known_mw[name]}"))
print(f"{len(DRUGS)} molecules; {len(bad)} problems")
for b in bad: print("  ", b)
