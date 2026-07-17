"""DESI DR2 BAO, Table IV of arXiv:2503.14738 (Abdul-Karim et al. 2025).
Each tracer: z_eff, DM/rd, sigma(DM/rd), DH/rd, sigma(DH/rd), r (corr DM-DH).
BGS and QSO are DV/rd only (isotropic) -> stored separately.
Values as published (DR2)."""
import numpy as np
# tracer: z, DM/rd, err, DH/rd, err, r_MH
DESI_DR2 = [
    # BGS z=0.295 is DV/rd only: DV/rd=7.944+/-0.075
    ("LRG1", 0.510, 13.588, 0.167, 21.863, 0.425, -0.459),
    ("LRG2", 0.706, 17.351, 0.177, 19.455, 0.330, -0.404),
    ("LRG3+ELG1", 0.934, 21.576, 0.152, 17.641, 0.193, -0.416),
    ("ELG2", 1.321, 27.601, 0.318, 14.176, 0.221, -0.434),
    ("QSO", 1.484, 30.512, 0.760, 12.817, 0.516, -0.500),
    ("Lya", 2.330, 38.988, 0.531, 8.632, 0.101, -0.482),
]
# isotropic DV/rd
DESI_DV = [("BGS",0.295,7.944,0.075)]

def get_arrays():
    z=np.array([t[1] for t in DESI_DR2])
    DM=np.array([t[2] for t in DESI_DR2]); sDM=np.array([t[3] for t in DESI_DR2])
    DH=np.array([t[4] for t in DESI_DR2]); sDH=np.array([t[5] for t in DESI_DR2])
    r =np.array([t[6] for t in DESI_DR2])
    # build block-diagonal covariance (each tracer 2x2 DM,DH correlated by r)
    n=len(z); C=np.zeros((2*n,2*n))
    for i in range(n):
        C[2*i,2*i]=sDM[i]**2; C[2*i+1,2*i+1]=sDH[i]**2
        C[2*i,2*i+1]=r[i]*sDM[i]*sDH[i]; C[2*i+1,2*i]=r[i]*sDM[i]*sDH[i]
    return z,DM,DH,C
if __name__=="__main__":
    z,DM,DH,C=get_arrays()
    print("DESI DR2:",len(z),"anisotropic tracers z=",z)
    print("DM/rd frac err:",np.round([DM[i]/DM[i]*0+ (np.sqrt(C[2*i,2*i])/DM[i]) for i in range(len(z))],4))
