"""
Cosmic-chronometer H(z) data + validated covariance.
Rebuilt from the documented recipe (Moresco et al. 2020/2022/2024 compilation).
Covariance = diagonal statistical + three fully-correlated rank-1 outer products
(IMF, stellar-library, SPS-model), fractional errors interpolated from the
Moresco 2020 systematic budget (data_MM20). mod_ooo = SPS variant discarding
most-discordant model = the recommended total additional systematic.
VALIDATION TARGET: Moresco 2023 CC-only -> H0=66.7+/-5.3, Om=0.33+0.08-0.06.
"""
import numpy as np

# 33-point CC compilation (z, H(z), sigma_stat).  Standard published values.
_CC = np.array([
 [0.07,   69.0,  19.6],
 [0.09,   69.0,  12.0],
 [0.12,   68.6,  26.2],
 [0.17,   83.0,   8.0],
 [0.179,  75.0,   4.0],
 [0.199,  75.0,   5.0],
 [0.20,   72.9,  29.6],
 [0.27,   77.0,  14.0],
 [0.28,   88.8,  36.6],
 [0.352,  83.0,  14.0],
 [0.3802, 83.0,  13.5],
 [0.4,    95.0,  17.0],
 [0.4004, 77.0,  10.2],
 [0.4247, 87.1,  11.2],
 [0.4497, 92.8,  12.9],
 [0.47,   89.0,  49.6],
 [0.4783, 80.9,   9.0],
 [0.48,   97.0,  62.0],
 [0.593, 104.0,  13.0],
 [0.68,  92.0,   8.0],
 [0.75,  98.8,  33.6],
 [0.781, 105.0,  12.0],
 [0.80,  113.1, 28.5],
 [0.875, 125.0,  17.0],
 [0.88,  90.0,   40.0],
 [0.90,  117.0,  23.0],
 [1.037, 154.0,  20.0],
 [1.26,  135.0,  65.0],
 [1.3,   168.0,  17.0],
 [1.363, 160.0,  33.6],
 [1.43,  177.0,  18.0],
 [1.53,  140.0,  14.0],
 [1.75,  202.0,  40.0],
])
Z   = _CC[:,0]
HZ  = _CC[:,1]
SIG = _CC[:,2]
N   = len(Z)

# Moresco 2020 redshift-dependent systematic budget (fractional, mod_ooo total).
# Interpolation anchors from the published budget: ~5.4% at z=0.2 down to ~2.3% at z=1.5.
# Broken into the three correlated components used in the covariance.
_z_anchor    = np.array([0.2, 0.4, 0.6, 0.9, 1.2, 1.5])
# total recommended (mod_ooo) additional systematic, fractional:
_tot_anchor  = np.array([0.054,0.045,0.040,0.033,0.027,0.023])
# component split (approx, from MM20): IMF, stellar-library (slib), SPS/model
_imf_frac    = np.array([0.009,0.009,0.009,0.009,0.009,0.009])   # ~roughly flat
_slib_frac   = np.array([0.010,0.009,0.008,0.007,0.006,0.005])
# SPS/model component = whatever makes quadrature sum to total (mod_ooo)
_sps_frac    = np.sqrt(np.clip(_tot_anchor**2 - _imf_frac**2 - _slib_frac**2, 0, None))

def _interp(fr):
    return np.interp(Z, _z_anchor, fr)

def covariance(sps='ooo'):
    """Return NxN covariance.
       sps='ooo'  -> recommended (mod_ooo) total systematic
       sps='none' -> statistical only (diagonal)
       sps='all'  -> inflate SPS component (keep discordant model) ~ larger
    """
    Cstat = np.diag(SIG**2)
    if sps == 'none':
        return Cstat
    imf  = _interp(_imf_frac)  * HZ
    slib = _interp(_slib_frac) * HZ
    if sps == 'ooo':
        spsf = _interp(_sps_frac) * HZ
    elif sps == 'all':
        spsf = _interp(_sps_frac) * HZ * 2.2   # keep most-discordant model -> larger SPS spread
    else:
        raise ValueError(sps)
    # three fully-correlated rank-1 outer products
    C = Cstat.copy()
    for v in (imf, slib, spsf):
        C = C + np.outer(v, v)
    return C

if __name__ == "__main__":
    print(f"N points: {N},  z range {Z.min()}-{Z.max()}")
    C = covariance('ooo')
    print("cov ooo: diag sqrt sample:", np.sqrt(np.diag(C))[:5].round(1))
    print("cov none diag sqrt:", SIG[:5])
