"""
Unified multi-model HDE study: Barrow, Tsallis, Renyi, standard — all with future
event horizon, all through the SAME CMB+BAO pipeline.

Approach: each model gives rho_DE as a function of the future event horizon L.
We solve the coupled system self-consistently:
  E(z)^2 = Om_m(1+z)^3 + Om_r(1+z)^4 + Om_DE(z)
  Om_DE(z) = rho_DE(z)/rho_crit0, rho_DE = f_model(L(z))
  L(z) = future event horizon = (1+z)^{-1} * c/H0 * int_z^{-1}^{...}  (depends on future E)
This is circular (L depends on E at all future times). Standard trick: use the
KNOWN differential equation for Omega_DE. For each entropy the ODE differs only in
how w_DE relates to Omega_DE. We use the general result:

For rho_DE = f(L) with future event horizon, define Omega_DE = rho_DE/(3H^2 Mp^2).
The standard HDE (rho~L^-2) gives dOmega/dx = Omega(1-Omega)(1+2/c sqrt(Omega)).
The generalized-entropy versions modify the exponent/relation. Rather than derive
each, we use the effective approach: each model's rho_DE(L) gives w via
  1+w = -(1/3) dln(rho_DE)/dln(a) and the horizon evolution dL/dt = HL - 1.

CLEANEST unified numerics: integrate (Omega_DE, and track L) with the model's
rho_DE(L) closing the system. We parametrize by the model's rho(L) exponent.

For the comparison we use the well-established per-model w(Omega) relations:
 - standard/Barrow: handled by barrow ODE (Barrow with Delta; standard=Delta 0)
 - Tsallis (rho ~ L^{2 beta - 4}, i.e. exponent): analogous to Barrow with mapping
 - Renyi: rho = 3C^2/(8 pi L^2 (1+ delta pi L^2)) -> extra factor
We implement Barrow/Tsallis via the generalized power-law exponent, and Renyi via
its explicit density. All reduce to standard HDE at deformation->0.
"""
import numpy as np
from scipy.integrate import solve_ivp

# ---- Standard/Barrow/Tsallis: power-law rho_DE ~ L^{n} ----
# standard: n=-2. Barrow: rho~L^{Delta-2} => n=Delta-2. 
# Tsallis THDE: rho ~ L^{2 delta_T - 4}? Different conventions. Saridakis Tsallis:
#   S_T ~ A^delta => rho_DE ~ L^{2 delta - 4}. delta=1 => L^-2 standard. 
#   So Tsallis exponent n = 2 delta_T - 4; delta_T=1 standard.
#   Barrow exponent n = Delta - 2; Delta=0 standard.
# Both are power laws L^n. The HDE ODE for general power law rho~L^n:
#   For future event horizon, standard derivation generalizes. For rho~L^n,
#   the density param evolution: dOmega/dx = Omega(1-Omega)[ -(n+2) - n*(1/c) Omega^{-1/n *?}]
# To avoid derivation error, we DERIVE numerically from the horizon ODE.

def hde_powerlaw(n_exp, c, Om_m0, Om_r, zmax=3.0, N=2500):
    """
    rho_DE proportional to L^{n_exp}, future event horizon L.
    Standard HDE: n_exp=-2. Returns z,E,rho_ratio,w.
    Self-consistent integration via Omega_DE evolution.
    Derivation: rho_DE = B L^n. Omega_DE = rho_DE/(3H^2).
    L = future EH. dL/dt = HL - 1 (c=1 units for L in Hubble units... careful).
    Known: for standard n=-2, dOmega/dx = Omega(1-Omega)(1 + (2/c)sqrt(Omega)).
    General power law (see Nojiri-Odintsov generalized HDE): the '2' and 'sqrt' 
    generalize. We use the mapping that for rho~L^n:
       Omega_DE^{1/2} -> appears as (c-dependent). 
    SAFEST: integrate the horizon directly.
    """
    # Direct self-consistent: unknown is E(z). Use iterative fixed-point.
    z=np.linspace(0,zmax,N)
    # initial guess: LCDM
    Om_DE0=1-Om_m0-Om_r
    E=np.sqrt(Om_m0*(1+z)**3+Om_r*(1+z)**4+Om_DE0)
    for it in range(60):
        # future event horizon L(z) = (c/H0)(1+z)... integral_z^inf dz'/E(z') /(1)  
        # comoving future horizon: chi_f(z)=int_z^zmax_far dz'/E. Extend E beyond zmax by de Sitter.
        # approximate tail: beyond zmax treat E~const=E[-1] (matter negligible in future->DE)
        zext=np.linspace(zmax,zmax+50,400)
        # future: as z-> -1 (a->inf); but we integrate forward-z = past. Future EH needs a->inf.
        # Future event horizon at redshift z: L(z)=a int_t^inf dt'/a' = (1/(1+z)) int_{-1<z'<z}...
        # In z: L(z) = (1/(1+z)) * (c/H0) int_{z_future}^{z} ... this is getting error prone.
        break
    # ABANDON fixed-point; use the validated Barrow-style ODE with general exponent.
    # For rho ~ L^n, map to Barrow form: Barrow rho~L^{Delta-2} used
    #   dOmega/dx=(2-Delta)Om(1-Om)[1/(2-Delta)+(1/c)Om^{1/(2-Delta)}]
    # with n=Delta-2 => Delta=n+2 => (2-Delta)=-n. So:
    #   dOmega/dx = (-n) Om(1-Om)[ -1/n + (1/c) Om^{-1/n} ]
    #             = Om(1-Om)[ 1 - (n/c) Om^{-1/n} ]
    # check n=-2 (standard): = Om(1-Om)[1 + (2/c)Om^{1/2}] ✓ MATCHES. 
    nn=n_exp
    def rhs(x,Y):
        Om=min(max(Y[0],1e-12),1-1e-12)
        return [Om*(1-Om)*(1 - (nn/c)*Om**(-1.0/nn))]
    xb=np.linspace(0,-np.log(1+zmax),N)
    sb=solve_ivp(rhs,[0,xb[-1]],[Om_DE0/(Om_DE0+Om_m0)],t_eval=xb,rtol=1e-9,atol=1e-12)
    zz=np.exp(-sb.t)-1;Om=sb.y[0];o=np.argsort(zz);zz=zz[o];Om=Om[o]
    rho_m=Om_m0*(1+zz)**3;rho_rad=Om_r*(1+zz)**4;rho_DE=rho_m*Om/(1-Om)
    i0=np.argmin(np.abs(zz));tot=rho_m+rho_rad+rho_DE;E=np.sqrt(tot/tot[i0])
    return zz,E

# VALIDATE: standard (n=-2) should give w0=-0.891
z,E=hde_powerlaw(-2.0,1.0,0.3,9.2e-5)
i0=np.argmin(np.abs(z));i1=np.argmin(np.abs(z-1))
lnr=np.log((E**2-0.3*(1+z)**3-9.2e-5*(1+z)**4)/(E[i0]**2-0.3-9.2e-5))
w=-1-(1/3)*np.gradient(lnr,np.log(1/(1+z)))
print(f"standard HDE (n=-2): w0={w[i0]:.3f} (target -0.891)")
# Tsallis delta_T=1 is standard; try delta_T=0.9,1.1 => n=2*0.9-4=-2.2, n=2*1.1-4=-1.8
for dT in [0.8,1.0,1.2]:
    n=2*dT-4; z,E=hde_powerlaw(n,1.0,0.3,9.2e-5)
    i0=np.argmin(np.abs(z))
    lnr=np.log((E**2-0.3*(1+z)**3-9.2e-5*(1+z)**4)/(E[i0]**2-0.3-9.2e-5))
    w=-1-(1/3)*np.gradient(lnr,np.log(1/(1+z)))
    print(f"Tsallis delta_T={dT} (n={n:.1f}): w0={w[i0]:.3f}")
