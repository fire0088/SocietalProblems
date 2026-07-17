"""
Renyi HDE via correct coupled integration.
Variables in e-folds N=ln a. Track the dimensionless horizon u = H0 * L.
Future EH: dL/dt = HL - 1  =>  dL/dN = L - 1/H  (H in H0 units, L in c/H0 units)
  => du/dN = u - 1/E   where E=H/H0, u=L/(c/H0)... consistent dimensionless.
Density: rho_DE/rho_c0 = Om_DE0_factor / (u^2 (1+ delta_hat u^2))  where delta_hat=delta*pi*(c/H0)^2
  We fold constants: rho_DE = K/(u^2(1+ d u^2)), d = effective Renyi param (dimensionless).
Friedmann: E^2 = Om_m a^-3 + Om_r a^-4 + rho_DE/rho_c0.
Normalize K so Om_DE(today)=Om_DE0.
Integrate from FUTURE (large a, u->?) backward, OR shoot. Simplest: integrate forward
in N from a small a with u guessed, iterate u boundary so EH is consistent.
Actually EH is a future integral -> integrate BACKWARD from far future where we set
the de Sitter attractor: as a->inf, matter/rad ->0, DE dominates, E->const=E_inf,
and future EH L-> c/H_inf (de Sitter horizon), so u_inf = 1/E_inf.
"""
import numpy as np
from scipy.integrate import solve_ivp

def renyi_E(d_hat, Om_m0, Om_r, zmax=3.0):
    # de Sitter future: E_inf^2 = rho_DE_inf. In dS, u=1/E (horizon=1/H).
    # rho_DE_inf = K/(u^2(1+d u^2)), u=1/E_inf => E_inf^2=K/( (1/E_inf^2)(1+d/E_inf^2))
    #   => E_inf^2 = K E_inf^2/(1+d/E_inf^2) => 1 = K/(1+d/E_inf^2) => E_inf^2=d/(K-1)
    # We'll solve for K by shooting to Om_DE(a=1)=Om_DE0. Parametrize by K, find E_inf.
    Om_DE0=1-Om_m0-Om_r
    def integrate(K):
        if K<=1: return None
        E_inf=np.sqrt(d_hat/(K-1)) if d_hat>0 else np.sqrt(K)  # d=0: dS with rho=K/u^2, u=1/E->E^2=K
        if d_hat==0: E_inf=np.sqrt(Om_DE0)  # not used for d=0 path
        # state: y=[u, ...]; integrate dU/dN backward from N_far to N=ln(a_lo)
        Nf=np.log(300.0); Nlo=np.log(1/(1+zmax))
        def rhs(N,Y):
            u=Y[0]; a=np.exp(N)
            rho_DE=K/(u**2*(1+d_hat*u**2))
            E=np.sqrt(max(Om_m0*a**-3+Om_r*a**-4+rho_DE,1e-12))
            return [u-1.0/E]
        u_far=1.0/np.sqrt(max(Om_m0*np.exp(-3*Nf)+Om_r*np.exp(-4*Nf)+K/((1/ (d_hat>0 and np.sqrt(d_hat/(K-1)) or 1))**2),1e-6)) if False else 1.0/E_inf
        sol=solve_ivp(rhs,[Nf,Nlo],[1.0/E_inf],dense_output=True,rtol=1e-8,atol=1e-11,max_step=0.05)
        return sol
    # shoot on K to hit Om_DE0 at a=1
    from scipy.optimize import brentq
    def resid(K):
        sol=integrate(K)
        if sol is None: return 1e3
        u1=sol.sol(0.0)[0]; a=1.0
        rho_DE=K/(u1**2*(1+d_hat*u1**2))
        E1=np.sqrt(Om_m0+Om_r+rho_DE)
        Om_DE_today=rho_DE/E1**2
        return Om_DE_today-Om_DE0
    try:
        Kbest=brentq(resid,1.001,20,xtol=1e-6)
    except Exception as e:
        return None
    sol=integrate(Kbest)
    Ngrid=np.linspace(0,np.log(1/(1+zmax)),400)
    us=np.array([sol.sol(N)[0] for N in Ngrid]); a=np.exp(Ngrid); z=1/a-1
    rho_DE=Kbest/(us**2*(1+d_hat*us**2))
    E=np.sqrt(Om_m0*a**-3+Om_r*a**-4+rho_DE)
    o=np.argsort(z)
    return z[o],E[o],Kbest

# validate d_hat=0 -> standard HDE w0=-0.891
for d in [0.0, 0.5, 1.0, 2.0]:
    r=renyi_E(d,0.3,9.2e-5)
    if r is None: print(f"  d={d}: FAILED"); continue
    z,E,K=r; i0=np.argmin(np.abs(z))
    Om_DE=E**2-0.3*(1+z)**3-9.2e-5*(1+z)**4
    lnr=np.log(Om_DE/Om_DE[i0]); w=-1-(1/3)*np.gradient(lnr,np.log(1/(1+z)))
    print(f"  Renyi d_hat={d}: K={K:.3f}, w0={w[i0]:.3f}" + (" <- target -0.891" if d==0 else ""))
