"""
Barrow HDE, future event horizon — clean validated implementation.
ODE: dOmega_DE/dx = (2-Delta) Omega_DE(1-Omega_DE)[1/(2-Delta) + (1/c)Omega_DE^{1/(2-Delta)}]
  reduces to Li(2004) standard HDE at Delta=0. x=ln a.
Then: w_DE = -1 - (1/3) dln rho_DE/dx, rho_DE from Omega_DE and rho_m.
Delta=0 must give w0=-0.891 (c=1,Om_DE0=0.7). Validated.
"""
import numpy as np
from scipy.integrate import solve_ivp

def bhde(Delta, c, Om_m0=0.3, zmax=3.0, N=3000):
    Om_DE0 = 1-Om_m0
    p=1.0/(2.0-Delta)
    def rhs(x,Y):
        Om=min(max(Y[0],1e-12),1-1e-12)
        return [(2-Delta)*Om*(1-Om)*(1.0/(2-Delta)+(1.0/c)*Om**p)]
    xb=np.linspace(0,-np.log(1+zmax),N)
    sb=solve_ivp(rhs,[0,xb[-1]],[Om_DE0],t_eval=xb,rtol=1e-10,atol=1e-12)
    x=sb.t[::-1]; Om_DE=sb.y[0][::-1]; z=np.exp(-x)-1   # ascending z
    rho_m=Om_m0*(1+z)**3
    rho_DE=rho_m*Om_DE/(1-Om_DE)
    i0=np.argmin(np.abs(z))
    rr=rho_DE/rho_DE[i0]
    E=np.sqrt((rho_m+rho_DE)/(rho_m[i0]+rho_DE[i0]))
    w=-1-(1/3)*np.gradient(np.log(rr),x)
    return z,w,rr,E

print("Barrow HDE validation & scan (c=1, Om_m0=0.3):")
print(f"{'Delta':>6} {'w0':>8} {'w(z=1)':>8} {'rhoDE(1)/0':>11}")
for Delta in [0.0,0.1,0.2,0.4,0.6,0.8]:
    z,w,rr,E=bhde(Delta,1.0)
    i0=np.argmin(np.abs(z)); i1=np.argmin(np.abs(z-1))
    print(f"{Delta:6.1f} {w[i0]:+8.3f} {w[i1]:+8.3f} {rr[i1]:11.3f}")
print("\n[target] Delta=0 w0 should be -0.891")

# save the Delta=0.6 case shape for inspection
z,w,rr,E=bhde(0.6,1.0)
np.savez('/home/claude/bhde_shape.npz',z=z,w=w,rr=rr,E=E)
