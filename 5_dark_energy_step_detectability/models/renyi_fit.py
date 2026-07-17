import numpy as np, sys
exec(open('multimodel_final.py').read().split('print("FINAL')[0])
# Renyi only, coarser grid for speed
print("Renyi HDE fits (coarse grid):")
for dv in [0.3,0.5,1.0,2.0]:
    b,bp=best('renyi',dv,H0range=np.linspace(65,74,5))
    print(f"  Renyi d={dv}: total={b:.1f} at Om={bp[0]:.3f} H0={bp[1]:.1f} c={bp[2]:.2f}")
