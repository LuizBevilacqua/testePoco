import matplotlib.pyplot as plt
import numpy as np
import mpmath as mp
from mpmath import *
from scipy import special as sp
from invlaplace import invlaplace

N=10000
red = np.array([100, 200, 300, 400, 500, 600,10000])
tD = np.linspace(1e3, 1e6, N)
G = np.zeros((len(tD),len(red)))
rD = 1
# Definindo s como uma variável simbólicas

for i in range(len(tD)):
    for r in range(len(red)):
        f = lambda u: (
                (sp.i0(rD * np.sqrt(u)) * sp.k1(red[r] * np.sqrt(u)))
                + (sp.k0(rD * np.sqrt(u)) * sp.i1(red[r] * np.sqrt(u)))
            ) / (
                (
                    (u ** (3 / 2))
                    * (
                        sp.i1(red[r] * np.sqrt(u)) * sp.k1(np.sqrt(u))
                        - (sp.i1(np.sqrt(u)) * sp.k1(red[r] * np.sqrt(u)))
                    )
                )
            )
        var = invlaplace(tD[i],f, 10)
        G[i][r] = var.gavsteh(f,tD[i],10)

Gtrans = G.T

fig ,ax = plt.subplots()

for i in range(len(red)):
    if i == len(red)-1:
        ax.semilogx(tD,Gtrans[i][:],label=r'$\frac{r_e}{r_w} \rightarrow \infty$')
    else:    
        ax.semilogx(tD, Gtrans[i][:], label=r'$\frac{re}{rw} =' + f'{red[i]}$')
    

plt.legend()
plt.xlabel(r'$\mathbf{log(t_D)}$')
plt.ylabel(r'$\mathbf{p_D}$')
plt.grid(True, which="both", ls="-.")
plt.axis([1e3,1e6,4,11])
plt.show()
