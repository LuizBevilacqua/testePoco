import matplotlib.pyplot as plt
import numpy as np
import mpmath as mp
from mpmath import *
from scipy import special as sp
from invlaplace import invlaplace

N=100000
rd = np.array([1, 1.1, 1.2, 1.3, 1.4, 1.5, 2, 3, 5, 10, 15, 20])
tD = np.linspace(10**(-2), 10**3, N)
G = np.zeros((len(tD),len(rd)))
trd = np.zeros((len(tD),len(rd)))
# Definindo s como uma variável simbólicas

for i in range(len(tD)):
    for r in range(len(rd)):
        f = lambda u: sp.k0(rd[r] * np.sqrt(u)) / (u *np.sqrt(u) * sp.k1(np.sqrt(u)))
        var = invlaplace(tD[i],f, 16)
        G[i][r] = var.gavsteh(f,tD[i],16)
        trd[i][r] = tD[i] / rd[r]**2

Gtrans = G.T
trd = trd.T
for i in range(len(rd)):
    plt.loglog(trd[i][:],Gtrans[i][:])


plt.axis([1e-2, 1e3, 1e-2, 10])

plt.show()
