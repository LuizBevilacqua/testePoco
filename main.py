import matplotlib.pyplot as plt
import numpy as np
import mpmath as mp
from mpmath import *
from scipy import special as sp
from invlaplace import invlaplace

#mp.dps = 100
#mp.pretty = True

N=1000
rd = np.array([1, 1.1, 1.2, 1.3, 1.4, 1.5, 2, 3, 5, 10, 15, 20])
tD= np.linspace(1e-2,1e3,1000)
G = np.zeros((len(tD),len(rd)))
trd = np.zeros(len(tD))
# Definindo s como uma variável simbólica
for i in range(len(tD)):
    print(i)
    for r in range(len(rd)):
        f = lambda s: sp.k0(rd[r] * np.sqrt(s)) / (s ** (3 / 2) * sp.k1(np.sqrt(s)))
        var = invlaplace(tD[i],f, 8)
        G[i][r] = var.gavsteh_param()
        trd[i] = tD[i]/(rd[r]**2)

pdPlot = np.zeros((len(rd),len(tD)))
pd0 = [linha[0] for linha in G]
pd1 = [linha[1] for linha in G]
pd2 = [linha[2] for linha in G]
pd3 = [linha[3] for linha in G]
pd4 = [linha[4] for linha in G]
pd5 = [linha[5] for linha in G]
pd6 = [linha[6] for linha in G]
pd7 = [linha[7] for linha in G]
pd8 = [linha[8] for linha in G]
pd9 = [linha[9] for linha in G]
pd10 = [linha[10] for linha in G]
pd11 = [linha[11] for linha in G]

plt.loglog(trd,pd0)
plt.loglog(trd,pd1)
plt.loglog(trd,pd2)
plt.loglog(trd,pd3)
plt.loglog(trd,pd4)
plt.loglog(trd,pd5)
plt.loglog(trd,pd6)
plt.loglog(trd,pd7)
plt.loglog(trd,pd8)
plt.loglog(trd,pd9)
plt.loglog(trd,pd10)
plt.loglog(trd,pd11)

plt.show()