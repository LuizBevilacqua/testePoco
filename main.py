import matplotlib.pyplot as plt
import numpy as np
import mpmath as mp
from mpmath import *

mp.dps = 30
mp.pretty = True

N=1000
rd = np.array([1, 1.1, 1.2, 1.3, 1.4, 1.5, 2, 3, 5, 10, 15, 20])
tD= np.linspace(1e-2,1e3,50)
G = np.zeros((len(tD),len(rd)))
trd = np.zeros(len(tD))
# Definindo s como uma variável simbólica
for r in range(len(rd)):
    print(r)
    for i in range(len(tD)):
        f = lambda s: besselk(0, rd[r] * sqrt(s)) / (s ** (3 / 2) * besselk(1, sqrt(s)))
        G[i][r] = mp.invertlaplace(f, tD[i])
        trd[i] = tD[i]/(rd[r]**2)

print(G[0][0])
# plt.loglog(trd,G[:][:])
# plt.show()