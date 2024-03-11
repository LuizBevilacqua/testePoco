import matplotlib.pyplot as plt
import numpy as np
import mpmath as mp
from mpmath import *

mp.dps = 15
mp.pretty = True

N=1000
rd = np.array([1, 1.1, 1.2, 1.3, 1.4, 1.5, 2, 3, 5, 10, 15, 20])
tD= np.linspace(1e-2,1e3,N)
G = np.zeros((len(tD),len(rd)))

# Definindo s como uma variável simbólica

for i in range(len(tD)):
    for r in range(len(rd)):
        f = lambda s: besselk(0, rd[r] * sqrt(s)) / (s ** (3 / 2) * besselk(1, sqrt(s)))
        G[i][r] = mp.invertlaplace(f, tD[i], method='Stehfest', degree=18)

print(G.shape)


