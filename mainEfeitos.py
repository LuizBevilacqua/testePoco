import matplotlib.pyplot as plt
import numpy as np
import mpmath as mp
from mpmath import *
from scipy import special as sp
from invlaplace import invlaplace

N=10000
red = 10
tD = np.linspace(1e2, 1e8, N)
G = np.zeros((len(tD),len(red)))
rD = 1
# Definindo s como uma variável simbólicas
cD = [1e2,1e3,1e4,1e5]
s = [-5, 0,  5, 10, 20]
for i in range(len(tD)):
    for j in range(len(cD)):
        for k in range(len(s)):
            f = lambda u: (
                (sp.k0(red * np.sqrt(u)))/
                u*(np.sqrt(u)*sp.k1(np.sqrt(u)) + cD[j]*u*(sp.k0(np.sqrt(u)) + s[k]*np.sqrt(u) * sp.k1(np.sqrt(u))))
            )
            #ACHAR UM JEITO DE GUARDAR PARA CADA VALOR DE S
            var = invlaplace(tD[i],f, 10)
            G[i][j] = var.gavsteh(f,tD[i],10)

Gtrans = G.T

fig ,ax = plt.subplots()

for i in range(len(red)):
    ax.semilogx(tD,Gtrans[i][:])

plt.grid(True, which="both", ls="-.")
plt.axis([1e2,1e8,1e-1,1e2])
plt.show()
