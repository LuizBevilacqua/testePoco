import matplotlib.pyplot as plt
import numpy as np
import mpmath as mp
from mpmath import *
from scipy import special as sp
from invlaplace import invlaplace

N=10000
red = 10
tD = np.linspace(1e2, 1e8, N)
rD = 1
# Definindo s como uma variável simbólicas
cD = [0, 1e2, 1e3, 1e4, 1e5]
s = [-1, 0,  5, 10, 20]
pws = np.array([24,36,47,58,70,81,92,103,114,215,307,389,464,531,592,648,698,744,1048,1172,1232,1266,1288,1304,1316,1326,1335,1349,1370,1386,1413])
pws = pws*0.008
deltaT = np.array([0.0109,0.0164,0.0218,0.0273,0.0328,0.0382,0.0437,0.0491,0.0546,0.109,0.164,0.218,0.273,0.328,0.382,0.437,0.491,0.546,1.09,1.64,2.18,2.73,3.28,3.82,4.37,4.91,5.46,6.55,8.74,10.9,16.4])
deltaT = deltaT*20000

G = np.zeros((len(tD),len(s),len(cD)))
for i in range(len(tD)):
    for j in range(len(cD)):
        for k in range(len(s)):
            f = lambda u: (
                (
                    sp.k0(np.sqrt(u)) + s[k]*np.sqrt(u)*sp.k1(np.sqrt(u))
                ) /
                (
                    u*(
                        np.sqrt(u)*sp.k1(np.sqrt(u)) + cD[j]*u*(
                            sp.k0(np.sqrt(u)) + s[k]*np.sqrt(u)*sp.k1(np.sqrt(u))
                        )
                    )
                )
            )
            var = invlaplace(tD[i],f, 12)
            G[i][j][k] = var.gavsteh(f,tD[i],12)

fig ,ax = plt.subplots()

Gc0s0 = [G[t][0][0] for t in range(len(G))]
Gc0s1 = [G[t][0][1] for t in range(len(G))]
Gc0s2 = [G[t][0][2] for t in range(len(G))]
Gc0s3 = [G[t][0][3] for t in range(len(G))]
Gc0s4 = [G[t][0][4] for t in range(len(G))]

Gc1s0 = [G[t][1][0] for t in range(len(G))]
Gc1s1 = [G[t][1][1] for t in range(len(G))]
Gc1s2 = [G[t][1][2] for t in range(len(G))]
Gc1s3 = [G[t][1][3] for t in range(len(G))]
Gc1s4 = [G[t][1][4] for t in range(len(G))]

Gc2s0 = [G[t][2][0] for t in range(len(G))]
Gc2s1 = [G[t][2][1] for t in range(len(G))]
Gc2s2 = [G[t][2][2] for t in range(len(G))]
Gc2s3 = [G[t][2][3] for t in range(len(G))]
Gc2s4 = [G[t][2][4] for t in range(len(G))]

Gc3s0 = [G[t][3][0] for t in range(len(G))]
Gc3s1 = [G[t][3][1] for t in range(len(G))]
Gc3s2 = [G[t][3][2] for t in range(len(G))]
Gc3s3 = [G[t][3][3] for t in range(len(G))]
Gc3s4 = [G[t][3][4] for t in range(len(G))]

Gc4s0 = [G[t][4][0] for t in range(len(G))]
Gc4s1 = [G[t][4][1] for t in range(len(G))]
Gc4s2 = [G[t][4][2] for t in range(len(G))]
Gc4s3 = [G[t][4][3] for t in range(len(G))]
Gc4s4 = [G[t][4][4] for t in range(len(G))]

ax.loglog(tD,Gc0s0, label =r'$C_D = $' + f'{cD[0]}' , color='blue')
ax.loglog(tD,Gc0s1, color='blue')
ax.loglog(tD,Gc0s2, color='blue')
ax.loglog(tD,Gc0s3, color='blue')
ax.loglog(tD,Gc0s4, color='blue')

ax.loglog(tD,Gc1s0, label =r'$C_D = $' + f'{cD[1]}', color='green')
ax.loglog(tD,Gc1s1, color='green')
ax.loglog(tD,Gc1s2, color='green')
ax.loglog(tD,Gc1s3, color='green')
ax.loglog(tD,Gc1s4, color='green')

ax.loglog(tD,Gc2s0, label =r'$C_D = $' + f'{cD[2]}', color='red')
ax.loglog(tD,Gc2s1, color='red')
ax.loglog(tD,Gc2s2, color='red')
ax.loglog(tD,Gc2s3, color='red')
ax.loglog(tD,Gc2s4, color='red')

ax.loglog(tD,Gc3s0, label =r'$C_D = $' + f'{cD[3]}', color='orange')
ax.loglog(tD,Gc3s1, color='orange') 
ax.loglog(tD,Gc3s2, color='orange')
ax.loglog(tD,Gc3s3, color='orange')
ax.loglog(tD,Gc3s4, color='orange')

ax.loglog(tD,Gc4s0, label =r'$C_D = $' + f'{cD[4]}', color='black')
ax.loglog(tD,Gc4s1, color='black')
ax.loglog(tD,Gc4s2, color='black')
ax.loglog(tD,Gc4s3, color='black')
ax.loglog(tD,Gc4s4, color='black')

ax.loglog(deltaT,pws,'o',color='purple',label='Dados de Produção', markersize=3)


plt.xlabel(r'$\mathbf{log  t_D}$')
plt.ylabel(r'$\mathbf{log  p_D}$')
plt.legend(title = r"$C_D$ value:")
plt.grid(True, which="both", ls="-.")
plt.axis([1e2,1e8,1e-3,1e2])
plt.show()
