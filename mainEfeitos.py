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
s = [-5, 0,  5, 10, 20]
G = np.zeros((len(tD),len(s),len(cD)))
for i in range(len(tD)):
    for j in range(len(cD)):
        for k in range(len(s)):
            # f = lambda u: (
            #     (sp.k0(np.sqrt(u)) + s[k]*np.sqrt(u)*sp.k1(np.sqrt(u)))  /
            #     (
            #         u*(
            #             np.sqrt(u)*sp.k1(np.sqrt(u)) + cD[j]*u*(
            #                 sp.k0(np.sqrt(u)) + s[k]*np.sqrt(u)*sp.k1(np.sqrt(u))
            #             )
            #         )
            #     )
            # )

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
            var = invlaplace(tD[i],f, 4)
            G[i][j][k] = var.gavsteh(f,tD[i],4)

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

ax.loglog(tD,Gc0s0, label =f'C = {cD[0]}' , color='blue')
ax.loglog(tD,Gc0s1, color='blue') 
ax.loglog(tD,Gc0s2, color='blue')
ax.loglog(tD,Gc0s3, color='blue')
ax.loglog(tD,Gc0s4, color='blue')

ax.loglog(tD,Gc1s0, label =f'C = {cD[1]}', color='green')
ax.loglog(tD,Gc1s1, color='green') 
ax.loglog(tD,Gc1s2, color='green')
ax.loglog(tD,Gc1s3, color='green')
ax.loglog(tD,Gc1s4, color='green')

ax.loglog(tD,Gc2s0, label =f'C = {cD[2]}', color='red')
ax.loglog(tD,Gc2s1, color='red') 
ax.loglog(tD,Gc2s2, color='red')
ax.loglog(tD,Gc2s3, color='red')
ax.loglog(tD,Gc2s4, color='red')

ax.loglog(tD,Gc3s0, label =f'C = {cD[3]}', color='orange')
ax.loglog(tD,Gc3s1, color='orange') 
ax.loglog(tD,Gc3s2, color='orange')
ax.loglog(tD,Gc3s3, color='orange')
ax.loglog(tD,Gc3s4, color='orange')

ax.loglog(tD,Gc4s0, label =f'C = {cD[4]}', color='black')
ax.loglog(tD,Gc4s1, color='black') 
ax.loglog(tD,Gc4s2, color='black')
ax.loglog(tD,Gc4s3, color='black')
ax.loglog(tD,Gc4s4, color='black')

plt.xlabel(r'$\mathbf{Tempo}$ ($\mathbf{log(t_D)}$)')
plt.ylabel(r'$\mathbf{Pressão}$ ($\mathbf{log(p_D)}$)')
plt.legend(title = "C value:")
plt.grid(True, which="both", ls="-.")
plt.axis([1e2,1e8,1e-1,1e2])
plt.show()
