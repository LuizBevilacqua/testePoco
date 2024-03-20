import math
import numpy as np

class invlaplace:
    def __init__(self,t,func,l):
        self.t = t
        self.func = func
        self.l = l

    def gavsteh_param(self):

        n = int(self.l / 2)
        v = []

        for j in range(1, self.l+1):
            z = 0.0
            for k in range(((j+1)//2), min(j, n)+1):
                z = z + ((k**n) * math.factorial(2*k)) / (math.factorial(n-k) * math.factorial(k) * math.factorial(k-1) *
                                                        math.factorial(j-k) * math.factorial((2*k) - j))
            v.append((-1)**(j+n) * z)

        somme = 0.
        ln2_on_t = np.log(2.) / self.t
        for j in range(1, self.l+1):
            p = j * ln2_on_t
            somme += (v[j-1] * self.func(p))

        ilt = somme * ln2_on_t
        return ilt