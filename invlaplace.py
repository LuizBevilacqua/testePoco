import math
import numpy as np
from scipy import special as sp

class invlaplace:
    def __init__(self,t,func,l):
        self.t = t
        self.func = func
        self.l = l

    def gavsteh(self,funname, t, L):
        nn2 = L // 2
        nn21 = nn2 + 1
        v = np.zeros(L)
        
        for n in range(1, L+1):
            z = 0.0
            for k in range(math.floor((n + 1) / 2), min(n, nn2) + 1):
                z += ((k**nn2) * math.factorial(2*k)) / \
                    (math.factorial(nn2-k) * math.factorial(k) * math.factorial(k-1) * \
                    math.factorial(n-k) * math.factorial(2*k - n))
            v[n-1] = (-1)**(n+nn2) * z
        
        sum_val = 0.0
        ln2_on_t = math.log(2.0) / t
        
        for n in range(1, L+1):
            p = n * ln2_on_t
            sum_val += v[n-1] * funname(p)
        
        ilt = sum_val * ln2_on_t
        
        return ilt
    

