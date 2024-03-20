def gavsteh_param(l, func, temp):
    """
    Calcula uma transformada inversa de Laplace usando o método de Gaver-Stehfest.

    :param l: O número de termos a serem usados na aproximação.
    :param func: Uma função que representa uma operação matemática que aceita um único parâmetro.
    :param temp: O valor no qual a função será avaliada.
    :return: O resultado da transformada inversa de Laplace.
    """
    import math
    import numpy as np

    # Inicializa variáveis
    t = temp
    n = int(l / 2)
    v = []

    # Calcula os coeficientes v[j] usando o método de Gaver-Stehfest
    for j in range(1, l+1):
        z = 0.0
        for k in range(((j+1)//2), min(j, n)+1):
            z = z + ((k**n) * math.factorial(2*k)) / (math.factorial(n-k) * math.factorial(k) * math.factorial(k-1) *
                                                      math.factorial(j-k) * math.factorial((2*k) - j))
        v.append((-1)**(j+n) * z)

    # Calcula a transformada inversa de Laplace usando os coeficientes
    somme = 0.
    ln2_on_t = np.log(2.) / t
    for j in range(1, l+1):
        p = j * ln2_on_t
        somme += (v[j-1] * func(p))

    ilt = somme * ln2_on_t
    return ilt
