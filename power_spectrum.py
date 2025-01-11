import numpy as np
import matplotlib.pyplot as plt


A = 3.16e4
k0 = 0.05

def power_spectrum(k,k0=k0,A=A, ns = 0.5, nb=-2):
    if type(k) == float:
        if k < k0:
            return A*(k/k0)**ns
        else :
            return A*(k/k0)**nb
    else:
        k = np.array(k)
        return np.where(k<k0,A*(k/k0)**ns, A*(k/k0)**nb)


def density_fluctuation(k):

    return np.sqrt((2*np.pi)**3*power_spectrum(k))


