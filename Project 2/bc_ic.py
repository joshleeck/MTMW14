import parameters as pm
import math
import numpy as np

def bce():
    return 0

def bcw():
    return 0

def bcn():
    return 0

def bcs():
    return 0

def ic():
    u_on_u = np.zeros((pm.pm['nx']-1, pm.pm['nx']))
    v_on_u = np.zeros((pm.pm['nx']-1, pm.pm['nx']))

    v_on_v = np.zeros((pm.pm['nx'], pm.pm['nx']-1))
    u_on_v = np.zeros((pm.pm['nx'], pm.pm['nx']-1))

    n_on_n = np.zeros((pm.pm['nx']-1, pm.pm['nx']-1))
    u_on_n = np.zeros((pm.pm['nx']-1, pm.pm['nx']-1))
    v_on_n = np.zeros((pm.pm['nx']-1, pm.pm['nx']-1))
    return u_on_u, v_on_v, n_on_n, u_on_n, v_on_n, v_on_u, u_on_v
