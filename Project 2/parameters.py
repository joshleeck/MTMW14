"""
This script describes all the parameters (pm) required for the model.
"""
import math
import numpy as np

#functions for parameters
def calc_Rd(g, H, f0):
    """
    This function calculates the Rossby deformation radius
    """
    return math.sqrt(g*H)/float(f0)

def tau(tau0, y, L):
    """
    This function calculates the wind forcing tau as a function of y.
    """
    if y < 0:
        y = 0
    if y > L:
        y = L
    tau = np.zeros(2)
    tau[0] = -tau0*math.cos(math.pi*y/float(L))
    return tau

#dictionary for storing constant parameters
pm = {}

#model parameters
#constants
pm['L'] = 1e6
pm['f0'] = 1e-4
pm['beta'] = 1e-11
pm['g'] = 10
pm['gamma'] = 1e-6
pm['roe'] = 1e3
pm['H'] = 1e3
pm['tau0'] = 0.2

#spatial values
pm['nx'] = 26 #26 #41 #51
pm['d'] = pm['L']/float(pm['nx']-1)

#temporal values
pm['dt'] = 150 #26 #100 #100
pm['nt'] = 10 #864*30 #864*30