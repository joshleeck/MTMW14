"""
This script describes all the parameters (pm) required for the model.
"""
import math
import numpy as np

#functions for parameters
def calc_Rd(g, H, f0):
    return math.sqrt(g*H)/float(f0)

def tau(tau0, y, L):
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
pm['nx'] = 51 #101 #26 #11
pm['d'] = pm['L']/float(pm['nx']-1)


#temporal values
pm['dt'] = 140   #70      #280    #280
pm['nt'] = 10#*40#618*80 #618*20  #618*20

#initial conditions (do a function that sets all values to 0?)