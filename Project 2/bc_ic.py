"""
This script contains the functions for boundary conditions and initial conditions.
"""

import parameters as pm
import numpy as np

def bce():
    """
    This function defines the no normal flow on the east boundary for u. Values can be adjusted if BCs change.
    """
    return 0

def bcw():
    """
    This function defines the no normal flow on the west boundary for u. Values can be adjusted if BCs change.
    """
    return 0

def bcn():
    """
    This function defines the no normal flow on the north boundary for v. Values can be adjusted if BCs change.
    """
    return 0

def bcs():
    """
    This function defines the no normal flow on the south boundary for v. Values can be adjusted if BCs change.
    """
    return 0

def ic():
    """
    This function defines the initial conditions and returns the arrays for u, v and n on various grids.
    """
    #u and v on u grid
    u_on_u = np.zeros((pm.pm['nx']-1, pm.pm['nx']))
    v_on_u = np.zeros((pm.pm['nx']-1, pm.pm['nx']))

    #u and v on v grid
    v_on_v = np.zeros((pm.pm['nx'], pm.pm['nx']-1))
    u_on_v = np.zeros((pm.pm['nx'], pm.pm['nx']-1))

    #u, v and n on n grid
    n_on_n = np.zeros((pm.pm['nx']-1, pm.pm['nx']-1))
    u_on_n = np.zeros((pm.pm['nx']-1, pm.pm['nx']-1))
    v_on_n = np.zeros((pm.pm['nx']-1, pm.pm['nx']-1))
    return u_on_u, v_on_v, n_on_n, u_on_n, v_on_n, v_on_u, u_on_v

def bcsl(j_min, j_max, i_min, i_max, array = None, u = False, v = False, n = False, \
         f_var = False, func = None, f_array = None):
    """
    First variant: This function defines the boundary conditions for the non-linear model. Departure points outside of grid are
    handled here, returning a, b, c, d which correspond to the 4 corners of the grid based on the
    departure point. a, b, c, d can be used for interpolation.
    Second variant: This function allows an input function (such as for calculating gradient) to be input to return
    a, b, c, d. This is used when array values need to be calculated (eg. du_dx on n grid).
    """
    #determine if grid boundaries are for u, v or n grid
    if n == True:
        jlim = 2
        ilim = 2
    elif u == True:
        jlim = 2
        ilim = 1
    elif v == True:
        jlim = 1
        ilim = 2
    #first variant, when grid corner values can be taken directly from input array values
    if f_var == False:
        #south boundary
        if j_min < 0:
            if i_min < 0:
                a, b, c, d = array[j_max,i_max], array[j_max,i_max], array[j_max,i_max], array[j_max,i_max]
            elif i_max > pm.pm['nx']-ilim:
                a, b, c, d = array[j_max,i_min], array[j_max,i_min], array[j_max,i_min], array[j_max,i_min]
            else:
                a, c = array[j_max,i_min], array[j_max,i_min]
                b, d = array[j_max,i_max], array[j_max,i_max]
        #west boundary
        elif i_min < 0:
            if j_max > pm.pm['nx']-jlim:
                a, b, c, d = array[j_min,i_max], array[j_min,i_max], array[j_min,i_max], array[j_min,i_max]
            else:
                a, b = array[j_min,i_max], array[j_min,i_max]
                c, d = array[j_max,i_max], array[j_max,i_max]
        #north boundary
        elif j_max > pm.pm['nx']-jlim:
            if i_max > pm.pm['nx']-ilim:
                a, b, c, d = array[j_min,i_min], array[j_min,i_min], array[j_min,i_min], array[j_min,i_min]
            else:
                a, c = array[j_min,i_min], array[j_min,i_min]
                b, d = array[j_min,i_max], array[j_min,i_max]
        #east boundary
        elif i_max > pm.pm['nx']-ilim:
            a, b = array[j_min,i_min], array[j_min,i_min]
            c, d = array[j_max,i_min], array[j_max,i_min]
        else:
            a, b, c, d = array[j_min,i_min], array[j_min,i_max], array[j_max,i_min], array[j_max,i_max]
        return a, b, c, d
    #second variant, when grid corner values need to be calculated from input array values
    if f_var == True:
        f = func
        #south boundary
        if j_min < 0:
            if i_min < 0:
                a, b, c, d = f(j_max,i_max,f_array), f(j_max,i_max,f_array), f(j_max,i_max,f_array), f(j_max,i_max,f_array)
            elif i_max > pm.pm['nx']-ilim:
                a, b, c, d = f(j_max,i_min,f_array), f(j_max,i_min,f_array), f(j_max,i_min,f_array), f(j_max,i_min,f_array)
            else:
                a, c = f(j_max,i_min,f_array), f(j_max,i_min,f_array)
                b, d = f(j_max,i_max,f_array), f(j_max,i_max,f_array)
        #west boundary
        elif i_min < 0:
            if j_max > pm.pm['nx']-jlim:
                a, b, c, d = f(j_min,i_max,f_array), f(j_min,i_max,f_array), f(j_min,i_max,f_array), f(j_min,i_max,f_array)
            else:
                a, b = f(j_min,i_max,f_array), f(j_min,i_max,f_array)
                c, d = f(j_max,i_max,f_array), f(j_max,i_max,f_array)
        #north boundary
        elif j_max > pm.pm['nx']-jlim:
            if i_max > pm.pm['nx']-ilim:
                a, b, c, d = f(j_min,i_min,f_array), f(j_min,i_min,f_array), f(j_min,i_min,f_array), f(j_min,i_min,f_array)
            else:
                a, c = f(j_min,i_min,f_array), f(j_min,i_min,f_array)
                b, d = f(j_min,i_max,f_array), f(j_min,i_max,f_array)
        #east boundary
        elif i_max > pm.pm['nx']-ilim:
            a, b = f(j_min,i_min,f_array), f(j_min,i_min,f_array)
            c, d = f(j_max,i_min,f_array), f(j_max,i_min,f_array)
        else:
            a, b, c, d = f(j_min,i_min,f_array), f(j_min,i_max,f_array), f(j_max,i_min,f_array), f(j_max,i_max,f_array)
        return a, b, c, d
