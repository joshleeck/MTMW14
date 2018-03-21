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

def bcsl(j_min, j_max, i_min, i_max, array = None, u = False, v = False, n = False, \
         f_var = False, func = None, f_array = None):
    if n == True:
        jlim = 2
        ilim = 2
    elif u == True: #actually u and v supposed to have 0 at bc? but this works okay
        jlim = 2
        ilim = 1
    elif v == True:
        jlim = 1
        ilim = 2
    if f_var == False:
        if j_min < 0:
            if i_min < 0:
                a, b, c, d = array[j_max,i_max], array[j_max,i_max], array[j_max,i_max], array[j_max,i_max]
            elif i_max > pm.pm['nx']-ilim:
                a, b, c, d = array[j_max,i_min], array[j_max,i_min], array[j_max,i_min], array[j_max,i_min]
            else:
                a, c = array[j_max,i_min], array[j_max,i_min]
                b, d = array[j_max,i_max], array[j_max,i_max]
        elif i_min < 0:
            if j_max > pm.pm['nx']-jlim:
                a, b, c, d = array[j_min,i_max], array[j_min,i_max], array[j_min,i_max], array[j_min,i_max]
            else:
                a, b = array[j_min,i_max], array[j_min,i_max]
                c, d = array[j_max,i_max], array[j_min,i_max]
        elif j_max > pm.pm['nx']-jlim:
            if i_max > pm.pm['nx']-ilim:
                a, b, c, d = array[j_min,i_min], array[j_min,i_min], array[j_min,i_min], array[j_min,i_min]
            else:
                a, c = array[j_min,i_min], array[j_min,i_min]
                b, d = array[j_min,i_max], array[j_min,i_max]
        elif i_max > pm.pm['nx']-ilim:
            a, b = array[j_min,i_min], array[j_min,i_min]
            c, d = array[j_max,i_min], array[j_max,i_min]
        else:
            a, b, c, d = array[j_min,i_min], array[j_min,i_max], array[j_max,i_min], array[j_max,i_max]
        return a, b, c, d
    if f_var == True:
        f = func
        if j_min < 0:
            if i_min < 0:
                a, b, c, d = f(j_max,i_max,f_array), f(j_max,i_max,f_array), f(j_max,i_max,f_array), f(j_max,i_max,f_array)
            elif i_max > pm.pm['nx']-ilim:
                a, b, c, d = f(j_max,i_min,f_array), f(j_max,i_min,f_array), f(j_max,i_min,f_array), f(j_max,i_min,f_array)
            else:
                a, c = f(j_max,i_min,f_array), f(j_max,i_min,f_array)
                b, d = f(j_max,i_max,f_array), f(j_max,i_max,f_array)
        elif i_min < 0:
            if j_max > pm.pm['nx']-jlim:
                a, b, c, d = f(j_min,i_max,f_array), f(j_min,i_max,f_array), f(j_min,i_max,f_array), f(j_min,i_max,f_array)
            else:
                a, b = f(j_min,i_max,f_array), f(j_min,i_max,f_array)
                c, d = f(j_max,i_max,f_array), f(j_max,i_max,f_array)
        elif j_max > pm.pm['nx']-jlim:
            if i_max > pm.pm['nx']-ilim:
                a, b, c, d = f(j_min,i_min,f_array), f(j_min,i_min,f_array), f(j_min,i_min,f_array), f(j_min,i_min,f_array)
            else:
                a, c = f(j_min,i_min,f_array), f(j_min,i_min,f_array)
                b, d = f(j_min,i_max,f_array), f(j_min,i_max,f_array)
        elif i_max > pm.pm['nx']-ilim:
            a, b = f(j_min,i_min,f_array), f(j_min,i_min,f_array)
            c, d = f(j_max,i_min,f_array), f(j_max,i_min,f_array)
        else:
            a, b, c, d = f(j_min,i_min,f_array), f(j_min,i_max,f_array), f(j_max,i_min,f_array), f(j_max,i_max,f_array)
        return a, b, c, d
