"""
This script contains all the utility functions required for the model.
"""
import parameters as pm
import bc_ic as bc
import math

def du_dx_on_n(j,i,u_on_u):
    """
    This function calculates du_dx on the n grid given u grid values.
    """
    du_dx = (u_on_u[j,i+1] - u_on_u[j,i])/float(pm.pm['d'])
    return du_dx

def dv_dy_on_n(j,i,v_on_v):
    """
    This function calculates dv_dy on the n grid given v grid values.
    """
    dv_dy = (v_on_v[j+1,i] - v_on_v[j,i])/float(pm.pm['d'])
    return dv_dy

def dn_dx_on_u(j,i,n_on_n):
    """
    This function calculates dn_dx on the u grid given n grid values.
    """
    if (i-1 < 0):
        a = n_on_n[j,i]
        b = n_on_n[j,i]
    elif (i > pm.pm['nx']-2):
        a = n_on_n[j,i-1]
        b = n_on_n[j,i-1]
    else:
        a = n_on_n[j,i]
        b = n_on_n[j,i-1]
    return (a - b)/float(pm.pm['d'])

def dn_dy_on_v(j,i,n_on_n):
    """
    This function calculates dn_dy on the v grid given n grid values.
    """
    if (j-1 < 0):
        a = n_on_n[j,i]
        b = n_on_n[j,i]
    elif (j > pm.pm['nx']-2):
        a = n_on_n[j-1,i]
        b = n_on_n[j-1,i]
    else:
        a = n_on_n[j,i]
        b = n_on_n[j-1,i]
    return (a - b)/float(pm.pm['d'])

def v_on_u(j,i,v_on_v):
    """
    This function calculates v on the u grid given v grid values.
    """
    if (i-1 < 0):
        a = v_on_v[j,i]
        b = v_on_v[j+1,i]
        c = bc.bcw()
        d = bc.bcw()
    elif (i > pm.pm['nx']-2):
        a = bc.bce()
        b = bc.bce()
        c = v_on_v[j,i-1]
        d = v_on_v[j+1,i-1]
    else:
        a = v_on_v[j,i]
        b = v_on_v[j+1,i]
        c = v_on_v[j,i-1]
        d = v_on_v[j+1,i-1]
    return (a + b + c + d)/4.

def u_on_v(j,i,u_on_u):
    """
    This function calculates u on the v grid given u grid values.
    """
    if (j-1 < 0):
        a = u_on_u[j,i]
        b = u_on_u[j,i+1]
        c = bc.bcn()
        d = bc.bcn()
    elif (j > pm.pm['nx']-2):
        a = bc.bcs()
        b = bc.bcs()
        c = u_on_u[j-1,i]
        d = u_on_u[j-1,i+1]
    else:
        a = u_on_u[j,i]
        b = u_on_u[j,i+1]
        c = u_on_u[j-1,i]
        d = u_on_u[j-1,i+1]
    return (a + b + c + d)/4.

def u_on_n(j,i,u_on_u):
    """
    This function calculates u on the n grid given u grid values.
    """
    u_centred = (u_on_u[j,i+1] + u_on_u[j,i])/2.
    return u_centred

def v_on_n(j,i,v_on_v):
    """
    This function calculates v on the n grid given v grid values.
    """
    v_centred = (v_on_v[j+1,i] + v_on_v[j,i])/2.
    return v_centred

def f1(x):
    """
    This function calculates f1 required for the analytical solution.
    """
    eps = pm.pm['gamma']/float(pm.pm['L']*pm.pm['beta'])
    a = (-1 - math.sqrt(1+(2*math.pi*eps)**2))/float(2*eps)
    b = (-1 + math.sqrt(1+(2*math.pi*eps)**2))/float(2*eps)
    num = (math.exp(a)-1)*math.exp(b*x) + (1-math.exp(b))*math.exp(a*x)
    den = math.exp(b) - math.exp(a)
    return math.pi*(1 + num/float(den))

def f2(x):
    """
    This function calculates f2 required for the analytical solution.
    """
    eps = pm.pm['gamma']/float(pm.pm['L']*pm.pm['beta'])
    a = (-1 - math.sqrt(1+(2*math.pi*eps)**2))/float(2*eps)
    b = (-1 + math.sqrt(1+(2*math.pi*eps)**2))/float(2*eps)
    num = (math.exp(a)-1)*b*math.exp(b*x) + (1-math.exp(b))*a*math.exp(a*x)
    den = math.exp(b) - math.exp(a)
    return num/float(den)

#work in n grid index space
def find_d_point_on_n(j,i,v_on_n_old,u_on_n_old,u_on_n,v_on_n):
    """
    This function calculates the departure point given velocity fields on n grid at 2 previous time levels.
    The input is a grid point (j,i) and returns (xd_j, xd_i) with respect to the n grid.
    """
    #find midpoint
    xstar_j = j - v_on_n[j,i]*pm.pm['dt']/float(2*pm.pm['d'])
    xstar_i = i - u_on_n[j,i]*pm.pm['dt']/float(2*pm.pm['d'])
    #find u at midpoint using previous time levels
    u_at_half_dt = 1.5*bi_inter_arr(xstar_j,xstar_i,u_on_n,var='n') - 0.5*bi_inter_arr(xstar_j,xstar_i,u_on_n_old,var='n')
    v_at_half_dt = 1.5*bi_inter_arr(xstar_j,xstar_i,v_on_n,var='n') - 0.5*bi_inter_arr(xstar_j,xstar_i,v_on_n_old,var='n')
    #find departure point (in index space)
    xd_j = j - v_at_half_dt*pm.pm['dt']/float(pm.pm['d'])
    xd_i = i - u_at_half_dt*pm.pm['dt']/float(pm.pm['d'])
    return xd_j, xd_i

def find_four_corners(y_index,x_index,n=False,u=False,v=False):
    """
    This function finds the four corners on the grid (of n, u or v) a given a departure point (y_index, x_index)
    with respect to the ocean basin.
    This is used for interpolation.
    """
    #fix basin grid index to corresponding u, v or n grid index
    if n == True:
        jd = y_index - 0.5
        id = x_index - 0.5
    if u == True:
        jd = y_index - 0.5
        id = x_index
    if v == True:
        jd = y_index
        id = x_index - 0.5
    i_min = int(math.floor(id))
    i_max = int(math.ceil(id))
    j_min = int(math.floor(jd))
    j_max = int(math.ceil(jd))
    return j_min, j_max, i_min, i_max, jd, id

def interpolate(a, b, c, d, j_min, j_max, i_min, i_max, y_index, x_index):
    """
    This function interpolates given values a, b, c, d corresponding to values at the 4 corners of the grid of
    a departure point.
    """
    #interpolate on y-direction
    if j_min == j_max:
        y_alfa = 0.5
        y_beta = 0.5
    else:
        y_alfa = (y_index - j_min)/1.
        y_beta = (j_max - y_index)/1.
    lb = y_beta*a + y_alfa*c
    rb = y_beta*b + y_alfa*d
    #interpolate on x-direction
    if i_min == i_max:
        x_alfa = 0.5
        x_beta = 0.5
    else:
        x_alfa = (x_index - i_min)/1.
        x_beta = (i_max - x_index)/1.
    return x_beta*lb + x_alfa*rb

def bi_inter_arr(y_index,x_index,array,var = None):
    """
    This function interpolates for u, v and n given the respective arrays. Takes values from the array and
    interpolates using the interpolate function.
    """
    i_min = int(math.floor(x_index))
    i_max = int(math.ceil(x_index))
    j_min = int(math.floor(y_index))
    j_max = int(math.ceil(y_index))
    if var == 'n':
        a, b, c, d = bc.bcsl(j_min, j_max, i_min, i_max, array=array, n=True)
    elif var == 'u':
        a, b, c, d = bc.bcsl(j_min, j_max, i_min, i_max, array=array, u=True)
    elif var == 'v':
        a, b, c, d = bc.bcsl(j_min, j_max, i_min, i_max, array=array, v=True)
    return interpolate(a, b, c, d, j_min, j_max, i_min, i_max, y_index, x_index)

def assign_corners(corners):
    """
    This function assigns the tuple from find_four_corners to the relevant indexes. This is used for
    performing functions at each corner grid point.
    """
    i_min = corners[2]
    i_max = corners[3]
    j_min = corners[0]
    j_max = corners[1]
    y_index = corners[4]
    x_index = corners[5]
    return i_min, i_max, j_min, j_max, y_index, x_index

def bi_inter_du_dx_on_narr(corners, array, function):
    """
    This function calculates du_dx on the n grid for the 4 corners of the grid of a departure point and interpolates
    to find the value for du_dx at the departure point.
    """
    i_min, i_max, j_min, j_max, y_index, x_index = assign_corners(corners)
    a, b, c, d = bc.bcsl(j_min, j_max, i_min, i_max, n=True, f_var=True, func=function, f_array=array)
    return interpolate(a, b, c, d, j_min, j_max, i_min, i_max, y_index, x_index)

def bi_inter_dv_dy_on_narr(corners, array, function):
    """
    This function calculates dv_dy on the n grid for the 4 corners of the grid of a departure point and interpolates
    to find the value for dv_dy at the departure point.
    """
    i_min, i_max, j_min, j_max, y_index, x_index = assign_corners(corners)
    a, b, c, d = bc.bcsl(j_min, j_max, i_min, i_max, n=True, f_var=True, func=function, f_array=array)
    return interpolate(a, b, c, d, j_min, j_max, i_min, i_max, y_index, x_index)

def bi_inter_v_on_uarr(corners, array, function):
    """
    This function calculates v on the u grid for the 4 corners of the grid of a departure point and interpolates
    to find the value for v at the departure point.
    """
    i_min, i_max, j_min, j_max, y_index, x_index = assign_corners(corners)
    a, b, c, d = bc.bcsl(j_min, j_max, i_min, i_max, u=True, f_var=True, func=function, f_array=array)
    return interpolate(a, b, c, d, j_min, j_max, i_min, i_max, y_index, x_index)

def bi_inter_u_on_varr(corners, array, function):
    """
    This function calculates u on the v grid for the 4 corners of the grid of a departure point and interpolates
    to find the value for u at the departure point.
    """
    i_min, i_max, j_min, j_max, y_index, x_index = assign_corners(corners)
    a, b, c, d = bc.bcsl(j_min, j_max, i_min, i_max, v=True, f_var=True, func=function, f_array=array)
    return interpolate(a, b, c, d, j_min, j_max, i_min, i_max, y_index, x_index)

def bi_inter_dn_dx_on_uarr(corners, array, function):
    """
    This function calculates dn_dx on the u grid for the 4 corners of the grid of a departure point and interpolates
    to find the value for dn_dx at the departure point.
    """
    i_min, i_max, j_min, j_max, y_index, x_index = assign_corners(corners)
    a, b, c, d = bc.bcsl(j_min, j_max, i_min, i_max, u=True, f_var=True, func=function, f_array=array)
    return interpolate(a, b, c, d, j_min, j_max, i_min, i_max, y_index, x_index)

def bi_inter_dn_dy_on_varr(corners, array, function):
    """
    This function calculates dn_dy on the v grid for the 4 corners of the grid of a departure point and interpolates
    to find the value for dn_dy at the departure point.
    """
    i_min, i_max, j_min, j_max, y_index, x_index = assign_corners(corners)
    a, b, c, d = bc.bcsl(j_min, j_max, i_min, i_max, v=True, f_var=True, func=function, f_array=array)
    return interpolate(a, b, c, d, j_min, j_max, i_min, i_max, y_index, x_index)