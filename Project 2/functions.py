import parameters as pm
import bc_ic as bc
import math


def du_dx_on_n(j,i,u_on_u):
    du_dx = (u_on_u[j,i+1] - u_on_u[j,i])/float(pm.pm['d'])
    return du_dx

def dv_dy_on_n(j,i,v_on_v):
    dv_dy = (v_on_v[j+1,i] - v_on_v[j,i])/float(pm.pm['d'])
    return dv_dy

def dn_dx_on_u(j,i,n_on_n):
    if (i-1 < 0):
        a = n_on_n[j,i]
        b = bc.bcw()
    elif (i > pm.pm['nx']-2):
        a = bc.bce()
        b = n_on_n[j,i-1]
    else:
        a = n_on_n[j,i]
        b = n_on_n[j,i-1]
    return (a - b)/float(pm.pm['d'])

def dn_dy_on_v(j,i,n_on_n):
    if (j-1 < 0):
        a = n_on_n[j,i]
        b = bc.bcn()
    elif (j > pm.pm['nx']-2):
        a = bc.bcs()
        b = n_on_n[j-1,i]
    else:
        a = n_on_n[j,i]
        b = n_on_n[j-1,i]
    return (a - b)/float(pm.pm['d'])

def v_on_u(j,i,v_on_v):
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
    u_centred = (u_on_u[j,i+1] + u_on_u[j,i])/2.
    return u_centred

def v_on_n(j,i,v_on_v):
    v_centred = (v_on_v[j+1,i] + v_on_v[j,i])/2.
    return v_centred

def f1(x):
    eps = pm.pm['gamma']/float(pm.pm['L']*pm.pm['beta'])
    a = (-1 - math.sqrt(1+(2*math.pi*eps)**2))/float(2*eps)
    b = (-1 + math.sqrt(1+(2*math.pi*eps)**2))/float(2*eps)
    num = (math.exp(a)-1)*math.exp(b*x) + (1-math.exp(b))*math.exp(a*x)
    den = math.exp(b) - math.exp(a)
    return math.pi*(1 + num/float(den))

def f2(x):
    eps = pm.pm['gamma']/float(pm.pm['L']*pm.pm['beta'])
    a = (-1 - math.sqrt(1+(2*math.pi*eps)**2))/float(2*eps)
    b = (-1 + math.sqrt(1+(2*math.pi*eps)**2))/float(2*eps)
    num = (math.exp(a)-1)*b*math.exp(b*x) + (1-math.exp(b))*a*math.exp(a*x)
    den = math.exp(b) - math.exp(a)
    return num/float(den)

#########################################################
#work in n-index space
def find_d_point_on_n(j,i,v_on_n_old,u_on_n_old,u_on_n,v_on_n):
    #find midpoint
    xstar_j = j - v_on_n[j,i]*pm.pm['dt']/float(2*pm.pm['d'])
    xstar_i = i - u_on_n[j,i]*pm.pm['dt']/float(2*pm.pm['d'])
    #print xstar_j, xstar_i
    #find u at midpoint using previous time levels
    u_at_half_dt = 1.5*bi_inter_narr(xstar_j,xstar_i,u_on_n) - 0.5*bi_inter_narr(xstar_j,xstar_i,u_on_n_old)
    v_at_half_dt = 1.5*bi_inter_narr(xstar_j,xstar_i,v_on_n) - 0.5*bi_inter_narr(xstar_j,xstar_i,v_on_n_old)
    #find departure point (in index space)
    xd_j = j - v_at_half_dt*pm.pm['dt']/float(pm.pm['d'])
    xd_i = i - u_at_half_dt*pm.pm['dt']/float(pm.pm['d'])
    #print xd_j, xd_i
    return xd_j, xd_i

def find_four_corners(y_index,x_index,n=False,u=False,v=False):
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

def bi_inter_narr(y_index,x_index,array):
    i_min = int(math.floor(x_index))
    i_max = int(math.ceil(x_index))
    j_min = int(math.floor(y_index))
    j_max = int(math.ceil(y_index))
    try:
        a = array[j_min,i_min]
    except IndexError:
        a = bc.bce()
    try:
        b = array[j_min,i_max]
    except IndexError:
        b = bc.bcw()
    try:
        c = array[j_max,i_min]
    except IndexError:
        c = bc.bcs()
    try:
        d = array[j_max,i_max]
    except IndexError:
        d = bc.bcn()
    #interpolate on y-direction
    y_alfa = (y_index - j_min)/1.
    y_beta = (j_max - y_index)/1.
    lb = y_beta*a + y_alfa*c
    rb = y_beta*b + y_alfa*d
    #interpolate on x-direction
    x_alfa = (x_index - i_min)/1.
    x_beta = (i_max - x_index)/1.
    return x_beta*lb + x_alfa*rb

def bi_inter_uarr(y_index,x_index,array):
    i_min = int(math.floor(x_index))
    i_max = int(math.ceil(x_index))
    j_min = int(math.floor(y_index))
    j_max = int(math.ceil(y_index))
    try:
        a = array[j_min,i_min]
    except IndexError:
        a = bc.bce()
    try:
        b = array[j_min,i_max]
    except IndexError:
        b = bc.bcw()
    try:
        c = array[j_max,i_min]
    except IndexError:
        c = bc.bcs()
    try:
        d = array[j_max,i_max]
    except IndexError:
        d = bc.bcn()
    #interpolate on y-direction
    y_alfa = (y_index - j_min)/1.
    y_beta = (j_max - y_index)/1.
    lb = y_beta*a + y_alfa*c
    rb = y_beta*b + y_alfa*d
    #interpolate on x-direction
    x_alfa = (x_index - i_min)/1.
    x_beta = (i_max - x_index)/1.
    return x_beta*lb + x_alfa*rb

def bi_inter_varr(y_index,x_index,array):
    i_min = int(math.floor(x_index))
    i_max = int(math.ceil(x_index))
    j_min = int(math.floor(y_index))
    j_max = int(math.ceil(y_index))
    if j_min < 0:
        a = bc.bcs()
        b = bc.bcs()
        c = array[j_max,i_min]
        d = array[j_max,i_max]
    elif i_min < 0:
        a = bc.bcw()
        b = array[j_min,i_max]
        c = bc.bcw()
        d = array[j_max,i_max]
    elif j_max > pm.pm['nx']-1:
        a = array[j_min,i_min]
        b = array[j_min,i_max]
        c = bc.bcn()
        d = bc.bcn()
    elif i_max > pm.pm['nx']-2:
        a = array[j_min,i_min]
        b = bc.bce()
        c = array[j_max,i_min]
        d = bc.bce()
    else:
        a = array[j_min,i_min]
        b = array[j_min,i_max]
        c = array[j_max,i_min]
        d = array[j_max,i_max]
    #interpolate on y-direction
    y_alfa = (y_index - j_min)/1.
    y_beta = (j_max - y_index)/1.
    lb = y_beta*a + y_alfa*c
    rb = y_beta*b + y_alfa*d
    #interpolate on x-direction
    x_alfa = (x_index - i_min)/1.
    x_beta = (i_max - x_index)/1.
    return x_beta*lb + x_alfa*rb

def bi_inter_du_dx_on_narr(corners,array):
    if corners[0] < 0:
        a = bc.bcs()
        b = bc.bcs()
        c = du_dx_on_n(corners[1],corners[2],array)
        d = du_dx_on_n(corners[1],corners[3],array)
    elif corners[2] < 0:
        a = bc.bcw()
        b = du_dx_on_n(corners[0],corners[3],array)
        c = bc.bcw()
        d = du_dx_on_n(corners[1],corners[3],array)
    elif corners[1] > pm.pm['nx']-2:
        a = du_dx_on_n(corners[0], corners[2], array)
        b = du_dx_on_n(corners[0], corners[3], array)
        c = bc.bcn()
        d = bc.bcn()
    elif corners[3] > pm.pm['nx']-2:
        a = du_dx_on_n(corners[0], corners[2], array)
        b = bc.bce()
        c = du_dx_on_n(corners[1], corners[2], array)
        d = bc.bce()
    else:
        a = du_dx_on_n(corners[0],corners[2],array)
        b = du_dx_on_n(corners[0],corners[3],array)
        c = du_dx_on_n(corners[1],corners[2],array)
        d = du_dx_on_n(corners[1],corners[3],array)
    #interpolate on y-direction
    y_alfa = (corners[4] - corners[0])/1.
    y_beta = (corners[1] - corners[4])/1.
    lb = y_beta*a + y_alfa*c
    rb = y_beta*b + y_alfa*d
    #interpolate on x-direction
    x_alfa = (corners[5] - corners[2])/1.
    x_beta = (corners[3] - corners[5])/1.
    return x_beta*lb + x_alfa*rb

def bi_inter_dv_dy_on_narr(corners,array):
    if corners[0] < 0:
        a = bc.bcs()
        b = bc.bcs()
        c = dv_dy_on_n(corners[1],corners[2],array)
        d = dv_dy_on_n(corners[1],corners[3],array)
    elif corners[2] < 0:
        a = bc.bcw()
        b = dv_dy_on_n(corners[0],corners[3],array)
        c = bc.bcw()
        d = dv_dy_on_n(corners[1],corners[3],array)
    elif corners[1] > pm.pm['nx']-2:
        a = dv_dy_on_n(corners[0], corners[2], array)
        b = dv_dy_on_n(corners[0], corners[3], array)
        c = bc.bcn()
        d = bc.bcn()
    elif corners[3] > pm.pm['nx']-2:
        a = dv_dy_on_n(corners[0], corners[2], array)
        b = bc.bce()
        c = dv_dy_on_n(corners[1], corners[2], array)
        d = bc.bce()
    else:
        a = dv_dy_on_n(corners[0],corners[2],array)
        b = dv_dy_on_n(corners[0],corners[3],array)
        c = dv_dy_on_n(corners[1],corners[2],array)
        d = dv_dy_on_n(corners[1],corners[3],array)
    #interpolate on y-direction
    y_alfa = (corners[4] - corners[0])/1.
    y_beta = (corners[1] - corners[4])/1.
    lb = y_beta*a + y_alfa*c
    rb = y_beta*b + y_alfa*d
    #interpolate on x-direction
    x_alfa = (corners[5] - corners[2])/1.
    x_beta = (corners[3] - corners[5])/1.
    return x_beta*lb + x_alfa*rb

def bi_inter_v_on_uarr(corners,array):
    if corners[0] < 0:
        a = bc.bcs()
        b = bc.bcs()
        c = v_on_u(corners[1],corners[2],array)
        d = v_on_u(corners[1],corners[3],array)
    elif corners[2] < 0:
        a = bc.bcw()
        b = v_on_u(corners[0],corners[3],array)
        c = bc.bcw()
        d = v_on_u(corners[1],corners[3],array)
    elif corners[1] > pm.pm['nx']-2:
        a = v_on_u(corners[0], corners[2], array)
        b = v_on_u(corners[0], corners[3], array)
        c = bc.bcn()
        d = bc.bcn()
    elif corners[3] > pm.pm['nx']-1:
        a = v_on_u(corners[0], corners[2], array)
        b = bc.bce()
        c = v_on_u(corners[1], corners[2], array)
        d = bc.bce()
    else:
        a = v_on_u(corners[0],corners[2],array)
        b = v_on_u(corners[0],corners[3],array)
        c = v_on_u(corners[1],corners[2],array)
        d = v_on_u(corners[1],corners[3],array)
    #interpolate on y-direction
    y_alfa = (corners[4] - corners[0])/1.
    y_beta = (corners[1] - corners[4])/1.
    lb = y_beta*a + y_alfa*c
    rb = y_beta*b + y_alfa*d
    #interpolate on x-direction
    x_alfa = (corners[5] - corners[2])/1.
    x_beta = (corners[3] - corners[5])/1.
    return x_beta*lb + x_alfa*rb

def bi_inter_u_on_varr(corners,array):
    if corners[0] < 0:
        a = bc.bcs()
        b = bc.bcs()
        c = u_on_v(corners[1],corners[2],array)
        d = u_on_v(corners[1],corners[3],array)
    elif corners[2] < 0:
        a = bc.bcw()
        b = u_on_v(corners[0],corners[3],array)
        c = bc.bcw()
        d = u_on_v(corners[1],corners[3],array)
    elif corners[1] > pm.pm['nx']-1:
        a = u_on_v(corners[0], corners[2], array)
        b = u_on_v(corners[0], corners[3], array)
        c = bc.bcn()
        d = bc.bcn()
    elif corners[3] > pm.pm['nx']-2:
        a = u_on_v(corners[0], corners[2], array)
        b = bc.bce()
        c = u_on_v(corners[1], corners[2], array)
        d = bc.bce()
    else:
        a = u_on_v(corners[0],corners[2],array)
        b = u_on_v(corners[0],corners[3],array)
        c = u_on_v(corners[1],corners[2],array)
        d = u_on_v(corners[1],corners[3],array)
    #interpolate on y-direction
    y_alfa = (corners[4] - corners[0])/1.
    y_beta = (corners[1] - corners[4])/1.
    lb = y_beta*a + y_alfa*c
    rb = y_beta*b + y_alfa*d
    #interpolate on x-direction
    x_alfa = (corners[5] - corners[2])/1.
    x_beta = (corners[3] - corners[5])/1.
    return x_beta*lb + x_alfa*rb

def bi_inter_dn_dx_on_uarr(corners,array):
    if corners[0] < 0:
        a = bc.bcs()
        b = bc.bcs()
        c = dn_dx_on_u(corners[1],corners[2],array)
        d = dn_dx_on_u(corners[1],corners[3],array)
    elif corners[2] < 0:
        a = bc.bcw()
        b = dn_dx_on_u(corners[0],corners[3],array)
        c = bc.bcw()
        d = dn_dx_on_u(corners[1],corners[3],array)
    elif corners[1] > pm.pm['nx']-2:
        a = dn_dx_on_u(corners[0], corners[2], array)
        b = dn_dx_on_u(corners[0], corners[3], array)
        c = bc.bcn()
        d = bc.bcn()
    elif corners[3] > pm.pm['nx']-1:
        a = dn_dx_on_u(corners[0], corners[2], array)
        b = bc.bce()
        c = dn_dx_on_u(corners[1], corners[2], array)
        d = bc.bce()
    else:
        a = dn_dx_on_u(corners[0],corners[2],array)
        b = dn_dx_on_u(corners[0],corners[3],array)
        c = dn_dx_on_u(corners[1],corners[2],array)
        d = dn_dx_on_u(corners[1],corners[3],array)
    #interpolate on y-direction
    y_alfa = (corners[4] - corners[0])/1.
    y_beta = (corners[1] - corners[4])/1.
    lb = y_beta*a + y_alfa*c
    rb = y_beta*b + y_alfa*d
    #interpolate on x-direction
    x_alfa = (corners[5] - corners[2])/1.
    x_beta = (corners[3] - corners[5])/1.
    return x_beta*lb + x_alfa*rb

def bi_inter_dn_dy_on_varr(corners,array):
    if corners[0] < 0:
        a = bc.bcs()
        b = bc.bcs()
        c = dn_dy_on_v(corners[1],corners[2],array)
        d = dn_dy_on_v(corners[1],corners[3],array)
    elif corners[2] < 0:
        a = bc.bcw()
        b = dn_dy_on_v(corners[0],corners[3],array)
        c = bc.bcw()
        d = dn_dy_on_v(corners[1],corners[3],array)
    elif corners[1] > pm.pm['nx']-1:
        a = dn_dy_on_v(corners[0], corners[2], array)
        b = dn_dy_on_v(corners[0], corners[3], array)
        c = bc.bcn()
        d = bc.bcn()
    elif corners[3] > pm.pm['nx']-2:
        a = dn_dy_on_v(corners[0], corners[2], array)
        b = bc.bce()
        c = dn_dy_on_v(corners[1], corners[2], array)
        d = bc.bce()
    else:
        a = dn_dy_on_v(corners[0],corners[2],array)
        b = dn_dy_on_v(corners[0],corners[3],array)
        c = dn_dy_on_v(corners[1],corners[2],array)
        d = dn_dy_on_v(corners[1],corners[3],array)
    #interpolate on y-direction
    y_alfa = (corners[4] - corners[0])/1.
    y_beta = (corners[1] - corners[4])/1.
    lb = y_beta*a + y_alfa*c
    rb = y_beta*b + y_alfa*d
    #interpolate on x-direction
    x_alfa = (corners[5] - corners[2])/1.
    x_beta = (corners[3] - corners[5])/1.
    return x_beta*lb + x_alfa*rb