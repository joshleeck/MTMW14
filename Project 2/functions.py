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
