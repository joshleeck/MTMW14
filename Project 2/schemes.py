import parameters as pm
import functions as f
import math


def fbts_odd(n_on_n, u_on_u, v_on_v):
    for j in range(len(n_on_n)):
        for i in range(len(n_on_n[0])):
            du_dx = f.du_dx_on_n(j,i,u_on_u)
            dv_dy = f.dv_dy_on_n(j,i,v_on_v)
            n_on_n[j,i] = n_on_n[j,i] - pm.pm['H']*pm.pm['dt']*(du_dx + dv_dy)
    for j in range(len(u_on_u)):
        for i in range(1,len(u_on_u[0])-1):
            y = (0.5+j)*pm.pm['d']
            dn_dx = f.dn_dx_on_u(j,i,n_on_n)
            u_on_u[j,i] = u_on_u[j,i] + (pm.pm['f0'] + pm.pm['beta']*y)*pm.pm['dt']*f.v_on_u(j,i,v_on_v)          \
                          - pm.pm['g']*pm.pm['dt']*dn_dx - pm.pm['gamma']*pm.pm['dt']*u_on_u[j,i]                 \
                          + pm.tau(pm.pm['tau0'], y, pm.pm['L'])[0]*pm.pm['dt']/(float(pm.pm['roe']*pm.pm['H']))
    for j in range(1,len(v_on_v)-1):
        for i in range(len(v_on_v[0])):
            y = j*pm.pm['d']
            dn_dy = f.dn_dy_on_v(j,i,n_on_n)
            v_on_v[j,i] = v_on_v[j,i] - (pm.pm['f0'] + pm.pm['beta']*y)*pm.pm['dt']*f.u_on_v(j,i,u_on_u)          \
                          - pm.pm['g']*pm.pm['dt']*dn_dy - pm.pm['gamma']*pm.pm['dt']*v_on_v[j,i]                 \
                          + pm.tau(pm.pm['tau0'], y, pm.pm['L'])[1]*pm.pm['dt']/(float(pm.pm['roe']*pm.pm['H']))

def fbts_even(n_on_n, u_on_u, v_on_v):
    for j in range(len(n_on_n)):
        for i in range(len(n_on_n[0])):
            du_dx = f.du_dx_on_n(j, i, u_on_u)
            dv_dy = f.dv_dy_on_n(j,i,v_on_v)
            n_on_n[j,i] = n_on_n[j,i] - pm.pm['H']*pm.pm['dt']*(du_dx + dv_dy)
    for j in range(1,len(v_on_v)-1):
        for i in range(len(v_on_v[0])):
            y = j*pm.pm['d']
            dn_dy = f.dn_dy_on_v(j,i,n_on_n)
            v_on_v[j,i] = v_on_v[j,i] - (pm.pm['f0'] + pm.pm['beta']*y)*pm.pm['dt']*f.u_on_v(j,i,u_on_u)          \
                          - pm.pm['g']*pm.pm['dt']*dn_dy - pm.pm['gamma']*pm.pm['dt']*v_on_v[j,i]                 \
                          + pm.tau(pm.pm['tau0'], y, pm.pm['L'])[1]*pm.pm['dt']/(float(pm.pm['roe']*pm.pm['H']))
    for j in range(len(u_on_u)):
        for i in range(1,len(u_on_u[0])-1):
            y = (0.5+j)*pm.pm['d']
            dn_dx = f.dn_dx_on_u(j,i,n_on_n)
            u_on_u[j,i] = u_on_u[j,i] + (pm.pm['f0'] + pm.pm['beta']*y)*pm.pm['dt']*f.v_on_u(j,i,v_on_v)          \
                          - pm.pm['g']*pm.pm['dt']*dn_dx - pm.pm['gamma']*pm.pm['dt']*u_on_u[j,i]                 \
                          + pm.tau(pm.pm['tau0'], y, pm.pm['L'])[0]*pm.pm['dt']/(float(pm.pm['roe']*pm.pm['H']))

def inter_u_on_n(n_on_n,u_on_u,u_on_n):
    for j in range(len(n_on_n)):
        for i in range(len(n_on_n[0])):
            u_on_n[j,i] = f.u_on_n(j,i,u_on_u)
    return u_on_n

def inter_v_on_n(n_on_n,v_on_v,v_on_n):
    for j in range(len(n_on_n)):
        for i in range(len(n_on_n[0])):
            v_on_n[j,i] = f.v_on_n(j,i,v_on_v)
    return v_on_n

def inter_v_on_u(u_on_u, v_on_v, v_on_u):
    for j in range(len(u_on_u)):
        for i in range(len(u_on_u[0])):
            v_on_u[j,i] = f.v_on_u(j,i,v_on_v)
    return v_on_u

def inter_u_on_v(u_on_u, v_on_v, u_on_v):
    for j in range(len(v_on_v)):
        for i in range(len(v_on_v[0])):
            u_on_v[j,i] = f.u_on_v(j,i,u_on_u)
    return u_on_v

def calc_E(n_on_n, u_on_n, v_on_n):
    total = 0
    for j in range(len(n_on_n)):
        for i in range(len(n_on_n[0])):
            a = pm.pm['H']*((u_on_n[j,i]**2) + (v_on_n[j,i]**2))
            b = pm.pm['g']*(n_on_n[j,i]**2)
            total += (pm.pm['d']**2)*0.5*pm.pm['roe']*(a + b)
    return total

def analytic_solution(n_on_n, u_on_n, v_on_n, n0):
    a = pm.pm['tau0']/float(math.pi*pm.pm['gamma']*pm.pm['roe']*pm.pm['H'])
    b = pm.pm['f0']*pm.pm['L']/float(pm.pm['g'])
    c = pm.pm['gamma']/float(pm.pm['f0']*math.pi)
    for j in range(len(u_on_n)):
        for i in range(len(u_on_n[0])):
            x = (0.5+i)*pm.pm['d']
            y = (0.5+j)*pm.pm['d']
            u_on_n[j,i] = -a*f.f1(x/float(pm.pm['L']))*math.cos(math.pi*y/float(pm.pm['L']))
    for j in range(len(v_on_n)):
        for i in range(len(v_on_n[0])):
            x = (0.5+i)*pm.pm['d']
            y = (0.5+j)*pm.pm['d']
            v_on_n[j,i] = a*f.f2(x/float(pm.pm['L']))*math.sin(math.pi*y/float(pm.pm['L']))
    for j in range(len(n_on_n)):
        for i in range(len(n_on_n[0])):
            x = (0.5+i)*pm.pm['d']
            y = (0.5+j)*pm.pm['d']
            d = math.sin(math.pi*y/float(pm.pm['L']))*(1+pm.pm['beta']*y/float(pm.pm['f0']))
            e = pm.pm['beta']*pm.pm['L']/float(pm.pm['f0']*math.pi)*math.cos(math.pi*y/float(pm.pm['L']))
            g = f.f2(x/float(pm.pm['L']))*math.cos(math.pi*y/float(pm.pm['L']))
            n_on_n[j,i] = n0 + a*b*(c*g + f.f1(x/float(pm.pm['L']))*(d+e)/math.pi)
