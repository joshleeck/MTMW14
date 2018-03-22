"""
This script contains the schemes (linear and non-linear SWE model) which are called in main.py,
also contains the analytic solution. Functions for interpolating for u on v and n grid, v on u and n grid and
calculating energy is only used for plotting (and does not affect the model run), so it is
defined in this script.
"""

import parameters as pm
import functions as f
import math

def fbts_odd(n_on_n, u_on_u, v_on_v):
    """
    This function runs the forward-backward time scheme for one odd timestep for the linear SWE model.
    The order of solving the momentum equations are switched for even timesteps (below).
    """
    #n equation
    for j in range(len(n_on_n)):
        for i in range(len(n_on_n[0])):
            du_dx = f.du_dx_on_n(j,i,u_on_u)
            dv_dy = f.dv_dy_on_n(j,i,v_on_v)
            n_on_n[j,i] = n_on_n[j,i] - pm.pm['H']*pm.pm['dt']*(du_dx + dv_dy)
    #u equation with no normal flow BCs
    for j in range(len(u_on_u)):
        for i in range(1,len(u_on_u[0])-1):
            y = (0.5+j)*pm.pm['d']
            dn_dx = f.dn_dx_on_u(j,i,n_on_n)
            u_on_u[j,i] = u_on_u[j,i] + (pm.pm['f0'] + pm.pm['beta']*y)*pm.pm['dt']*f.v_on_u(j,i,v_on_v)          \
                          - pm.pm['g']*pm.pm['dt']*dn_dx - pm.pm['gamma']*pm.pm['dt']*u_on_u[j,i]                 \
                          + pm.tau(pm.pm['tau0'], y, pm.pm['L'])[0]*pm.pm['dt']/(float(pm.pm['roe']*pm.pm['H']))
    #v equation with no normal flow BCs
    for j in range(1,len(v_on_v)-1):
        for i in range(len(v_on_v[0])):
            y = j*pm.pm['d']
            dn_dy = f.dn_dy_on_v(j,i,n_on_n)
            v_on_v[j,i] = v_on_v[j,i] - (pm.pm['f0'] + pm.pm['beta']*y)*pm.pm['dt']*f.u_on_v(j,i,u_on_u)          \
                          - pm.pm['g']*pm.pm['dt']*dn_dy - pm.pm['gamma']*pm.pm['dt']*v_on_v[j,i]                 \
                          + pm.tau(pm.pm['tau0'], y, pm.pm['L'])[1]*pm.pm['dt']/(float(pm.pm['roe']*pm.pm['H']))

def fbts_even(n_on_n, u_on_u, v_on_v):
    """
    This function runs the forward-backward time scheme for one even timestep for the linear SWE model.
    The order of solving the momentum equations are switched for odd timesteps (above).
    """
    #n equation
    for j in range(len(n_on_n)):
        for i in range(len(n_on_n[0])):
            du_dx = f.du_dx_on_n(j, i, u_on_u)
            dv_dy = f.dv_dy_on_n(j,i,v_on_v)
            n_on_n[j,i] = n_on_n[j,i] - pm.pm['H']*pm.pm['dt']*(du_dx + dv_dy)
    #v equation with no normal flow BCs
    for j in range(1,len(v_on_v)-1):
        for i in range(len(v_on_v[0])):
            y = j*pm.pm['d']
            dn_dy = f.dn_dy_on_v(j,i,n_on_n)
            v_on_v[j,i] = v_on_v[j,i] - (pm.pm['f0'] + pm.pm['beta']*y)*pm.pm['dt']*f.u_on_v(j,i,u_on_u)          \
                          - pm.pm['g']*pm.pm['dt']*dn_dy - pm.pm['gamma']*pm.pm['dt']*v_on_v[j,i]                 \
                          + pm.tau(pm.pm['tau0'], y, pm.pm['L'])[1]*pm.pm['dt']/(float(pm.pm['roe']*pm.pm['H']))
    #u equation with no normal flow BCs
    for j in range(len(u_on_u)):
        for i in range(1,len(u_on_u[0])-1):
            y = (0.5+j)*pm.pm['d']
            dn_dx = f.dn_dx_on_u(j,i,n_on_n)
            u_on_u[j,i] = u_on_u[j,i] + (pm.pm['f0'] + pm.pm['beta']*y)*pm.pm['dt']*f.v_on_u(j,i,v_on_v)          \
                          - pm.pm['g']*pm.pm['dt']*dn_dx - pm.pm['gamma']*pm.pm['dt']*u_on_u[j,i]                 \
                          + pm.tau(pm.pm['tau0'], y, pm.pm['L'])[0]*pm.pm['dt']/(float(pm.pm['roe']*pm.pm['H']))

def inter_u_on_n(n_on_n,u_on_u,u_on_n):
    """
    This function interpolates for u on n grid using u on u grid values.
    """
    for j in range(len(n_on_n)):
        for i in range(len(n_on_n[0])):
            u_on_n[j,i] = f.u_on_n(j,i,u_on_u)
    return u_on_n

def inter_v_on_n(n_on_n,v_on_v,v_on_n):
    """
    This function interpolates for v on n grid using v on v grid values.
    """
    for j in range(len(n_on_n)):
        for i in range(len(n_on_n[0])):
            v_on_n[j,i] = f.v_on_n(j,i,v_on_v)
    return v_on_n

def inter_v_on_u(u_on_u, v_on_v, v_on_u):
    """
    This function interpolates for v on u grid using v on v grid values.
    """
    for j in range(len(u_on_u)):
        for i in range(len(u_on_u[0])):
            v_on_u[j,i] = f.v_on_u(j,i,v_on_v)
    return v_on_u

def inter_u_on_v(u_on_u, v_on_v, u_on_v):
    """
    This function interpolates for u on v grid using u on u grid values.
    """
    for j in range(len(v_on_v)):
        for i in range(len(v_on_v[0])):
            u_on_v[j,i] = f.u_on_v(j,i,u_on_u)
    return u_on_v

def calc_E(n_on_n, u_on_n, v_on_n):
    """
    This function calculates the the total energy perturbation from the resting ocean.
    The equation used is shown in the report.
    """
    #integrate numerically starting from 0
    total = 0
    for j in range(len(n_on_n)):
        for i in range(len(n_on_n[0])):
            a = pm.pm['H']*((u_on_n[j,i]**2) + (v_on_n[j,i]**2))
            b = pm.pm['g']*(n_on_n[j,i]**2)
            total += (pm.pm['d']**2)*0.5*pm.pm['roe']*(a + b)
    return total

def analytic_solution(n_on_n, u_on_n, v_on_n, n0):
    """
    This function calculates the steady state analytical solution. The equations are in Mushgrave (1985).
    """
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

def SL_fbts_odd(n_on_n_new, v_on_v_new, u_on_u_new, n_on_n, u_on_u, v_on_v,
                u_on_n, v_on_n, v_on_n_old, u_on_n_old):
    """
    This function runs the forward-backward time scheme for one odd timestep for the non-linear SWE model.
    The order of solving the momentum equations are switched for even timesteps (below).
    """
    #n equation
    for j in range(len(n_on_n)):
        for i in range(len(n_on_n[0])):
            #remember that 0,0 on n-grid is 0.5,0.5 on basin grid for departure point
            #remember that 0,0 on u-grid is 0.5,0 on basin grid for departure point
            #remember that 0,0 on v-grid is 0,0.5 on basin grid for departure point
            jd, id = f.find_d_point_on_n(j,i,v_on_n_old,u_on_n_old,u_on_n,v_on_n)
            jd_basin = jd + 0.5
            id_basin = id + 0.5
            #find the 4 grid corners of the departure point on the n grid
            corners_n = f.find_four_corners(jd_basin, id_basin, n = True)
            du_dx = f.du_dx_on_n(j,i,u_on_u)
            dv_dy = f.dv_dy_on_n(j,i,v_on_v)
            du_dx_at_dp = f.bi_inter_du_dx_on_narr(corners_n,u_on_u,f.du_dx_on_n)
            dv_dy_at_dp = f.bi_inter_dv_dy_on_narr(corners_n,v_on_v,f.dv_dy_on_n)
            n_on_n_new[j,i] = f.bi_inter_arr(jd,id,n_on_n,var='n') - pm.pm['dt']*pm.pm['H']*(du_dx + dv_dy
                              + du_dx_at_dp + dv_dy_at_dp)/2.
    #u equation with no normal flow BCs
    for j in range(len(u_on_u)):
        for i in range(1,len(u_on_u[0])-1):
            jd, id = f.find_d_point_on_n(j,i,v_on_n_old,u_on_n_old,u_on_n,v_on_n)
            jd_basin = jd + 0.5
            id_basin = id + 0.5
            #find the 4 grid corners of the departure point on the u grid
            corners_u = f.find_four_corners(jd_basin, id_basin, u = True)
            y = (0.5+j)*pm.pm['d']
            y_dp = jd_basin*pm.pm['d']
            dn_dx_new = f.dn_dx_on_u(j,i,n_on_n_new)
            #source terms at x_ij
            S_ij = (pm.pm['f0'] + pm.pm['beta']*y)*f.v_on_u(j,i,v_on_v)  \
                   - pm.pm['g']*dn_dx_new      \
                   + pm.tau(pm.pm['tau0'], y, pm.pm['L'])[0]/(float(pm.pm['roe']*pm.pm['H']))
            #source terms at departure point
            S_dp = (pm.pm['f0'] + pm.pm['beta']*y_dp)*f.bi_inter_v_on_uarr(corners_u,v_on_v,f.v_on_u) \
                   - pm.pm['g']*f.bi_inter_dn_dx_on_uarr(corners_u,n_on_n,f.dn_dx_on_u) \
                   - pm.pm['gamma']*f.bi_inter_arr(jd,id+0.5,u_on_u,var='u') \
                   + pm.tau(pm.pm['tau0'], y_dp, pm.pm['L'])[0]/(float(pm.pm['roe']*pm.pm['H']))
            u_on_u_new[j,i] = (f.bi_inter_arr(jd,id+0.5,u_on_u,var='u') + pm.pm['dt']*(S_ij + S_dp)/2.)/    \
                              ((pm.pm['gamma']*pm.pm['dt'])/2. + 1)
    #v equation with no normal flow BCs
    for j in range(1,len(v_on_v)-1):
        for i in range(len(v_on_v[0])):
            jd, id = f.find_d_point_on_n(j,i,v_on_n_old,u_on_n_old,u_on_n,v_on_n)
            jd_basin = jd + 0.5
            id_basin = id + 0.5
            #find the 4 grid corners of the departure point on the v grid
            corners_v = f.find_four_corners(jd_basin, id_basin, v = True)
            y = j*pm.pm['d']
            y_dp = jd_basin*pm.pm['d']
            dn_dy_new = f.dn_dy_on_v(j,i,n_on_n_new)
            #source terms at x_ij
            S_ij = -(pm.pm['f0'] + pm.pm['beta']*y)*f.u_on_v(j,i,u_on_u_new)  \
                   - pm.pm['g']*dn_dy_new      \
                   + pm.tau(pm.pm['tau0'], y, pm.pm['L'])[1]/(float(pm.pm['roe']*pm.pm['H']))
            #source terms at departure point
            S_dp = -(pm.pm['f0'] + pm.pm['beta']*y_dp)*f.bi_inter_u_on_varr(corners_v,u_on_u,f.u_on_v) \
                   - pm.pm['g']*f.bi_inter_dn_dy_on_varr(corners_v,n_on_n,f.dn_dy_on_v) \
                   - pm.pm['gamma']*f.bi_inter_arr(jd+0.5,id,v_on_v,var='v') \
                   + pm.tau(pm.pm['tau0'], y_dp, pm.pm['L'])[1]/(float(pm.pm['roe']*pm.pm['H']))
            v_on_v_new[j,i] = (f.bi_inter_arr(jd+0.5,id,v_on_v,var='v') + pm.pm['dt']*(S_ij + S_dp)/2.)/ \
                              ((pm.pm['gamma']*pm.pm['dt'])/2. + 1)


def SL_fbts_even(n_on_n_new, v_on_v_new, u_on_u_new, n_on_n, u_on_u, v_on_v,
                 u_on_n, v_on_n, v_on_n_old, u_on_n_old):
    #n equation
    for j in range(len(n_on_n)):
        for i in range(len(n_on_n[0])):
            #remember that 0,0 on n-grid is 0.5,0.5 on basin grid for departure point
            #remember that 0,0 on u-grid is 0.5,0 on basin grid for departure point
            #remember that 0,0 on v-grid is 0,0.5 on basin grid for departure point
            jd, id = f.find_d_point_on_n(j,i,v_on_n_old,u_on_n_old,u_on_n,v_on_n)
            jd_basin = jd + 0.5
            id_basin = id + 0.5
            #find the 4 grid corners of the departure point on the n grid
            corners_n = f.find_four_corners(jd_basin, id_basin, n = True)
            du_dx = f.du_dx_on_n(j,i,u_on_u)
            dv_dy = f.dv_dy_on_n(j,i,v_on_v)
            du_dx_at_dp = f.bi_inter_du_dx_on_narr(corners_n,u_on_u,f.du_dx_on_n)
            dv_dy_at_dp = f.bi_inter_dv_dy_on_narr(corners_n,v_on_v,f.dv_dy_on_n)
            n_on_n_new[j,i] = f.bi_inter_arr(jd,id,n_on_n,var='n') - pm.pm['dt']*pm.pm['H']*(du_dx + dv_dy
                              + du_dx_at_dp + dv_dy_at_dp)/2.
    #v equation with no normal flow BCs
    for j in range(1,len(v_on_v)-1):
        for i in range(len(v_on_v[0])):
            jd, id = f.find_d_point_on_n(j,i,v_on_n_old,u_on_n_old,u_on_n,v_on_n)
            jd_basin = jd + 0.5
            id_basin = id + 0.5
            #find the 4 grid corners of the departure point on the v grid
            corners_v = f.find_four_corners(jd_basin, id_basin, v = True)
            y = j*pm.pm['d']
            y_dp = jd_basin*pm.pm['d']
            dn_dy_new = f.dn_dy_on_v(j,i,n_on_n_new)
            #source terms at x_ij
            S_ij = -(pm.pm['f0'] + pm.pm['beta']*y)*f.u_on_v(j,i,u_on_u)  \
                   - pm.pm['g']*dn_dy_new      \
                   + pm.tau(pm.pm['tau0'], y, pm.pm['L'])[1]/(float(pm.pm['roe']*pm.pm['H']))
            #source terms at departure point
            S_dp = -(pm.pm['f0'] + pm.pm['beta']*y_dp)*f.bi_inter_u_on_varr(corners_v,u_on_u,f.u_on_v) \
                   - pm.pm['g']*f.bi_inter_dn_dy_on_varr(corners_v,n_on_n,f.dn_dy_on_v) \
                   - pm.pm['gamma']*f.bi_inter_arr(jd+0.5,id,v_on_v,var='v') \
                   + pm.tau(pm.pm['tau0'], y_dp, pm.pm['L'])[1]/(float(pm.pm['roe']*pm.pm['H']))
            v_on_v_new[j,i] = (f.bi_inter_arr(jd+0.5,id,v_on_v,var='v') + pm.pm['dt']*(S_ij + S_dp)/2.)/ \
                              ((pm.pm['gamma']*pm.pm['dt'])/2. + 1)
    #u equation with no normal flow BCs
    for j in range(len(u_on_u)):
        for i in range(1,len(u_on_u[0])-1):
            jd, id = f.find_d_point_on_n(j,i,v_on_n_old,u_on_n_old,u_on_n,v_on_n)
            jd_basin = jd + 0.5
            id_basin = id + 0.5
            #find the 4 grid corners of the departure point on the u grid
            corners_u = f.find_four_corners(jd_basin, id_basin, u = True)
            y = (0.5+j)*pm.pm['d']
            y_dp = jd_basin*pm.pm['d']
            dn_dx_new = f.dn_dx_on_u(j,i,n_on_n_new)
            #source terms at x_ij
            S_ij = (pm.pm['f0'] + pm.pm['beta']*y)*f.v_on_u(j,i,v_on_v_new)  \
                   - pm.pm['g']*dn_dx_new      \
                   + pm.tau(pm.pm['tau0'], y, pm.pm['L'])[0]/(float(pm.pm['roe']*pm.pm['H']))
            #source terms at departure point
            S_dp = (pm.pm['f0'] + pm.pm['beta']*y_dp)*f.bi_inter_v_on_uarr(corners_u,v_on_v,f.v_on_u) \
                   - pm.pm['g']*f.bi_inter_dn_dx_on_uarr(corners_u,n_on_n,f.dn_dx_on_u) \
                   - pm.pm['gamma']*f.bi_inter_arr(jd,id+0.5,u_on_u,var='u') \
                   + pm.tau(pm.pm['tau0'], y_dp, pm.pm['L'])[0]/(float(pm.pm['roe']*pm.pm['H']))
            u_on_u_new[j,i] = (f.bi_inter_arr(jd,id+0.5,u_on_u,var='u') + pm.pm['dt']*(S_ij + S_dp)/2.)/    \
                              ((pm.pm['gamma']*pm.pm['dt'])/2. + 1)