import numpy as np
import parameters as pm
import functions as f
import schemes as sch
import bc_ic as bc
import plotting as pt
import copy


#if task a = true do this, if task b = true do this etc. main(TaskA=True,TaskB=True...)
#initialise arrays with initial conditions
def main(TaskB = False, TaskC1 = False, TaskC2 = False, TaskD = False, TaskD1 = False):
    u_on_u, v_on_v, n_on_n, u_on_n, v_on_n, v_on_u, u_on_v = bc.ic()

    E_num = []
    for timestep in range(pm.pm['nt']):
        if timestep % 2 != 0:
            sch.fbts_odd(n_on_n, u_on_u, v_on_v)
        elif timestep % 2 == 0:
            sch.fbts_even(n_on_n, u_on_u, v_on_v)
        #interpolate for u and v on n grid
        sch.inter_u_on_n(n_on_n, u_on_u, u_on_n)
        sch.inter_v_on_n(n_on_n, v_on_v, v_on_n)
        #calculate energy on n grid
        E_num.append(sch.calc_E(n_on_n, u_on_n, v_on_n))
        if timestep%100 == 0:
            print timestep

    v_on_u = sch.inter_v_on_u(u_on_u, v_on_v, v_on_u)
    u_on_v = sch.inter_u_on_v(u_on_u, v_on_v, u_on_v)

    if TaskB == True:
        u_on_vsouth = u_on_v[0,:].tolist()
        u_on_vsouth.insert(0,0)
        u_on_vsouth.append(0)

        v_on_uwest = v_on_u[:,0].tolist()
        v_on_uwest.insert(0,0)
        v_on_uwest.append(0)

        n_on_nmiddle = n_on_n[(pm.pm['nx']-1)/2,:].tolist()

        pt.plot_TaskB(u_on_vsouth, v_on_uwest, n_on_nmiddle, n_on_n)


    if TaskC1 == True:
        pt.plot_TaskC(E_num, 100)
    if TaskC2 == True:
        pt.plot_TaskC(E_num, 576, days = True)
    # print u_on_u
    # print n_on_n
    # print v_on_v
    # print E_num

    n_on_n_ana, u_on_n_ana, v_on_n_ana = bc.ic()[2:5]
    sch.analytic_solution(n_on_n_ana, u_on_n_ana, v_on_n_ana, n_on_n[(pm.pm['nx']-1)/2,0])
    # print u_on_n_ana
    # print v_on_n_ana
    # print n_on_n_ana

    u_diff_on_n = np.subtract(u_on_n, u_on_n_ana)
    v_diff_on_n = np.subtract(v_on_n, v_on_n_ana)
    n_diff_on_n = np.subtract(n_on_n, n_on_n_ana)

    if TaskD1 == True:
        print "The energy calculated from the difference fields is " \
              + str(sch.calc_E(n_diff_on_n, u_diff_on_n, v_diff_on_n)) + " J"

    if TaskD == True:
        pt.plot_TaskD(n_on_n_ana, n_diff_on_n)

####################################################

def main2(TaskE1 = False, TaskE2 = False):
    u_on_u, v_on_v, n_on_n, u_on_n, v_on_n, v_on_u, u_on_v = bc.ic()
    u_on_u_old, v_on_v_old, n_on_n_old, u_on_n_old, v_on_n_old, v_on_u_old, u_on_v_old = bc.ic()
    u_on_u_new, v_on_v_new, n_on_n_new, u_on_n_new, v_on_n_new, v_on_u_new, u_on_v_new = bc.ic()
    for timestep in range(pm.pm['nt']):
        if timestep % 2 != 0:
            sch.SL_fbts_odd(n_on_n_new, v_on_v_new, u_on_u_new, n_on_n, u_on_u, v_on_v, \
                            u_on_u_old, v_on_v_old, u_on_n, v_on_n, v_on_n_old, u_on_n_old)
            u_on_n_old = copy.copy(u_on_n)
            v_on_n_old = copy.copy(v_on_n)
            u_on_u_old = copy.copy(u_on_u)
            v_on_v_old = copy.copy(v_on_v)
            n_on_n = copy.copy(n_on_n_new)
            v_on_v = copy.copy(v_on_v_new)
            u_on_u = copy.copy(u_on_u_new)
            sch.inter_u_on_n(n_on_n_new, u_on_u_new, u_on_n)
            sch.inter_v_on_n(n_on_n_new, v_on_v_new, v_on_n)
        elif timestep % 2 == 0:
            sch.SL_fbts_even(n_on_n_new, v_on_v_new, u_on_u_new, n_on_n, u_on_u, v_on_v, \
                             u_on_u_old, v_on_v_old, u_on_n, v_on_n, v_on_n_old, u_on_n_old)
            u_on_n_old = copy.copy(u_on_n)
            v_on_n_old = copy.copy(v_on_n)
            u_on_u_old = copy.copy(u_on_u)
            v_on_v_old = copy.copy(v_on_v)
            n_on_n = copy.copy(n_on_n_new)
            v_on_v = copy.copy(v_on_v_new)
            u_on_u = copy.copy(u_on_u_new)
            sch.inter_u_on_n(n_on_n_new, u_on_u_new, u_on_n)
            sch.inter_v_on_n(n_on_n_new, v_on_v_new, v_on_n)
        print timestep

    v_on_u = sch.inter_v_on_u(u_on_u, v_on_v, v_on_u)
    u_on_v = sch.inter_u_on_v(u_on_u, v_on_v, u_on_v)

    if TaskE1 == True:
        u_on_vsouth = u_on_v[0,:].tolist()
        u_on_vsouth.insert(0,0)
        u_on_vsouth.append(0)

        v_on_uwest = v_on_u[:,0].tolist()
        v_on_uwest.insert(0,0)
        v_on_uwest.append(0)

        n_on_nmiddle = n_on_n[(pm.pm['nx']-1)/2,:].tolist()

        pt.plot_TaskE1(u_on_vsouth, v_on_uwest, n_on_nmiddle, n_on_n)

    n_on_n_ana, u_on_n_ana, v_on_n_ana = bc.ic()[2:5]
    sch.analytic_solution(n_on_n_ana, u_on_n_ana, v_on_n_ana, n_on_n[(pm.pm['nx']-1)/2,0])

    if TaskE2 == True:
        pt.plot_TaskE2(n_on_n, n_on_n_ana)

#main(TaskD1 = True)
#main2(TaskE1 = True)