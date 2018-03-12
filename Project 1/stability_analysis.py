"""
This script describes all the functions required for the stability analysis.
"""
import parameters as pm
import cmath

def eig_stability(dt, gamma, r, b0, c, alfa, t=0):
    """
    This function performs the stability analysis of the 4th order Runge Kutta scheme for the coupled
    equations (3) and (4) in the report. |A| is returned for both eigenvalues.
    :param t: time t, mew is constant if mew_ann = 0
    :param dt: timestep of model
    :param gamma: gamma (constant)
    :param r: r (constant)
    :param b0: b0 (constant)
    :param c: c (constant)
    :param alfa: alfa (constant)
    :return: |A| for both eigenvalues
    """
    mew = pm.mew(t)
    b = pm.b(b0, mew)
    caps_r = gamma*b - c
    eig_1 = (caps_r - r + cmath.sqrt((r - caps_r)**2 - 4*(gamma*alfa*b - caps_r*r)))/2.
    eig_2 = (caps_r - r - cmath.sqrt((r - caps_r)**2 - 4*(gamma*alfa*b - caps_r*r)))/2.

    v_1 = eig_1*dt
    v_2 = eig_2*dt

    caps_a1 = 1 + v_1 + 0.5*v_1**2 + 1/6.*v_1**3 + 1/24.*v_1**4
    caps_a2 = 1 + v_2 + 0.5*v_2**2 + 1/6.*v_2**3 + 1/24.*v_2**4
    print "The amplification factor is " + str(abs(caps_a1)) #same value for abs(caps_a2)

if __name__ == "__main__":
    eig_stability(0, pm.pm['dt'], pm.pm['gamma'], pm.pm['r'], pm.pm['b0'], pm.pm['c'], pm.pm['alfa'])