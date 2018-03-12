"""
This script describes all the schemes required for the model. Only the Runge Kutta scheme is included in this case.
"""
import parameters as pm

def f(t, temp, h, gamma, eps, b0, c):
    """
    This function returns the Right-hand-side of equation (3) in the report evaluated at time t. Since the Runge Kutta
    scheme requires f evaluated at various times (t, t+0.5*dt, t+dt), only constant parameters are taken as arguments
    and other time dependent parameters (such as b, mew, caps_r, xi) are calculated at time t in this function.
    :param t: f evaluated at time t
    :param temp: east Pacific SST anomaly
    :param h: west Pacific ocean thermocline depth
    :param gamma: gamma (constant)
    :param eps: epsilon (constant)
    :param b0: b0 (constant)
    :param c: c (constant)
    :return: Right-hand-side of equation (3) in report
    """
    mew = pm.mew(t)
    b = pm.b(b0,mew)
    caps_r = gamma*b - c
    xi = pm.xi(t)
    return caps_r*temp + gamma*h - eps*(h+b*temp)**3 + gamma*xi

def g(t, temp, h, r, alfa, b0):
    """
    This function returns the Right-hand-side of equation (4) in the report evaluated at time t. Since the Runge Kutta
    scheme requires g evaluated at various times (t, t+0.5*dt, t+dt), only constant parameters are taken as arguments
    and other time dependent parameters (such as b, mew, xi) are calculated at time t in this function.
    :param t: g evaluated at time t
    :param temp: east Pacific SST anomaly
    :param h: west Pacific ocean thermocline depth
    :param r: r (constant)
    :param alfa: alpha (constant)
    :param b0: b0 (constant)
    :return: Right-hand-side of equation (4) in report
    """
    mew = pm.mew(t)
    b = pm.b(b0,mew)
    xi = pm.xi(t)
    return -r*h - alfa*b*temp - alfa*xi

#4th-order Runge Kutta
def rk4(dt, t, temp, h, gamma, eps, r, alfa, b0, c):
    """
    This function performs the Runge Kutta scheme for a single timestep. The functions f and g are based on the
    differential equations which require solving, and only constant parameters are required as arguments. Time dependent
    parameters are calculated in functions f  and g. The description of variables in the model is explained in the
    report.
    :param dt: timestep of model
    :param t: Runge Kutta done at time t
    :param temp: east Pacific SST anomaly
    :param h: west Pacific ocean thermocline depth
    :param gamma: gamma (constant) to pass into f
    :param eps: epsilon (constant) to pass into f
    :param r: r (constant) to pass into g
    :param alfa: alpha (constant) to pass into g
    :param b0: b0 (constant) to pass into f and g
    :param c: c (constant) to pass into f
    :return: temp and h at time t+dt using the 4th order Runge Kutta scheme.
    """
    k1 = dt*f(t, temp, h, gamma, eps, b0, c)
    m1 = dt*g(t, temp, h, r, alfa, b0)
    k2 = dt*f(t+0.5*dt, temp+0.5*k1, h+0.5*m1, gamma, eps, b0, c)
    m2 = dt*g(t+0.5*dt, temp+0.5*k1, h+0.5*m1, r, alfa, b0)
    k3 = dt*f(t+0.5*dt, temp+0.5*k2, h+0.5*m2, gamma, eps, b0, c)
    m3 = dt*g(t+0.5*dt, temp+0.5*k2, h+0.5*m2, r, alfa, b0)
    k4 = dt*f(t+dt, temp+k3, h+m3, gamma, eps, b0, c)
    m4 = dt*g(t+dt, temp+k3, h+m3, r, alfa, b0)

    tempnew = temp + (1/6.)*(k1+(2.*k2)+(2.*k3)+k4)
    hnew = h + (1/6.)*(m1+(2.*m2)+(2.*m3)+m4)

    return tempnew, hnew

