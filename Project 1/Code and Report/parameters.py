"""
This script describes all the parameters (pm) required for the model.
"""
import math
import random

#functions for time dependent parameters
def b(b0, mew):
    """
    This function calculates b at time t, since mew is a function of t for Tasks D and E.
    :param b0: b0 (constant)
    :param mew: mew (at time t)
    :return: b (at time t)
    """
    return b0*mew

def caps_r(gamma, b, c):
    """
    This function calculates R at time t, since b is a function of t for Tasks D and E.
    :param gamma: gamma (constant)
    :param b: b (at time t)
    :param c: c (constant)
    :return: R (at time t)
    """
    return gamma*b - c

def mew(t):
    """
    This function calculates mew at time t, which is required in Tasks D and E. Mew is constant for Task A to C where
    the time varying component is set to 0 (mew_ann = 0) and mew (constant) is varied by changing mew_0. For Tasks D
    and E, the time varying component is included.
    :param t: time t, mew is constant if mew_ann = 0
    :return: mew (at time t)
    """
    mew = pm['mew_0']*(1 + pm['mew_ann']*math.cos(2*math.pi*t/float(pm['tau'])-5*math.pi/6.))
    return mew

def xi(t):
    """
    This function calculates xi at time t, which is the case in Tasks D and E. Xi is set to 0 by setting f_ann = 0 and
    f_ran = 0 for Tasks A to C. These components are included (non-zero) for Tasks D and E.
    :param t: time t
    :return: xi (at time t)
    """
    caps_w = random.randint(-100,101)/100.
    xi = pm['f_ann']*math.cos(2*math.pi*t/float(pm['tau'])) + pm['f_ran']*caps_w*pm['tau_cor']/float(pm['dt'])
    return xi

#dictionary for storing constant parameters
pm = {}
#scaling values
pm['Tdim'] = 7.5 #1 model SST anomaly unit represents 7.5K
pm['hdim'] = 150 #1 model thermocline depth unit represents 150m
pm['tdim'] = 2 #1 model time unit represents 2 months

#time values
pm['dt'] = 1/60. #1 day after dimensionalising
pm['nt'] = 1260  #one period is 42 months after dimensionalising

#initial conditions
pm['temp0'] = 0.15 #1.125K after dimensionalising
pm['h0'] = 0 #0m after dimensionalising

#model parameters
#constants
pm['b0'] = 2.5
pm['gamma'] = 0.75
pm['c'] = 1
pm['r'] = 0.25
pm['alfa'] = 0.125
pm['eps'] = 0
#parameters which modulate mew
pm['mew_0'] = 2/3.
pm['mew_ann'] = 0
pm['tau'] = 6  #12 months after dimensionalising
#parameters which modulate xi
pm['f_ann'] = 0
pm['f_ran'] = 0
pm['tau_cor'] = 0
