"""
This script contains the function to run the model and plot the trajectory of T and h, and time series of T.
"""
import numpy as np
import matplotlib.pyplot as plt

import parameters as pm
import schemes


def main():
    """
    This function finds the trajectory of T and h, and time series of T. The description of variables in this model is
    explained in the report.
    """
    #initialise arrays
    temp = np.zeros(pm.pm['nt']+1)
    h = np.zeros(pm.pm['nt']+1)
    t = np.zeros(pm.pm['nt']+1)
    for i in range(pm.pm['nt']+1):
        t[i] = pm.pm['dt']*i

    temp[0] = pm.pm['temp0'] #assign initial conditions
    h[0] = pm.pm['h0'] #assign initial conditions

    #repeat for nt number of timesteps
    for i in range(pm.pm['nt']):
        temp[i+1], h[i+1] = schemes.rk4(pm.pm['dt'], t[i], temp[i], h[i], pm.pm['gamma'], \
                                        pm.pm['eps'], pm.pm['r'], pm.pm['alfa'], pm.pm['b0'], pm.pm['c'])

    tempD = temp*pm.pm['Tdim'] #redimensionalise
    hD = h*pm.pm['hdim'] #redimensionalise
    tD = t*pm.pm['tdim'] #redimensionalise
    return tempD, hD, tD

def plotgraph(line_label=''):
    """
    This function plots the trajectory of T and h, and time series of T. The description of variables in this model is
    explained in the report.
    """
    plt.subplot(1, 2, 1)
    plt.title("(a)")
    plt.plot(main()[0], main()[1], label=line_label, linewidth=0.3)
    plt.xlabel('T (K)')
    plt.axvline(x=0, color='k', linestyle='--', linewidth=0.2)
    plt.xticks(rotation=70)
    plt.ylabel('h (m)')
    plt.axhline(y=0, color='k', linestyle='--', linewidth=0.2)

    plt.subplot(1, 2, 2)
    plt.title("(b)")
    plt.xlabel('time (months)')
    plt.xticks(rotation=70)
    plt.ylabel('T (K)')
    plt.plot(main()[2], main()[0], label=line_label, linewidth=0.3)
    plt.tight_layout()

def task_f_plot(line_label='',line_width=0.3):
    """
    This function plots for task F only.
    """
    plt.xlabel('time (months)')
    plt.xticks(rotation=70)
    plt.ylabel('T (K)')
    plt.plot(main()[2], main()[0], label=line_label, linewidth=line_width)
    plt.tight_layout()