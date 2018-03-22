"""
This script describes all the plotting schemes for the various Tasks.
"""
import matplotlib.pyplot as plt
import parameters as pm
import numpy as np

def line_plots(u_on_vsouth, v_on_uwest, n_on_nmiddle):
    """
    This function allows for the plotting of the line plots required for Task B and E.
    """
    x_array = [(0.5+i)*pm.pm['d'] for i in range(pm.pm['nx']-1)]
    x_array.insert(0,0)
    x_array.append(pm.pm['L'])
    y_array = [(0.5+i)*pm.pm['d'] for i in range(pm.pm['nx']-1)]
    y_array.insert(0,0)
    y_array.append(pm.pm['L'])

    plt.plot(x_array, u_on_vsouth, label="$u$ on southern edge")
    plt.plot(y_array, v_on_uwest, label="$v$ on western edge")
    plt.xticks(rotation=90)
    plt.xlabel("Spatial distance (m)")
    plt.ylabel("$u$, $v$ and $\eta$ (ms$^{-1}$ and m)")
    x_array = [(0.5 + i) * pm.pm['d'] for i in range(pm.pm['nx'] - 1)]
    plt.plot(x_array, n_on_nmiddle, label="$\eta$ through middle")
    plt.legend(loc='upper right', bbox_to_anchor=(0.8, -0.4))
    plt.tight_layout()

def contour_plot(n_array):
    """
    This function plots the contour plots given an n array.
    """
    plt.xlabel("Spatial distance (m)")
    plt.xticks([0, 5, 10, 15, 20, 25], [0, 200000, 400000, 600000, 800000, 1000000], rotation=90)
    plt.ylabel("Spatial distance (m)")
    plt.yticks([0, 5, 10, 15, 20, 25], [0, 200000, 400000, 600000, 800000, 1000000])
    plt.contourf(n_array, cmap='bwr')
    plt.colorbar(label="$\eta$ (m)")
    plt.tight_layout()

def plot_TaskB(u_on_vsouth, v_on_uwest, n_on_nmiddle, n_on_n):
    """
    This function plots for Task B.
    """
    #line plots
    plt.subplot(1,2,1)
    plt.title("(a)")
    line_plots(u_on_vsouth, v_on_uwest, n_on_nmiddle)

    #contour plot
    plt.subplot(1,2,2)
    plt.title("(b)")
    contour_plot(n_on_n)
    plt.show()

def plot_TaskC(E_num, tick_space, days = False, hours = False):
    """
    This function plots for Task C.
    """
    t_array = range(pm.pm['nt'])
    plt.plot(t_array, E_num)
    plt.xlabel("Time (s)")
    #time axis in hours
    if hours == True:
        plt.xlabel("Time (hours)")
        plt.xticks(np.arange(0, pm.pm['nt'], tick_space), np.arange(0, 24), rotation=90)
    #time axis in days
    if days == True:
        plt.xlabel("Time (days)")
        plt.xticks(np.arange(0,pm.pm['nt'],tick_space),np.arange(0, 30),rotation=90)
    plt.ylabel("$E$ (J)")

def plot_TaskD(n_on_n_ana, n_diff_on_n):
    """
    This function plots for Task D.
    """
    #contour plot
    plt.subplot(1,2,1)
    plt.title("(a)")
    contour_plot(n_on_n_ana)

    #contour plot
    plt.subplot(1,2,2)
    plt.title("(b)")
    contour_plot(n_diff_on_n)

def plot_TaskD3(n_on_n_ana, n_on_n):
    """
    This function plots for Task D.
    """
    #contour plot
    plt.subplot(1,2,1)
    plt.title("(a)")
    contour_plot(n_on_n_ana)

    #contour plot
    plt.subplot(1,2,2)
    plt.title("(b)")
    contour_plot(n_on_n)

def plot_TaskD2(r1, r2, r3, error1, error2, error3):
    """
    This function plots for Task D.
    """
    #line plot
    plt.xlabel("Horizontal resolution (m)")
    plt.ylabel("Error in model (J)")
    plt.plot([r1,r2,r3],[error1,error2,error3])

def plot_TaskE1(u_on_vsouth, v_on_uwest, n_on_nmiddle, n_on_n):
    """
    This function plots for Task E.
    """
    #line plots
    plt.subplot(1,2,1)
    plt.title("(a)")
    line_plots(u_on_vsouth, v_on_uwest, n_on_nmiddle)

    #contour plot
    plt.subplot(1,2,2)
    plt.title("(b)")
    contour_plot(n_on_n)
    plt.show()

def plot_TaskE2(n_on_n, n_on_n_ana):
    """
    This function plots for Task E.
    """
    #contour plot
    plt.subplot(1,2,1)
    plt.title("(a)")
    contour_plot(n_on_n_ana)

    #contour plot
    plt.subplot(1,2,2)
    plt.title("(b)")
    contour_plot(n_on_n)