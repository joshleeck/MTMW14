import matplotlib.pyplot as plt
import parameters as pm
import numpy as np


def plot_TaskB(u_on_vsouth, v_on_uwest, n_on_nmiddle, n_on_n):
    x_array = [(0.5+i)*pm.pm['d'] for i in range(pm.pm['nx']-1)]
    x_array.insert(0,0)
    x_array.append(pm.pm['L'])
    y_array = [(0.5+i)*pm.pm['d'] for i in range(pm.pm['nx']-1)]
    y_array.insert(0,0)
    y_array.append(pm.pm['L'])

    plt.subplot(1,2,1)
    plt.title("(a)")
    plt.plot(x_array,u_on_vsouth, label="$u$ on southern edge")
    plt.plot(y_array, v_on_uwest, label="$v$ on western edge")
    plt.xticks(rotation=90)
    plt.xlabel("spatial distance ($m$)")
    plt.ylabel("$u$, $v$ and $\eta$ ($ms^{-1}$ and $m$)")
    x_array = [(0.5 + i)*pm.pm['d'] for i in range(pm.pm['nx']-1)]
    plt.plot(x_array, n_on_nmiddle, label="$\eta$ through middle")
    plt.legend(loc='upper right', bbox_to_anchor=(0.8, -0.4))
    plt.tight_layout()

    plt.subplot(1,2,2)
    plt.title("(b)")
    plt.xlabel("spatial distance ($m$)")
    plt.xticks([0, 5, 10, 15, 20, 25],[0,200000,400000,600000,800000,1000000], rotation=90)
    #plt.xlim(0,50)
    #plt.ylim(0,50)
    plt.ylabel("spatial distance ($m$)")
    plt.yticks([0, 5, 10, 15, 20, 25],[0,200000,400000,600000,800000,1000000])
    plt.contourf(n_on_n)
    plt.colorbar(label="$\eta$ (m)")
    plt.tight_layout()

def plot_TaskC(E_num, tick_space, days = False):
    t_array = range(pm.pm['nt'])
    plt.plot(t_array, E_num)
    plt.xlabel("time (s)")
    if days = True:
        plt.xlabel("time (days)")
        plt.xticks(np.arange(0,pm.pm['nt'],tick_space),np.arange(0, 30),rotation=90)
    plt.ylabel("$E$ (J)")

def plot_TaskD(n_on_n_ana, n_diff_on_n):
    plt.subplot(1,2,1)
    plt.title("a")
    plt.contourf(n_on_n_ana)
    plt.colorbar()

    plt.subplot(1,2,2)
    plt.title("b")
    plt.contourf(n_diff_on_n)
    plt.colorbar()

def plot_TaskE1(u_on_vsouth, v_on_uwest, n_on_nmiddle, n_on_n):
    x_array = [(0.5+i)*pm.pm['d'] for i in range(pm.pm['nx']-1)]
    x_array.insert(0,0)
    x_array.append(pm.pm['L'])
    y_array = [(0.5+i)*pm.pm['d'] for i in range(pm.pm['nx']-1)]
    y_array.insert(0,0)
    y_array.append(pm.pm['L'])

    plt.subplot(1,2,1)
    plt.title("(a)")
    plt.plot(x_array,u_on_vsouth)
    plt.plot(y_array, v_on_uwest)
    plt.xticks(rotation=90)
    plt.xlabel("spatial distance ($m$)")
    plt.ylabel("$u$, $v$ and $\eta$ ($ms^{-1}$ and $m$)")
    x_array = [(0.5 + i)*pm.pm['d'] for i in range(pm.pm['nx']-1)]
    plt.plot(x_array, n_on_nmiddle)
    plt.tight_layout()

    plt.subplot(1,2,2)
    plt.title("(b)")
    plt.contourf(n_on_n)
    plt.colorbar()
    plt.show()

def plot_TaskE2(n_on_n, n_on_n_ana):
    plt.subplot(1,2,1)
    plt.title("a")
    plt.contourf(n_on_n)
    plt.colorbar()

    plt.subplot(1,2,2)
    plt.title("b")
    plt.contourf(n_on_n_ana)
    plt.colorbar()