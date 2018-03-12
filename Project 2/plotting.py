import matplotlib.pyplot as plt
import parameters as pm


def plot_TaskB(u_on_vsouth, v_on_uwest, n_on_nmiddle, n_on_n):
    x_array = [(0.5+i)*pm.pm['d'] for i in range(pm.pm['nx']-1)]
    x_array.insert(0,0)
    x_array.append(pm.pm['L'])
    y_array = [(0.5+i)*pm.pm['d'] for i in range(pm.pm['nx']-1)]
    y_array.insert(0,0)
    y_array.append(pm.pm['L'])

    plt.subplot(1, 2, 1)
    plt.title("(a)")
    plt.plot(x_array,u_on_vsouth)
    plt.plot(y_array, v_on_uwest)
    plt.xticks(rotation=90)
    plt.xlabel("spatial distance ($m$)")
    plt.ylabel("$u$, $v$ and $\eta$ ($ms^{-1}$ and $m$)")
    x_array = [(0.5 + i)*pm.pm['d'] for i in range(pm.pm['nx']-1)]
    plt.plot(x_array, n_on_nmiddle)
    plt.tight_layout()

    plt.subplot(1, 2, 2)
    plt.title("(b)")
    plt.contourf(n_on_n)
    plt.colorbar()

def plot_TaskC(E_num):
    t_array = range(pm.pm['nt'])
    plt.plot(t_array, E_num)