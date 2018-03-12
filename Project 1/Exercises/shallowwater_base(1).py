import numpy as np
#import theano.tensor as T
from pylab import figure, imshow, title, colorbar
import SWE_base_parameters

# Initial Conditions
n = SWE_base_parameters.nx
u_start = np.zeros((n,n)) # velocity in x direction - still water
v_start = np.zeros((n,n)) # velocity in y direction - still water

# eta (like height) will be uniform with a perturbation in the center
eta_start = np.ones((n,n)) # pressure deviation (like height)
x,y = np.mgrid[:n,:n]
droplet_x, droplet_y = n/2, n/2
rr = (x-droplet_x)**2 + (y-droplet_y)**2
eta_start[rr<10**2] = 1.1 # add a perturbation in pressure surface

# Parameters describing simulation PLV modified to match what is in the iPython Notebook
box_size = SWE_base_parameters.xlen
xlen = SWE_base_parameters.xlen
grid_spacing =  SWE_base_parameters.dx
g = SWE_base_parameters.g
dt = SWE_base_parameters.dt
dx = SWE_base_parameters.dx
nx = SWE_base_parameters.nx
Hhat = SWE_base_parameters.Hhat

def init():
    line.set_data([],[])
    return line, 

def animate(i):
    line.set_data(x, Ht[:,i])
    return line,

def initialise(x,IC=1): 
    H = np.zeros(nx)
    U = np.zeros(nx)  
    if IC == 1:
        H[:nx/2]=3.
        H[nx/2:]=1.
    
    elif IC == 2:
        x0=xlen/2.
        eps=0.5
        H = 2.0 + eps*np.exp(-(x-x0)**2/0.5)
        
    else:
        return ("invalid IC")
    return H,U

def boundary_conditions(H,U,BC=1):
    if BC == 1:
        #Porous boundary conditions NEEDS WORK TO MAKE USE OF PAST TIME STEP
        H[0]=H[1]
        H[-1]=H[-2]
        U[0]=U[1]
        U[-1]=U[-2]
        
    elif BC == 2:
        #Reflective boundary conditions
        H[0]=H[1]
        H[-1]=H[-2]    
        U[0]=0.
        U[-1]=0.
    else:
        return ("invalid BC")
    return H,U
    
def roll(x, shift, axis):
    """
    A numpy-theano agnostic version of the numpy.roll operator
    calls either numpy.roll or theano.tensor.roll depending on class

    See numpy.roll for usage
    """
    if isinstance(x, np.ndarray):
        return np.roll(x, shift, axis)
    if isinstance(x, T.basic.TensorVariable):
        return T.roll(x, shift, axis)
    raise NotImplementedError()

def spatial_derivative(A, axis=0):
    """
    PLV's comment: CAREFUL!!! THIS ONLY WORKS 100% CORRECTLY WITH PERIODIC BCs
    Compute derivative of array A using balanced finite differences
    Axis specifies direction of spatial derivative (d/dx or d/dy)

    dA[i] =  A[i+1] - A[i-1]   / 2
    ... or with grid spacing included ...
    dA[i]/dx =  A[i+1] - A[i-1]   / 2dx

    Used By:
        d_dx
        d_dy
    """
    return (roll(A, -1, axis) - roll(A, 1, axis)) / (grid_spacing*2.)

def d_dx_2D_CS(A):
    return spatial_derivative(A,1)
def d_dy_2D_CS(A):
    return spatial_derivative(A,0)

def d_dx_1D_CS(A):
    return spatial_derivative(A,0)

def d_dx_1D_CS_explicit(A):
    d_dx_ex = np.zeros(len(A))
    for i in range(1,len(A)-1,1):
        d_dx_ex[i] = (A[i+1]-A[i-1]) / (grid_spacing*2.)
    return d_dx_ex

def d_dx_1D_CS_python_specific(A):
    d_dx_ps = np.zeros(len(A))
    d_dx_ps[1:-1] = (A[2:] - A[:-2]) / (grid_spacing*2.)
    return d_dx_ps

def d_dt(eta, u, v, g, b=0):
    """
    http://en.wikipedia.org/wiki/Shallow_water_equations#Non-conservative_form
    """
    du_dt = -g*d_dx(eta) - b*u
    dv_dt = -g*d_dy(eta) - b*v

    H = 0#eta.mean() - our definition of eta includes this term
    deta_dt = -d_dx(u * (H+eta)) - d_dy(v * (H+eta))

    return deta_dt, du_dt, dv_dt

def step(eta, u, v, g, dt=dt):
    """
    Step forward eta, u, v one step in time of duration dt

    See Also:
        d_dt
    """
    deta_dt, du_dt, dv_dt = d_dt(eta, u, v, g)

    eta = eta + deta_dt * dt
    u = u + du_dt * dt
    v = v + dv_dt * dt

    return (eta, u, v)
