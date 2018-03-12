# Program to produce graphs showing the analytic and numerical stability of the
# recharge oscillator 
from __future__ import division
import numpy as np
import pylab as pl
import cmath

def Euler_stability():
    "Stability of rk4 for the linear recharge oscillator"
    
    x = np.arange(-4,2,0.01)
    y = np.arange(-3,3,0.01)
    [X,Y] = np.meshgrid(x,y)
    z = X + 1j*Y

    #CPL FOR YOUR EXERCISE
    #A = ????
    Aabs = abs(A)

    pl.contourf(X, Y, Aabs, pl.linspace(0,1,11))
#    pl.xlabel('$\lambda \Delta t$'); pl.ylabel('$\omega \Delta t$ ')
    pl.xlabel('$\Re(\lambda) \Delta t$'); pl.ylabel('$\Im(\lambda) \Delta t$')
    cbar = pl.colorbar()
    cbar.set_label("|A|", rotation=0)
 
def rk4_stability():
    "Stability of rk4 for the linear recharge oscillator"
    
    x = np.arange(-4,2,0.01)
    y = np.arange(-3,3,0.01)
    [X,Y] = np.meshgrid(x,y)
    z = X + 1j*Y

    A = 1 + z + 0.5*z**2 + 1/6*z**3 + 1/24*z**4
    Aabs = abs(A)

    pl.contourf(X, Y, Aabs, pl.linspace(0,1,11))
#    pl.xlabel('Re[$\lambda$]'); pl.ylabel('Im[$\lambda$]')
    pl.xlabel('$\Re(\lambda) \Delta t$'); pl.ylabel('$\Im(\lambda) \Delta t$')
    cbar = pl.colorbar()
    cbar.set_label("|A|", rotation=0)
 
def rk2_stability():
    "Stability of rk2 for the linear recharge oscillator"
    
    x = np.arange(-4,2,0.01)
    y = np.arange(-3,3,0.01)
    [X,Y] = np.meshgrid(x,y)
    z = X + 1j*Y

    A = 1 + z + 0.5*z**2
    Aabs = abs(A)

    pl.contourf(X, Y, Aabs, pl.linspace(0,1,11))
#    pl.xlabel('Re[$\lambda$]'); pl.ylabel('Im[$\lambda$]')
    pl.xlabel('$\Re(\lambda) \Delta t$'); pl.ylabel('$\Im(\lambda) \Delta t$')
    cbar = pl.colorbar()
    cbar.set_label("|A|", rotation=0)
  
 
def analytic_stability(r=0.25, Y=0.75, b0=2.5, c=1, a=0.125):
    "Stability of the analytic solution for the linear recharge oscillator"

    mu = pl.linspace(0,1,50)
    b = b0*mu
    R = Y*b-c

    evals_pos = pl.zeros((50), dtype=complex)
    evals_neg = pl.zeros((50), dtype=complex)

    for i in xrange(0,50):
        evals_pos[i] = 0.5*((R[i]-r) + cmath.sqrt((r-R[i])**2-4*(Y*a*b[i]-R[i]*r)))
        evals_neg[i] = 0.5*((R[i]-r) - cmath.sqrt((r-R[i])**2-4*(Y*a*b[i]-R[i]*r)))

    pl.plot(mu, pl.real(evals_pos), 'k-o', markersize=5)
    pl.plot(mu, pl.real(evals_neg), 'k-o', markersize=5)
    
    pl.xlabel('$\mu$')
    pl.ylabel('$\Re(\lambda)$')
    
    
