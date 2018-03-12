import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
import cmath

def eigenvector_example(N,mu,A):

    #construct a Gaussian random matrix
    #A = np.random.randn(N, N)
    #rescale it so it's maximum value is one
    absmax = np.abs(A).max()
    A /= absmax

    #compute the eigenvalues and eigenvectors of A
    eigenvalues,eigenvectors = np.linalg.eig(A)

    #we want to plot the complex-valued eigenvalues.  we'll
    #consider their real part the x coordinate and the imaginary
    #part the y coordinate.

    plt.figure()
    #first plot the random matrix
    plt.subplot(2, 1, 1)
    plt.imshow(A, interpolation='nearest', aspect='auto', vmin=-1, vmax=1)
    plt.colorbar()
    plt.title('The A Matrix for the ENSO oscillator model')

    #then plot the eigenvalue spectrum
    plt.subplot(2, 1, 2)
    #plot the unit circle
    phase = np.linspace(-np.pi, np.pi, 200)
    xcirc = np.cos(phase)
    ycirc = np.sin(phase)
    plt.axhline(0.0, c='k')
    plt.axvline(0.0, c='k')
    plt.plot(xcirc, ycirc, 'k-')
    plt.plot(eigenvalues.real, eigenvalues.imag, 'ro')
    plt.axis('tight')
    plt.title('Eigenvalues, mu= '+str(mu))

#mu=0.01*i
nm=5
mu = pl.linspace(0,1,nm)
b0=2.5
c=1.0
r=0.25
b=b0*mu
gamma=0.75
R=gamma*b-c
alpha=0.125

d=R
e=gamma
f=-alpha*b
g=-r

A=np.empty([nm, 2, 2], float)
for i in range(0,nm,1):
    A[i,:,:]=np.array([[d[i],e],[f[i],g]])
#print A

eigenvalues,eigenvectors = np.linalg.eig(A)
Aeigen= abs(eigenvalues[:,0])

x = eigenvalues[:,0].real
y = eigenvalues[:,0].imag
[X,Y] = np.meshgrid(x,y)

pl.plot(mu, pl.real(eigenvalues), 'k-o', markersize=5)
#pl.plot(mu, pl.real(evals_neg), 'k-o', markersize=5)
pl.xlabel('$\mu$')
pl.ylabel('Re[$\lambda$]')
    
#eigenvector_example(2,mu,A)
#plt.show()
#plt.savefig('AllEigen.png')
    
#mu=np.arange(0., 1., 0.01)
#for i in range(1,100,1):
#for i in range(1,5,1):
#    mu=0.01*i
#    b0=2.5
#    c=1.0
#    r=0.25
#    b=b0*mu
#    gamma=0.75
#    R=gamma*b-c
#    alpha=0.125

#    d=R
#    e=gamma
#    f=-alpha*b
#    g=-r

#    A=np.empty([2, 2], float)
#    A=np.array([[d,e],[f,g]])

#    eigenvector_example(2,mu,A)
#    plt.show()
#    plt.savefig('/Users/vidale/Python_scripts/EigenValues_'+str(i)+'.png')
