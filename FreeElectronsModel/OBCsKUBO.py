import pylab
from mpl_toolkits.mplot3d import Axes3D
import numpy
from math import *

# Retrieving physical parameters:

particleDensity = int(input('Please insert the particle density N/L: '))
cutoff = int(input('Please insert the desired energy cut-off [in units of Fermi energy]: '))

# Defining the thermodynamic limit:

length = numpy.array([])     # This is for plotting 
DW = numpy.array([])         # Drude Weight Vs L 
z1 = numpy.array([[]])
z2=z1

fig1 = pylab.figure('SurfaceReal')
ax1 = fig1.add_subplot(111, projection='3d')
fig2 = pylab.figure('SurfaceImag')
ax2 = fig2.add_subplot(111, projection='3d')
"""
length = numpy.array([])     # This is for plotting 
fsum = numpy.array([])       # f-sum rule Vs L 
"""
for L in range(2,13,2):
    
    N = particleDensity * L     # Initial N value
    print('N: ', N)
    
    def eigenE(n): return (pi**2 * n**2)/(2 * L**2)
    
    # We need to compute Fermi Energy for the system (m=1 and hbar=1)
    
    if (N%2) == 0:   # even N
    
        EF = eigenE(N/2)
    
    else:            # odd N 
    
        EF = eigenE((N+1)/2)
    
    # So now we can define the actual energy cut-off
    
    EnergyCutOff = EF*cutoff
    
    # Now we define the causality infinitesimal as to be greater than average level spacing
    
    eta = 10*EF/L
    print('EF: ', EF)
        
    # KUBO Formula for conductivity (e=1): g = <<r|v>> / L, and we'll call <<r|v>> =: kuboRV
    
    from kuboRV import kuboRV
    
    def g(w): return 1/L * kuboRV(w,N,L,eigenE,eta,EF,EnergyCutOff)
    
    # Plotting conductivity
    x = numpy.linspace(-eta*L,eta*L,101) # 10 linearly spaced numbers
    Rg = numpy.array([])
    Ig = numpy.array([])
    
    print('Wait! I am computing matrix elements...')
    
    for i in range(0,101): # Be Careful! Indexing starts with zero...
        
        w = x[i]
        y = g(w)
        Rg = numpy.append(Rg,y.real)
        Ig = numpy.append(Ig,y.imag)
    #print('Drude weight: ', Rg[10])
    
    pylab.figure('TLReal')
    pylab.plot(x,Rg)
    pylab.figure('TLImag')
    pylab.plot(x,Ig)
    
    z1 = numpy.append(z1,Rg)
    z2 = numpy.append(z2,Ig)
    
    """

    #Plotting f-sum rule Vs L

    dx = float(x[100]-x[0])/101   # We want to integrate Rg(w), we'll use the simple midpoint rule: SUM{dx*f(x_midpoint)}
    lenght = numpy.append(lenght,L)
    fsum = numpy.append(fsum,dx*sum(Rg))

pylab.plot(lenght,2*fsum/pi) 
    """
    
    #Plotting Drude Weight Vs L

    length = numpy.append(length,L)
    #DW = numpy.append(DW,g(0).real)
    
#pylab.plot(lenght,DW)  
X,Y = numpy.meshgrid(x,numpy.transpose(length))
Z1 = numpy.reshape(z1, X.shape)
Z2 = numpy.reshape(z2, X.shape)
    
Axes3D.plot_surface(ax1,X,Y,Z1, rstride=1, cstride=1)
Axes3D.plot_surface(ax2,X,Y,Z2, rstride=1, cstride=1)
    

pylab.show()




