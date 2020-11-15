import pylab
import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from BraKet import BraKet
from math import sqrt, factorial, sin, cos, pi

# Retrieving physical parameters:

particleDensity = int(input('Please insert the particle density N/L: '))
L = int(input('Please insert the desired energy QW lenght: '))
N = particleDensity * L
dx = L/1000

# Defining the eigenproblem

A = sqrt(2/L) # Normalization factor for eigenfunctions   
def psi_n(x): return A*sin(pi/L * n * x)

# Ground-state density

x = 0
ascissa = numpy.array([])
result = numpy.array([])
while x <= L:

    sum = 0
    for n in range (1,N+1):
    
        sum += psi_n(x)**2
            
    ascissa = numpy.append(ascissa, x)
    result = numpy.append(result, sum)
    x += dx

pylab.plot(ascissa, result)

# Excited-state density

x = 0
ascissa = numpy.array([])
result = numpy.array([])
while x <= L:

    sum = 0
    for n in range (N+1,2*N+1):
    
        sum += psi_n(x)**2
            
    ascissa = numpy.append(ascissa, x)
    result = numpy.append(result, sum)
    x += dx

pylab.plot(ascissa, result)
"""
# Excited-state density

x = 0
ascissa = numpy.array([])
result = numpy.array([])
while x <= L:

    sum = 0
    for n in range (2*N+1,3*N+1):
    
        sum += psi_n(x)**2
            
    ascissa = numpy.append(ascissa, x)
    result = numpy.append(result, sum)
    x += dx

pylab.plot(ascissa, result)
"""
pylab.show()
