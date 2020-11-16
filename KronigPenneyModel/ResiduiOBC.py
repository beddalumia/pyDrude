import pylab
import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from BraKet import BraKet
from math import sqrt, factorial, sin, cos, pi
 
# Retrieving physical parameters:

particleDensity = int(input('Please insert the particle density N/L: '))
cutoff = float(input('Please insert the desired energy cut-off [in units of Fermi energy]: '))
s = 100  # Sampling length for integrals: dx = L/s

# Verifying the f-sum rule with pole residues

"""
 
# Defining the thermodynamic limit:

length = numpy.array([])     # This is for plotting 
Res = numpy.array([])        # Residue Sum Vs L 

 
for L in range(2,100,2):
    
    N = particleDensity * L     # Initial N value
    print('N: ', N)
    
    def eigenE(n): return (pi**2 * n**2)/(2 * L**2)
    A = sqrt(2/L) # Normalization factor for eigenfunctions   
    def psi_n(x): return A*sin(pi/L * n * x)
    def psi_m(x): return A*sin(pi/L * m * x)     
    def vpsi_n(x): return 1j*A*cos(pi/L * n * x)*(pi/L)*n # \hat{v} applied to psi_n(x) 
    def vpsi_m(x): return 1j*A*cos(pi/L * m * x)*(pi/L)*m # \hat{v} applied to psi_m(x)     
    
    # We need to compute Fermi Energy for the system (m=1 and hbar=1)
    
    if (N%2) == 0:   # even N
    
        EF = eigenE(N/2)
    
    else:            # odd N 
    
        EF = eigenE((N+1)/2)
    
    # So now we can define the actual energy cut-off
    
    EnergyCutOff = EF*cutoff
    print('EF: ', EF) 
 
    
    # Sum over all possibile excited (many-body) states compatible with Pauli principle
    
    n = 1
    result = 0
        
    while eigenE(n) <= EF: # "each occupied state..."
        
            m = n
            dw = eigenE(m) - eigenE(n)  # i.e. w_{mn}
        
            while dw < cutoff: #"...can be excited to a non-occupied state within the cutoff"
            
                if eigenE(m) > EF:
                
                    result += 2*(BraKet(psi_n,vpsi_m,0,L,s)*BraKet(psi_m,vpsi_n,0,L,s))/(dw*L)
                    result += 2*(BraKet(psi_m,vpsi_n,0,L,s)*BraKet(psi_n,vpsi_m,0,L,s))/(dw*L)
            
                m += 1

		dw = eigenE(m) - eigenE(n)
            
            n += 1
            
    print('Somma dei residui: ',result)
    
    length = numpy.append(length,L)
    Res = numpy.append(Res,result)

pylab.plot(length,Res)  

"""

# Showing how the residue weight flows to zero frequency with increasing L

 
# Defining the thermodynamic limit:

length = numpy.array([])     # This is for plotting *SINGLE* residue value Vs L
matrix = numpy.array([]) 
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')

for L in range(2,100,2):
    
    N = particleDensity * L     # Initial N value
    length = numpy.append(length,L)
    
    def eigenE(n): return (pi**2 * n**2)/(2 * L**2)
    A = sqrt(2/L) # Normalization factor for eigenfunctions   
    def psi_n(x): return A*sin(pi/L * n * x)
    def psi_m(x): return A*sin(pi/L * m * x)     
    def vpsi_n(x): return 1j*A*cos(pi/L * n * x)*(pi/L)*n # \hat{v} applied to psi_n(x) 
    def vpsi_m(x): return 1j*A*cos(pi/L * m * x)*(pi/L)*m # \hat{v} applied to psi_m(x)     
    
    # We need to compute Fermi Energy for the system (m=1 and hbar=1)
    
    if (N%2) == 0:   # even N
    
        EF = eigenE(N/2)
    
    else:            # odd N 
    
        EF = eigenE((N+1)/2)
    
    # So now we can define the actual energy cut-off
    
    EnergyCutOff = EF*cutoff
    print('EF: ', EF) 
 
    
    # Journey through all possibile excited (many-body) states compatible with Pauli principle
    
    n = 1
    freq = numpy.array([])
    singleRes = numpy.array([])
        
    while eigenE(n) <= EF: # "each occupied state..."
        
            m = n
            dw = eigenE(m) - eigenE(n)  # i.e. w_{mn}
        
            while dw < cutoff: #"...can be excited to a non-occupied state within the cutoff"
            
                if eigenE(m) > EF:
                
                    w = abs(dw)
                    res = 2*(BraKet(psi_n,vpsi_m,0,L,s)*BraKet(psi_m,vpsi_n,0,L,s))/(dw*L)
                    
                    print('Pole: ', w)     
                    print('Residue: ',res)
                    freq = numpy.append(freq, w*L)
                    singleRes = numpy.append(singleRes,res)
                    #ax.scatter(w,L,res.real, c='red', edgecolors='gray')
                    pylab.scatter(L,res, c='yellow', edgecolors='black')
                    
                m += 1

		dw = eigenE(m) - eigenE(n)

            n += 1
           
            #matrix = numpy.append(matrix, singleRes)
            #print('matrix.shape: ', matrix.shape)
            #print('freq.shape: ', freq.shape)
            #print('length.shape: ', length.shape)
    
    #pylab.hist(freq*singleRes)
    #pylab.plot(freq,singleRes)  

#X,Y = numpy.meshgrid(freq,numpy.transpose(length))
#Z = numpy.reshape(matrix, X.shape)
    
#Axes3D.plot_wireframe(ax,X,Y,Z, rstride=10, cstride=10)
    
plt.show()
