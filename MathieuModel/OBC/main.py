import pylab
import numpy
from math import sqrt, factorial, sin, cos, pi

## Retrieving physical parameters: ##

particleDensity = int(input('Please insert the particle density N/L: '))
cutoff = float(input('Please insert the desired energy cut-off [in units of Fermi energy]: '))
V0 = float(input('Please insert the desired perturbation strenght [in units of Fermi energy]:'))
#W = float(input('Please insert the desired perturbation spread [in units of lattice parameter]:'))

a = 1
#W *= a

for L in range(2*a,40*a,2):

    N = particleDensity * L     # Initial N value
    print('N: ', N)
    
    ## Imperturbed eigenvalues ##
    def eigenE0(n): return (pi**2 * n**2)/(2 * L**2)
    
    # We need to compute Fermi Energy for the system (m=1 and hbar=1)
    
    if (N%2) == 0:   # even N
    
        EF = eigenE0(N/2)
    
    else:            # odd N 
    
        EF = eigenE0((N+1)/2)
    
    # So now we can define the actual energy cut-off and perturbation strenght
    
    EFcutoff = EF*cutoff
    U0 = EF*V0
    print('HV step:', U0)
    
    ## Imperturbed eigenvectors ##
    A = sqrt(2/L) # Normalization factor 
    def psi0_n(x): return A*sin(pi/L * n * x)
    def psi0_m(x): return A*sin(pi/L * m * x)
    
    ## Solving the lattice model ##
    
    from LatticeSolutions import LatticeSolutions 
    eigenEs, eigenCs = LatticeSolutions(L,a,U0,10*EFcutoff)
    
    """print('Energy levels: ', eigenEs)
    print('HV wavefunctions projected on free electron basis:', eigenCs)"""

    ## Plotting the lattice solutions ##
    
    from DrawSolutions import DrawSolutions
    #DrawSolutions(eigenEs, eigenCs, L)
    
    from DrawOBCBand import DrawBands
    #DrawBands(eigenEs, N, L)

    ## Plotting the change in polarization ##
    
    from DrawPolarization import DrawPolarization
    #DrawPolarization(eigenCs, N, L)
    
    ## Residue Analysis ##
    
    from ComputeRes import ComputeRes
    ComputeRes(N,L,eigenEs,eigenCs,EFcutoff)
    
    ## Conductivity ##
    
    from Conductivity import Conductivity
    #Conductivity(N,L,eigenEs,eigenCs,EFcutoff)
    
    pylab.show()
