import pylab
import numpy
from math import sqrt, factorial, sin, cos, pi

## Retrieving physical parameters: ##

particleDensity = int(input('Please insert the particle density N/L: '))
cutoff = float(input('Please insert the desired energy cut-off [in units of Fermi energy]: '))
V0 = float(input('Please insert the desired perturbation strenght [in units of Fermi energy]:'))

a = 1

for L in range(2*a,15*a,2):
    
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
    
    EnergyCutOff = EF*cutoff
    U0 = EF*V0
    print('KP step:', U0)
    
    ## Imperturbed eigenvectors ##
    A = sqrt(2/L) # Normalization factor 
    def psi0_n(x): return A*sin(pi/L * n * x)
    def psi0_m(x): return A*sin(pi/L * m * x)
    
    ## Solving the Kronig-Penney model ##
    
    from SolutionsKP import SolutionsKP 
    eigenEs, eigenCs = SolutionsKP(L,a,U0, 99*cutoff)
    
    """print('Energy levels: ', eigenEs)
    print('KP wavefunctions projected on free electron basis:', eigenCs)"""

    ## Plotting the KP solutions ##
    
    from PrintKP import PrintKP
    #PrintKP(eigenEs, eigenCs, L)

    ## Plotting the KP change in polarization ##
    
    from KronigPenneyPolarization import PrintKronigPenneyPolarization
    #PrintKronigPenneyPolarization(eigenCs, N, L)
    
    ## Residue Analysis ##
    
    from ResCompute import ResCompute
    #ResCompute(N,L,eigenEs,eigenCs,EnergyCutOff)
    
    ## Conductivity ##
    
    from conductivity import conductivity
    conductivity(N,L,eigenEs,eigenCs,EnergyCutOff)
    
    pylab.show()