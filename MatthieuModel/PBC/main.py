import pylab
import numpy
from cmath import sqrt, exp, pi

## Retrieving physical parameters: ##

particleDensity = int(input('Please insert the particle density N/L: '))
cutoff = float(input('Please insert the desired energy cut-off [in units of Fermi energy]: '))
V0 = float(input('Please insert the desired perturbation strenght [in units of Fermi energy]:'))
W = float(input('Please insert the desired perturbation spread [in units of lattice parameter]:'))

a = 1
W *= a

for i in range(15*a,16*a,1):
    
    if particleDensity%2 == 0: L = 2*i+1
    else: L = 4*i+2
    N = particleDensity * L     # Initial N value
    print('N: ', N)
    
    ## Imperturbed eigenvalues: PLANE WAVES ##
    def eigenE0(n): return (4 * pi**2 * n**2)/(2 * L**2)
    
    # We need to compute Fermi Energy for the system (Be carefull of plane-wave DEGENERACY!)
    
    n = 0
    numParticles = 2 # Two particles on the |k=0> state!
    
    while numParticles < N:
        
        n += 1
        numParticles += 4 # Four particles on each |k!=0> state!
    
    nF = n
    EF = eigenE0(nF)
    
    # So now we can define the actual energy cut-off and perturbation strenght
    
    EFcutoff = EF*cutoff
    U0 = 6.168502750680849 # 5*EF for N/L=1 Free Electrons
    U0 = V0*EF
    print('HV step:', U0)
    
    ## Imperturbed eigenvectors: PLANE WAVES ##
    A = sqrt(1/L) # Normalization factor
    def q(n): return (2*pi*n)/L 
    def psi0_n(x): return A*exp(1j*q(n)*x)
    def psi0_m(x): return A*exp(1j*q(m)*x)
    
    ## Solving the lattice model ##
    
    from LatticeSolutions import LatticeSolutions 
    eigenEs, eigenCs = LatticeSolutions(L,a,U0,W,10*EFcutoff)
    
    """print('Energy levels: ', eigenEs)
    print('HV wavefunctions projected on free electron basis:', eigenCs)"""

    ## Plotting the lattice solutions ##
    
    from DrawSolutions import DrawSolutions
    #DrawSolutions(eigenEs, eigenCs, L)
    
    ## Plotting the electron density ##
    
    from DrawPolarization import DrawPolarization
    #DrawPolarization(eigenCs, N, L)

    ## Plotting the band structure ##
    
    from DrawBandDispersion import DrawBands
    #DrawBands(eigenEs, N, L, a)

    ## Drude weight as a Fermi-Surface integral ##
    
    from DrudePBC import ComputeDW
    Dw = ComputeDW(eigenEs, N, L)
    print('Normalized Drude Weight: ', Dw/pi)
    
    #pylab.scatter(L,Dw/pi, c='green', edgecolor='black')
    
    ## Regular Conductivty from Kubo Formula ##
    
    from Conductivity import PBCsConductivity
    #PBCsConductivity(N,L,eigenEs,eigenCs,EFcutoff)
    
    pylab.show()
