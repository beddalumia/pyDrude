def ComputeDW(eigenEs, N, L): 

    import numpy
    import pylab
    from math import pi
    
    def q(n): return (2*pi*n)/L 
    
    #First of all: create a 'good' Ek vector (forcing degeneration...)
    
    Ek = numpy.array([eigenEs[0]])
    
    for n in range(1,eigenEs.shape[0]-1,2): Ek = numpy.append(Ek, eigenEs[n+1])

    # We need to compute Fermi Energy for the *PERTURBED* system:
    
    n = 0
    numParticles = 2 # Two particles on the |k=0> state!
    
    while numParticles < N:
        
        n += 1
        numParticles += 4 # Four particles on each |abs(k)!=0> state!
    
    nF = n
    EF = Ek[nF]
    
    # Then we need the Fermi velocity: vF = dEk/dk @ kF [hbar = 1]
    
    if (N/L)%2 != 0: vF = (Ek[nF+1]-Ek[nF])/((2*pi)/L)
    else: vF = (Ek[nF]-Ek[nF-1])/((2*pi)/L)

    print('Debug on numerical differentiation of Ek: vF, kF = ', vF, q(nF))
    
    # And finally...
    
    Dw = 2*vF # From Eq.58 of RNC...[e = 1]
    
    if Dw > 5*N/L: Dw = 0
    
    return Dw