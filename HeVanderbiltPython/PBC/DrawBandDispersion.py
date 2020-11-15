def DrawBands(ks, eigenEs, N, L, a, cutoff): 

    print('Plotting band dispersion:')
    
    import numpy
    import pylab
    from cmath import sqrt, cos, exp, pi
    
    # Imperturbed eigenproblem: PLANE WAVES 
    A = sqrt(1/a)
    def q(n): return (2*pi*n)/L 
    def psi0_n(x): return A*exp(1j*q(n)*x)
    def psi0_m(x): return A*exp(1j*q(m)*x)
    def eigenE0(n): return (4 * pi**2 * n**2)/(2 * L**2) 
    
    n = 0

    while eigenE0(n) <= cutoff:
        
        n += 1
    
    basisDIM = n # We've established the basis dimension from the energy cut-off
    print('Basis dimension: ', 2*basisDIM+1)
    
    # We proceed to construct the band-structure from eigenEs vector:
    
    BandNum = 2*int((2*basisDIM+1)*a/L)
    print('Num. of Bands: ', BandNum)
    j = 0
    
    for BZindex in range(int(L/a)):
        
        k = ks[BZindex]
        #print('k: ', k)
        
        for i in range(BandNum-1):
            
            pylab.plot(k,eigenEs[i+j], c='black', linestyle='none', marker='o')
            
        j += BandNum-1
        
            

    ## EF Line ##

    # We need to compute Fermi Energy for the *PERTURBED* system:
    
    Ek = numpy.sort(eigenEs)
    
    n = 0
    numParticles = 2 
    
    while numParticles < N:
        
        n += 1
        numParticles += 2 
    
    nF = n
    print('N:',numParticles)
    
    if (N/L)%2 == 0: # Insulator!
    
        EF = Ek[nF] + 0.5 * (Ek[nF+1] - Ek[nF])
        
    else:            # Metal
    
        EF = Ek[nF]
    
    pylab.axhline(y=EF, color='r', linestyle='--')
    EF = Ek[nF]

    ## Graphical refinements ##
    
    pylab.axvline(x=pi/a, color='black', linestyle='-')
    pylab.axvline(x=0, color='black', linestyle='--')
    pylab.axvline(x=-pi/a, color='black', linestyle='-')
    pylab.xlim(-3.2/a, 3.2/a)
    pylab.xlabel('Crystal Momentum')
    pylab.ylabel('Energy [atomic units]')
    pylab.xticks([])
    
    pylab.show()
    
    return 0