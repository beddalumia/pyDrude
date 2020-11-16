def DrawBands(eigenEs, N, L, a): # eigenCs being an array of projection coordinates in the free electron basis

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
    
    # We proceed to construct the band-structure from eigenEs vector:
    
    ## Extended Zone Scheme ##
    
    """pylab.figure('Extended Bandstructure')
    
    k  = numpy.array([0])
    Ek = numpy.array([eigenEs[0]])  
    
    print('k: ', k[0])  
    
    for n in range(1,eigenEs.shape[0]-2):
        
        if n%2 != 0:
            
            k  = numpy.append(k,q(int(n/2)+1))
            print('k: ', q(int(n/2)+1))
            Ek = numpy.append(Ek, eigenEs[n+1])
            
        else: 
        
            k  = numpy.append(k,q(-int(n/2)))
            print('k: ', q(-int(n/2)))
            Ek = numpy.append(Ek, eigenEs[n])
        
        
    pylab.plot(k,Ek, c='black', linestyle='none', marker='o')"""
    
    ## Reduced Zone Scheme ##
    
    pylab.figure('Folded Bandstructure')
    
    Ncell = int(L/a)
    
    k  = numpy.array([0])
    Ek = numpy.array([eigenEs[0]])  
    
    print('k: ', k[0])  
    
    for n in range(1,Ncell):
        
        if n%2 != 0:
            
            k  = numpy.append(k,q(int(n/2)+1))
            print('k: ', q(int(n/2)+1))
            Ek = numpy.append(Ek, eigenEs[n+1])
            
        else: 
        
            k  = numpy.append(k,q(-int(n/2)))
            print('k: ', q(-int(n/2)))
            Ek = numpy.append(Ek, eigenEs[n])
            
    for n in range(Ncell+1,2*Ncell):
        
        if n%2 != 0:
            
            k  = numpy.append(k,2*pi/a - q(int(n/2)+1))
            print('k: ', q(int(n/2)+1))
            Ek = numpy.append(Ek, eigenEs[n+1])
            
        else: 
        
            k  = numpy.append(k, -q(-int(n/2)) - 2*pi/a)
            print('k: ', q(-int(n/2)))
            Ek = numpy.append(Ek, eigenEs[n])
            
    #pylab.plot(0, Ek[Ek.shape[0]-2], c='black', linestyle='none', marker='o')
        
    """for n in range(2*Ncell+1,3*Ncell):
        
        if n%2 != 0:
            
            k  = numpy.append(k,q(int(n/2)+1) - 2*pi/a)
            print('k: ', q(int(n/2)+1))
            Ek = numpy.append(Ek, eigenEs[n+1])
            
        else: 
        
            k  = numpy.append(k,q(-int(n/2)) + 2*pi/a)
            print('k: ', q(-int(n/2)))
            Ek = numpy.append(Ek, eigenEs[n])"""
        
        
    pylab.plot(k,Ek, c='black', linestyle='none', marker='o')

    ## EF Line ##

    # We need to compute Fermi Energy for the *PERTURBED* system:
    
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