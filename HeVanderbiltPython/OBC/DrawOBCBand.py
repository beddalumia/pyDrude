def DrawBands(eigenEs, N, L): # eigenCs being an array of projection coordinates in the free electron basis

    print('Plotting levels:')
    
    import numpy
    import pylab
    from math import pi
    
    pylab.figure('OBCsLevels')

    """# We proceed to construct a band-like structure of the OBCs levels (usefull for debugging PBCs..)
    
    n = numpy.linspace(1,eigenEs.shape[0], eigenEs.shape[0])
    pylab.plot(n,eigenEs, color='black', linestyle='none', marker='o')"""
    
    # We draw energy levels as h-lines
    
    for n in range(eigenEs.shape[0]):
        
        pylab.axhline(y=eigenEs[n], color='black', linestyle='-', alpha=0.3)
        
    
    # We draw EF as a h-line:
    

    if (N%2) == 0:   # even N
    
        nF = int(N/2)
    
    else:            # odd N 
    
        nF = int((N+1)/2)
        
    if (N/L)%2 == 0: # Insulator!
    
        EF = eigenEs[nF-1] + 0.5 * (eigenEs[nF] - eigenEs[nF-1])
        
    else:            # Metal
    
        EF = eigenEs[nF-1]
        
    pylab.axhline(y=EF, color='r', linestyle='--')
    
    Ek = eigenEs
    
    if (N/L)%2 != 0: vF = (Ek[nF+1]-Ek[nF])/(pi/L)
    else: vF = (Ek[nF]-Ek[nF-1])/(pi/L)

    # And finally...
    Dw = 2*vF # From Eq.58 of RNC...[e = 1]
    print('Normalized Drude weight', Dw/pi)
    
    ## Graphical refinements ##
    
    pylab.ylim(0, 15)
    pylab.ylabel('Energy Levels [atomic units]')
    pylab.xticks([])    
    
    pylab.show()

    return 0