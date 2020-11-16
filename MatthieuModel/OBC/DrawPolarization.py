def DrawPolarization(eigenCs, N, L): # eigenCs being an array of projection coordinates in the free electron basis

    print('Computing the electron density...')

    import pylab
    import numpy
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from BraKet import BraKet
    from math import sqrt, sin, cos, pi

    dx = L/1000 # Sampling...
    
    ## Defining the eigenproblem ##
    
    # Imperturbed basis
    A = sqrt(2/L)
    def psi0_n(x): return A*sin(pi/L * n * x)
    def psi0_m(x): return A*sin(pi/L * m * x)
    def eigenE0(n): return (pi**2 * n**2)/(2 * L**2) 
    
    # We proceed to construct the wavefunctions from projection coordinates...
    
    basisDIM = eigenCs.shape[1]
    SamplingPts = int(L/dx)
    x = numpy.linspace(0,L,SamplingPts)
    psi = [[]]
    
    for i in range(basisDIM):
            
        psi_i_x = numpy.array([])
        
        for k in range(SamplingPts):
            
            sum = 0
            
            for j in range(basisDIM):
            
                n, m = i+1, j+1 # Dumb python indexing...
            
                sum += psi0_m(x[k]) * eigenCs[i][j]
                
            psi_i_x = numpy.append(psi_i_x, sum)
    
        psi.append(psi_i_x)
    if len(psi) < 5*N:  print('WARNING! The diagonalization basis may be too small...')

    ## Ground-state density ##
    
    ascissa = numpy.array([])
    result = numpy.array([])
    for k in range(SamplingPts):
    
        sum = 0
        for n in range (1,N+1):
        
            sum += (psi[n][k])**2
                
        result = numpy.append(result, sum)
    
    pylab.plot(x, 5*result)
    
    ## Excited-state density ##
    
    ascissa = numpy.array([])
    result = numpy.array([])
    for k in range(SamplingPts):
    
        sum = 0
        for n in range (N+1,2*N+1):
        
            sum += (psi[n][k])**2            
                
        result = numpy.append(result, sum)
    
    pylab.plot(x, 5*result)
    
    
    """
    # Higher Excited-state density

    ascissa = numpy.array([])
    result = numpy.array([])
    for k in range(SamplingPts):
    
        sum = 0
        for n in range (2*N+1,3*N+1):
        
            sum += (psi[n][k])**2
                
        result = numpy.append(result, sum)
    
    pylab.plot(x, result)
    """

    return 0
