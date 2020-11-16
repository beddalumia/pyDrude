def DrawPolarization(eigenCs, N, L): # eigenCs being an array of projection coordinates in the free electron basis

    print('Computing the electron density...')

    import pylab
    import numpy
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from BraKet import BraKet
    from cmath import sqrt, exp, pi

    dx = L/100 # Sampling...
    
    ## Imperturbed eigenproblem: PLANE WAVES ##
    A = sqrt(1/a) # Normalization factor
    def q(n): return (2*pi*n)/L 
    def psi0_n(x): return A*exp(1j*q(n)*x)
    def psi0_m(x): return A*exp(1j*q(m)*x)
    def eigenE0(n): return (4 * pi**2 * n**2)/(2 * L**2) 
    
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
            
                n, m = i, j
            
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
        
            sum += numpy.conjugate(psi[n][k])*psi[n][k]
                
        result = numpy.append(result, sum)
    
    pylab.plot(x, 5*result)
    
    ## Excited-state density ##
    
    ascissa = numpy.array([])
    result = numpy.array([])
    for k in range(SamplingPts):
    
        sum = 0
        for n in range (N+1,2*N+1):
        
            sum += numpy.conjugate(psi[n][k])*psi[n][k]          
                
        result = numpy.append(result, sum)
    
    pylab.plot(x, 5*result)
    
    
    
    # Higher Excited-state density

    ascissa = numpy.array([])
    result = numpy.array([])
    for k in range(SamplingPts):
    
        sum = 0
        for n in range (2*N+1,3*N+1):
        
            sum += numpy.conjugate(psi[n][k])*psi[n][k]
                
        result = numpy.append(result, sum)
    
    pylab.plot(x, 5*result)
    

    return 0
