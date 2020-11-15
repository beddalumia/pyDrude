def PrintKP(eigenEs, eigenCs, L): # eigenCs being an array of projection coordinates in the free electron basis

    print('Plotting KP wavefunctions with energy levels:')
    
    import numpy
    import pylab
    from math import sqrt, cos, sin, pi

    # Imperturbed eigenproblem
    A = sqrt(2/L)
    def psi0_n(x): return A*sin(pi/L * n * x)
    def psi0_m(x): return A*sin(pi/L * m * x)
    def eigenE0(n): return (pi**2 * n**2)/(2 * L**2) 
    
    # We proceed to construct the wavefunctions from projection coordinates...
    
    basisDIM = eigenCs.shape[1]
    SamplingPts = 150
    x = numpy.linspace(0,L,SamplingPts)
    
    for i in range(3):
            
        psi_i_x = numpy.array([])
        level_i_x = numpy.array([])
        
        for k in range(SamplingPts):
            
            sum = 0
            
            for j in range(basisDIM):
            
                n, m = i+1, j+1 # Dumb python indexing...
            
                sum += psi0_m(x[k]) * eigenCs[i][j]
                
            psi_i_x = numpy.append(psi_i_x, sum)
            level_i_x = numpy.append(level_i_x, eigenEs[i])
    
        #y = level_i_x
        #pylab.plot(x,y)
        y = psi_i_x
        pylab.plot(x,y)
        print(eigenEs[i])
        
    return 0
            