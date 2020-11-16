def DrawSolutions(eigenEs, eigenCs, L): # eigenCs being an array of projection coordinates in the free electron basis

    print('Plotting given wavefunctions with energy levels:')
    
    import numpy
    import pylab
    from cmath import sqrt, cos, exp, pi

    # Imperturbed eigenproblem: PLANE WAVES 
    A = sqrt(1/a)
    def q(n): return (2*pi*n)/L 
    def psi0_n(x): return A*exp(1j*q(n)*x)
    def psi0_m(x): return A*exp(1j*q(m)*x)
    def eigenE0(n): return (4 * pi**2 * n**2)/(2 * L**2) 
    
    # We proceed to construct the wavefunctions from projection coordinates...
    
    basisDIM = eigenCs.shape[1]
    SamplingPts = 500
    x = numpy.linspace(0,L,SamplingPts)
    
    for i in range(2):
            
        psi_i_x = numpy.array([])
        level_i_x = numpy.array([])
        
        for k in range(SamplingPts):
            
            sum = 0
            
            for j in range(basisDIM):
            
                n, m = i, j
            
                sum += psi0_m(x[k]) * eigenCs[i][j]
                
            psi_i_x = numpy.append(psi_i_x, sum)
            level_i_x = numpy.append(level_i_x, eigenEs[i])
    
        y1 = level_i_x + 1j*level_i_x
        pylab.plot(x,y1, c='gray')
        y2 = y1 + psi_i_x
        pylab.plot(x,y2.real, c='red')
        pylab.plot(x,y2.imag, c='orange')
        print(eigenEs[i])
        
    return 0
            
