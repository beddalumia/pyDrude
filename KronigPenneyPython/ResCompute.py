def ResCompute(N,L,eigenEs,eigenCs,cutoff):
    
    print('Computing poles and residues...')

    import pylab
    import numpy
    from numBraKet import numBraKet as BraKet
    from math import sqrt, sin, cos, pi
    
    dx = L/100 # Sampling...
    
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
    
    # We need to compute Fermi Energy for the *PERTURBED* system:
    
    if (N%2) == 0:   # even N
    
        EF = eigenEs[int(N/2)-1]
    
    else:            # odd N 
    
        EF = eigenEs[int((N+1)/2) - 1]

    
    ## Applying the velocity operator to the eigenfunctions ##
    
    """
    We need to apply the velocity operator \hat{v} = i * d/dx to the numerical
    eigenfunctions of the model Hamiltonian, i.e. we need to numerical evaluate
    the spatial derivative.
    """
    
    vpsi = [[]]
    
    print('Wait! I'm applying the velocity operator to the energy eigenstates...')
    
    for n in range(1, basisDIM+1):
        
        vpsi_i_x = numpy.array([])
    
        for k in range(1,SamplingPts):
            
            dpsi_dx = (psi[n][k] - psi[n][k-1]) / (x[k] - x[k-1]) 
            vpsi_i_x = numpy.append(vpsi_i_x, 1j*dpsi_dx)
            
        vpsi_i_x = numpy.append(vpsi_i_x, 0) # Border trick!
        
        vpsi.append(vpsi_i_x)
        
    if len(vpsi[1]) != len(psi[1]): print('ERROR! [in vpsi computation]')

    ## Scavenging poles and computing residues ##
    
    # Journey through all possibile excited (many-body) states compatible with Pauli principle
    
    n = 1
    freq = numpy.array([])
    singleRes = numpy.array([])
        
    while eigenEs[n-1] <= EF: # "each occupied state..."
        
            m = n
            dw = eigenEs[m-1] - eigenEs[n-1]  # i.e. w_{mn}
        
            while dw < cutoff: #"...can be excited to a non-occupied state within the cutoff"
            
                if eigenEs[m-1] > EF:
                
                    w = abs(dw)
                    print('Pole: ', w)
                    
                    print('Wait! I am computing matrix elements...')
                    
                    res = 2*(BraKet(psi[n],vpsi[m],0,L,x)*BraKet(psi[m],vpsi[n],0,L,x))/(dw*L)
                    print('Residue: ',res)
                    
                    freq = numpy.append(freq, w*L)
                    singleRes = numpy.append(singleRes,res)
                    pylab.scatter(L,w, c='yellow', edgecolors='black')
                    
                m += 1

		dw = eigenE(m) - eigenE(n)

            n += 1

    return 0
