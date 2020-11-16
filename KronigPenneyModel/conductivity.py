def conductivity(N,L,eigenEs,eigenCs,cutoff):

    print('Kubo linear response running...')

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
        
    if len(psi) < 5*N:  print('WARNING! The diagonalization basis may be too small.')

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

    ## KUBO FORMULA ##

    # We define the causality infinitesimal as to be greater than average level spacing
    
    eta = EF*10/L
        
    # KUBO Formula for conductivity (e=1): g = <<r|v>> / L, and we'll call <<r|v>> =: kuboRV
    
    from kuboRV import numericalKuboRV as kuboRV
    
    def g(w): return 1/L * kuboRV(w,N,L,eigenEs,psi,vpsi,x,eta,EF,cutoff)
    
    # Plotting conductivity
    #w = numpy.linspace(-3*eta,3*eta,101) # 100 linearly spaced numbers
    w = numpy.linspace(0,20,101) # 100 linearly spaced numbers
    Rg = numpy.array([])
    Ig = numpy.array([])
    
    print('Wait! I am computing matrix elements...')
    
    for i in range(0,101): # Be Careful! Indexing starts with zero...
        
        Rg = numpy.append(Rg,g(w[i]).real)
        #Ig = numpy.append(Ig,g(w[i]).imag)
        
    pylab.plot(w,Rg)
    #pylab.plot(w,Ig)
    
    return 0