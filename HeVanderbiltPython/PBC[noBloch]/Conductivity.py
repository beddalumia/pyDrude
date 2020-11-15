def PBCsConductivity(N,L,eigenEs,eigenCs,cutoff):

    print('Kubo linear response running...')

    import pylab
    import numpy
    from numBraKet import numBraKet as BraKet
    from cmath import sqrt, exp, pi
    
    dx = 1/100 # Sampling...
    
    ## Defining the eigenproblem ##
    
    # Imperturbed eigenproblem: PLANE WAVES 
    A = sqrt(1/L)
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
            
                n, m = i, j # |n=0> is perimetted in PBCs!
            
                sum += psi0_m(x[k]) * eigenCs[i][j]
                
            psi_i_x = numpy.append(psi_i_x, sum)
    
        psi.append(psi_i_x)
        
    if len(psi) < 5*N:  print('WARNING! The diagonalization basis may be too small.')

    # We need to compute Fermi Energy for the *PERTURBED* system:

    Ek = numpy.array([])
    for n in range(0,eigenEs.shape[0],2): Ek = numpy.append(Ek, eigenEs[n])

    n = 0
    numParticles = 2 # Two particles on the |k=0> state!
    
    while numParticles < N:
        
        n += 1
        numParticles += 4 # Four particles on each |k!=0> state!
    
    nF = n
    EF = Ek[nF]

    ## Applying the velocity operator to the eigenfunctions ##
    
    """
    We need to apply the velocity operator \hat{v} = i * d/dx to the numerical
    eigenfunctions of the model Hamiltonian, i.e. we need to numerical evaluate
    the spatial derivative.
    """
    
    vpsi = [[]]
    
    print('Wait! I am applying the velocity operator to the energy eigenstates...')
    
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
    
    eta = 0.1*EF/L

    # KUBO Formula for conductivity (e=1): g = <<r|v>> / L, and we'll call <<r|v>> =: kuboRV
    
    from kuboRV import numericalKuboRV as kuboRV
    
    def g(w): return 2/L * kuboRV(w,N,L,eigenEs,psi,vpsi,x,eta,EF,cutoff)
    
    # Plotting conductivity
    #w = numpy.linspace(-3*eta,3*eta,101) # 100 linearly spaced numbers
    w = numpy.linspace(0,40,101) # 100 linearly spaced numbers
    Rg = numpy.array([])
    Ig = numpy.array([])
    
    print('Wait! I am computing matrix elements...')
    
    for i in range(0,101): # Be Careful! Indexing starts with zero...
        
        Rg = numpy.append(Rg,g(w[i]).real)
        #Ig = numpy.append(Ig,g(w[i]).imag)
    
    pylab.figure('RegularConductivity')
        
    pylab.plot(w,Rg)
    #pylab.plot(w,Ig)
    
    ## Reg-sum on Re(g) ##
    
    from numBraKet import numIntegration as integrate
    Rsum = integrate(Rg,0, 40, w)/(pi/2)
    #pylab.scatter(L,Rsum, c='orange', edgecolor='black')
    print('Reg-sum', Rsum)
    
    
    return Rsum
