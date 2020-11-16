def PBCsConductivity(N,L,a,ks,eigenEs,eigenCs,cutoff):

    print('Kubo linear response running...')

    import pylab
    import numpy
    from numBraKet import numBraKet as BraKet
    from cmath import sqrt, exp, pi
    
    dx = 1/100 # Sampling...
    
    ## Defining the eigenproblem ##
    
    # Imperturbed eigenproblem: PLANE WAVES 
    A = sqrt(1/a)
    def q(n): return (2*pi*n)/L 
    def psi0_n(x): return A*exp(1j*q(n)*x)
    def psi0_m(x): return A*exp(1j*q(m)*x)
    def eigenE0(n): return (4 * pi**2 * n**2)/(2 * L**2)   
    
    # We proceed to construct the wavefunctions from projection coordinates...
    
    BandNum = eigenCs.shape[1]
    SamplingPts = int(L/dx)
    x = numpy.linspace(0,L,SamplingPts)
    psi = [[]]
    
    for BZindex in range(int(L/a)):
        
        k = ks[BZindex]
        #print('k',k)
    
        for j in range(BandNum):
            
            psi_k_x = numpy.array([])
        
            for h in range(SamplingPts):
                
                sum = 0
                
                for i in range(BandNum):
                    
                    n = i*int(L/a) # Reciprocal lattice...
                    sum += psi0_n(x[h])*eigenCs[i+j][i] # uk in Fourier form...
                    
                sum *= exp(1j*k*x[h]) # Giving Psi_k...
                    
                psi_k_x = numpy.append(psi_k_x, sum)
        
            psi.append(psi_k_x)
        
    if len(psi) < 5*N:  print('WARNING! The diagonalization basis may be too small.')

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
        
    ## Applying the velocity operator to the eigenfunctions ##
    
    """
    We need to apply the velocity operator \hat{v} = i * d/dx to the numerical
    eigenfunctions of the model Hamiltonian, i.e. we need to numerical evaluate
    the spatial derivative.
    """
    
    vpsi = psi
    
    print('Wait! I am applying the velocity operator to the energy eigenstates...')
    
    """for n in range(1, len(psi)):
        
        vpsi_i_x = numpy.array([])
    
        for k in range(1,SamplingPts):
            
            dpsi_dx = (psi[n][k] - psi[n][k-1]) / (x[k] - x[k-1]) 
            vpsi_i_x = numpy.append(vpsi_i_x, 1j*dpsi_dx)
            
        vpsi_i_x = numpy.append(vpsi_i_x, 0) # Border trick!
        
        vpsi.append(vpsi_i_x)"""
        
    for n in range(1, len(psi)):
        
        vpsi[n] = 1j*numpy.gradient(psi[n])
    
    for i in range(len(psi)):   
        if len(vpsi[i]) != len(psi[i]): print('ERROR! [in vpsi computation]')

    ## KUBO FORMULA ##

    # We define the causality infinitesimal as to be greater than average level spacing
    
    eta = 0.1*EF/L

    # KUBO Formula for conductivity (e=1): g = <<r|v>> / L, and we'll call <<r|v>> =: kuboRV
    
    from kuboRV import numericalKuboRV as kuboRV
    
    def g(w): return 2/L * kuboRV(w,N,L,a,eigenEs,psi,vpsi,x,eta,EF,BandNum,cutoff)
    
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
    print('Reg-sum', Rsum)
    
    
    return Rsum
