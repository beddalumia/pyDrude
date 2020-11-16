def ComputeRes(N,L,eigenEs,eigenCs,cutoff):
    
    print('Computing poles and residues...')

    import pylab
    import numpy
    from numBraKet import numBraKet as BraKet
    from math import sqrt, sin, cos, pi
    
    dx = 1/100 # Sampling...
    
    gap = 7
    
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
    
       nmax = int(N/2)
    
    else:            # odd N 
    
        nmax = int((N+1)/2)
        
    EF = eigenEs[nmax-1]

    
    ## Applying the velocity operator to the eigenfunctions ##
    
    """
    We need to apply the velocity operator \hat{v} = i * d/dx to the numerical
    eigenfunctions of the model Hamiltonian, i.e. we need to numerical evaluate
    the spatial derivative.
    """
    
    vpsi = [[]]
    
    #print('Wait! I am applying the velocity operator to the energy eigenstates...')
    
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
    Dvect = numpy.array([])
    Dsum = 0
    Rsum = 0
    res = 0
    flag = 0
    
    color = numpy.array(['yellow', 'orange', 'red', 'brown', 'blue', 'green', 'black'])
    #for i in range(int(nmax)): color = numpy.append(color,color)
        
    while eigenEs[n-1] <= EF: # "each occupied state..."
        
        m = n
        dw = eigenEs[m-1] - eigenEs[n-1]  # i.e. w_{mn}
        
        while dw < cutoff: #"...can be excited to a non-occupied state within the cutoff"
            
            if eigenEs[m-1] > EF:
                
                w = abs(dw)
                #print('Pole: ', w)
                
                #print('Wait! I am computing matrix elements...')
                    
                M1 = BraKet(psi[n],vpsi[m],0,L,x)
                M2 = BraKet(vpsi[m],psi[n],0,L,x)
                #print('M1, M2 = ', M1, M2)
                
                if abs(M1) < 0.01 or abs(M2) < 0.01:
                    res = 0
                else:
                    res = 4*M1*M2/(dw*L)
                    if res > 1: res=0
                #print('Residue: ',res)
                colorindex = m-n-1
                colorindex = colorindex % color.shape[0]
                pylab.figure('Pole Frequencies')
                pylab.scatter(L,w, c = color[colorindex], edgecolors='black')#, alpha = (1 - (nmax-n)/(nmax/5)))
                pylab.figure('Pole Families')
                if w < gap: 
                    pylab.scatter(L,w, c = 'red', edgecolors='black')
                else: 
                    pylab.scatter(L,w, c='orange', edgecolors='black')
                pylab.figure('Pole Residues')
                pylab.scatter(L,res, c = color[colorindex], edgecolors='black')
                
                """
                How to use color and alpha:                
                - alpha = (1 - (nmax-n)/nmax) shows how deep in the Fermi sea is the excitation
                - c = (m-n)/cutoff identifies the "thermodynamic limit families" of poles
                """
                
                if flag == 0: 
                
                    FirstRes = res
                    res = 0
                    flag = 1
                    
                if w < gap: Dvect = numpy.append(Dvect,res)
                
                elif w < cutoff: Rsum += res
                    
            m += 1
                
            dw = eigenEs[m-1] - eigenEs[n-1] 

        n += 1
    
    Dsum = numpy.sum(Dvect)
    #Dsum += FirstRes
    Fsum = Dsum+Rsum
    Error = abs(1-Fsum)/Fsum
    
    """pylab.figure('Sum Rules')
    pylab.scatter(L,Dsum, c='red', marker='o', edgecolors='black')
    pylab.scatter(L,Rsum, c='orange', marker='o', edgecolors='black')
    pylab.scatter(L,Fsum, c='green', marker='s', edgecolors='black')
    pylab.figure('Relative Error on f-sum')
    pylab.scatter(L,Error, c='yellow', marker='s', edgecolors='black')"""
    print('Normalized Dw: ', Dsum)
    print('f-sum: ', Fsum)
    
    pylab.figure('Pole Frequencies')
    pylab.scatter(L,cutoff, c = 'black', edgecolors='none', marker = 'X', alpha = 0.3)
    pylab.scatter(L,gap, c = 'black', edgecolors='none', marker = 'X', alpha = 0.3)
    pylab.figure('Pole Families')
    pylab.scatter(L,cutoff, c = 'black', edgecolors='none', marker = 'X', alpha = 0.3)
    pylab.scatter(L,gap, c = 'black', edgecolors='none', marker = 'X', alpha = 0.3)

    return Dsum,Rsum, Fsum, Error