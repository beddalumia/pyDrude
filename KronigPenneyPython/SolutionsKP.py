def SolutionsKP(L, a, U0, cutoff): # L: Total length | a: Lattice parameter | U0: Perturbation strenght | cutoff: Energy cut-off
    
    print('Diagonalizing KP Hamiltonian...')
    
    import numpy
    from math import sqrt, sin, cos, pi
    from BraKet import BraKet
    
    s = 10  # Sampling grain for integrals
    
    ## Imperturbed eigenproblem ##
    A = sqrt(2/L)
    def psi0_n(x): return A*sin(pi/L * n * x)
    def psi0_m(x): return A*sin(pi/L * m * x)
    def eigenE0(n): return (pi**2 * n**2)/(2 * L**2)    

    n = 1

    while eigenE0(n) <= cutoff:
        
        n += 1
    
    basisDIM = n # We've established the basis dimension from the energy cut-off
    print('Basis dimension: ', basisDIM)
    
    ## Now we proceed constructing the hamiltionian matrix H_{nm} ##
    
    H = [[0 for i in range(basisDIM)] for j in range(basisDIM)]
    
    for i in range(basisDIM):
    
        for j in range(basisDIM):
        
            n, m = i+1, j+1 # Dumb python indexing...

            ## Kinetic (diagonal) part ##
            if i==j: H[i][i] += eigenE0(n)
            
            ## Kronig-Penney part ##
            H[i][j] += U0*BraKet(psi0_m, psi0_n, 0, a/4, s)         # Left Surface
            
            k = 3/4 * a                                         
            while k < (L - 1/4 * a):      
                
                 H[i][j] += U0*BraKet(psi0_m, psi0_n, k, k+a/2, s)  # Inner Lattice
                 
                 k = k + a 
                                                              
            H[i][j] += U0*BraKet(psi0_m, psi0_n, L-a/4, L, s)       # Right Surface        
            
            """NB: This is a modified implemention of a KP-lattice.
                   It's good because of the 'realistic' surface
                   behaviour of the solutions."""
                   
    ## Perturbed solutions ##
    eigenvalues, eigenvectors = numpy.linalg.eigh(H)
    
    return eigenvalues, eigenvectors
 