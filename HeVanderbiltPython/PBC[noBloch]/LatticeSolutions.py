def LatticeSolutions(L, a, U0, W, cutoff): # L: Total length | a: Lattice parameter | U0: Perturbation strenght | Spread: Perturbation width | cutoff: Energy cut-off
    
    print('Diagonalizing the perturbed Hamiltonian...')
    
    import numpy
    import pylab
    from cmath import sqrt, sin, cos, exp, pi
    from numBraKet import numBraKet as BraKet
    
    s = 10*L  # Sampling grain for integrals
    
    ## Imperturbed eigenproblem: PLANE WAVES ##
    A = sqrt(1/L) # Normalization factor
    def q(n): return (2*pi*n)/L 
    def psi0_n(x): return A*exp(1j*q(n)*x)
    def psi0_m(x): return A*exp(1j*q(m)*x)
    def eigenE0(n): return (4 * pi**2 * n**2)/(2 * L**2)   

    n = 0

    while eigenE0(n) <= cutoff:
        
        n += 1
    
    basisDIM = n # We've established the basis dimension from the energy cut-off
    print('Basis dimension: ', 2*basisDIM+1)
    
    ## Periodic potential definition ##
    
    def U_lat(x): 
    
        return U0*exp(-x**2/W**2) 
        
        """NB: x is centered on a lattice *barrier*,
               not on the atomic site!!!"""
    
    def U(x):
        
        U = 0
        
        for x_lat in range(0,L+a,a): 
            
            U += U_lat(x-x_lat)
            
        return U
        
    """NB: We have choosed this particular implemention
           for the lattice because of the ''realistic'' 
           surface behaviour of the solutions.
    
    x = numpy.linspace(0,L,s)
    f = numpy.vectorize(U, otypes=[numpy.complex])
    y = f(x)
    pylab.fill_between(x,y, facecolor='gray', edgecolor='black', alpha=0.3)"""
    
    ## Now we proceed constructing the hamiltionian matrix H_{nm} ##
    
    H = [[0 for i in range(2*basisDIM+1)] for j in range(2*basisDIM+1)]
    
    def Ket_m(x): return U(x)*psi0_m(x)
    
    Bra_n = numpy.vectorize(psi0_n, otypes=[numpy.complex])
    Ket_m = numpy.vectorize(Ket_m, otypes=[numpy.complex])
    
    for i in range(-basisDIM, basisDIM+1):
    
        for j in range(-basisDIM, basisDIM+1):
        
            n, m = i+1, j+1 # Dumb python indexing...

            ## Kinetic (diagonal) part ##
            
            if i==j: H[i][i] += eigenE0(n)
            
            ## Periodic potential part ##
            
            x = numpy.linspace(0,L,s)
            H[i][j] += BraKet(Bra_n(x), Ket_m(x), 0, L, x)    
            
    ## Perturbed solutions ##
    eigenvalues, eigenvectors = numpy.linalg.eigh(H)
    
    return eigenvalues, eigenvectors
 
