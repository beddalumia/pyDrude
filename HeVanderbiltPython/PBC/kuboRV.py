## Analitical Free Electrons ##

def kuboRV(w,N,L,eigenE,eta,EF,cutoff):
    
    from BraKet import BraKet
    from math import sqrt, factorial, sin, cos, pi
 
    if eta != 0: # Just to made sure we are safe...
 
        A = sqrt(2/L) # Normalization factor for eigenfunctions   
        def psi_n(x): return A*sin(pi/L * n * x)
        def psi_m(x): return A*sin(pi/L * m * x)     
        def vpsi_n(x): return 1j*A*cos(pi/L * n * x)*(pi/L)*n # \hat{v} applied to psi_n(x) 
        def vpsi_m(x): return 1j*A*cos(pi/L * m * x)*(pi/L)*m # \hat{v} applied to psi_m(x) 
    
        s = 100  # Sampling lenght for integrals: dx = L/s
    
        # Sum over all possibile excited (many-body) states compatible with Pauli principle
    
        n = 1
        result = 0
        
        while eigenE(n) <= EF: # "each occupied state..."
        
            m = n
            dw = eigenE(m) - eigenE(n)  # i.e. w_{mn}
        
            while dw < cutoff: #"...can be excited to a non-occupied state within the cutoff"
        
                dw = eigenE(m) - eigenE(n)
            
                if eigenE(m) > EF:
                
                    result += (BraKet(psi_n,vpsi_m,0,L,s)*BraKet(psi_m,vpsi_n,0,L,s))/(dw*(w-dw+1j*eta))
                    result += (BraKet(psi_m,vpsi_n,0,L,s)*BraKet(psi_n,vpsi_m,0,L,s))/(dw*(w+dw+1j*eta))
            
                m += 1
            
            n += 1
            
    return result*1j
 
## Numerical wavefunctions ## 

def numericalKuboRV(w,N,L,a,eigenE,psi,vpsi,domain,eta,EF,BandNum,cutoff):
    
    import numpy
    from numBraKet import numBraKet as BraKet
    from math import sqrt, factorial, sin, cos, pi
    
    if eta != 0: # Just to made sure we are safe...
    
        """newindex = numpy.argsort(eigenE)
        tempE = eigenE
        tempPsi = psi
        tempVpsi = vpsi
        for i in range(eigenE.shape[0]):
            eigenE[i]=tempE[newindex[i]]
            psi[i] = tempPsi[newindex[i]+1]
            vpsi[i] = tempVpsi[newindex[i]+1]
        print('E: ', eigenE)
        E = numpy.sort(eigenE)
        print('E: ', E)"""
        
        #print('E: ', eigenE)
        
        result = 0
        
        """print('Transizioni oblique disabled')
        kvec = numpy.linspace(-pi/a,pi/a,int(L/a))
        for BZindex in range(int(L/a)):
            
            k = kvec[BZindex]
            n = 1 + BZindex*BandNum
            
        # Sum over all possibile excited (many-body) states compatible with Pauli principle
    
        while eigenE[n-1] <= EF: # "each occupied state..."
        
            m = n
            dw = eigenE[m-1] - eigenE[n-1]  # i.e. w_{mn}
            
            while dw < cutoff: #"...can be excited to a non-occupied state within the cutoff"
                dw = eigenE[m-1] - eigenE[n-1]

                if eigenE[m-1] > EF: # and abs(M) > 0.1 :
                    
                    M = BraKet(psi[n],vpsi[m],0,L,domain)
                    result += abs(M)**2/(dw*(w-dw+1j*eta))
                    result += abs(M)**2/(dw*(w+dw+1j*eta))
            
                m += 1
            
            n += 1"""
        
        for i in range(eigenE.shape[0]):
            
            for j in range(eigenE.shape[0]):
                
                if eigenE[i] <= EF and eigenE[j] > EF:
                
                    dw = eigenE[j]-eigenE[i]
                    
                    if dw > 9.5 and dw < cutoff:
                        
                        M = BraKet(psi[i+1],vpsi[j+1],0,L,domain)
                        
                        if abs(M) > 0:
                            
                            result += abs(M)**2/((w-dw+1j*eta))
                            result += abs(M)**2/((w+dw+1j*eta))
    return result*1j    