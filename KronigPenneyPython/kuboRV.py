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
 
## Numerical KP wavefunctions ## 

def numericalKuboRV(w,N,L,eigenE,psi,vpsi,domain,eta,EF,cutoff):
    
    from numBraKet import numBraKet as BraKet
    from math import sqrt, factorial, sin, cos, pi
 
    if eta != 0: # Just to made sure we are safe...

        # Sum over all possibile excited (many-body) states compatible with Pauli principle
    
        n = 1
        result = 0
        
        while eigenE[n-1]<= EF: # "each occupied state..."
        
            m = n
            dw = eigenE[m-1] - eigenE[n-1]  # i.e. w_{mn}
        
            while dw < cutoff: #"...can be excited to a non-occupied state within the cutoff"
        
                dw = eigenE[m-1] - eigenE[n-1]
            
                if eigenE[m-1] > EF:
                
                    result += (BraKet(psi[n],vpsi[m],0,L,domain)*BraKet(psi[m],vpsi[n],0,L,domain))/(dw*(w-dw+1j*eta))
                    result += (BraKet(psi[m],vpsi[n],0,L,domain)*BraKet(psi[n],vpsi[m],0,L,domain))/(dw*(w+dw+1j*eta))
            
                m += 1
            
            n += 1
            
    return result*1j    