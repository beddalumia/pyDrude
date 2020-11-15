## Matrix element numerical evaluation ##

def numBraKet(bra, ket, a, b, domain): # [a,b]: integration interval | domain: complete x sampling vector

    import numpy
    from numBraKet import numIntegration as integral
    
    integrand = numpy.conjugate(bra)*ket
    
    matrixElement = integral(integrand,a,b,domain)
    
    return matrixElement

## Trapezoidal Riemannian Sum ##
    
def numIntegration(f, a, b, x):  # f[]: Sampled function | [a,b]: integration interval | x[]: Sampled domain

    SamplingPts = x.shape[0]
    dx = float(x[SamplingPts-1] - x[0]) / SamplingPts
    
    # Sampling the interval...
    
    i = 0
    while x[i] < (a - dx): i += 1
    i_a = i
    while x[i] < (b - dx): i += 1
    i_b = i

    # Integrating...
    
    result = 0.5*f[i_a] + 0.5*f[i_b]
    for j in range(1,SamplingPts): result += f[i_a+j]
    result *= dx
    
    return result

