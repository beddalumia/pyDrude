def BraKet(bra, ket, a, b, n):
    from trapezoidalRiemannSum import trapezoidalRiemannSum
    def integrand(x): return bra(x).conjugate()*ket(x)
    matrixElement = trapezoidalRiemannSum(integrand,a,b,n)
    return matrixElement