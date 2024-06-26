# Estrellas/EDOs

def systemNS(r, yV, arg):
    """ 
    [g, g', N, N', p, p', \[Rho]] -> [g0, g1, n0, n1, p0, p1, rho]
    """
    g0, n0, p0 = yV

    feps, = arg
    # eps = feps(p0)  # warning
    if r==0:
        f0, f1, f2 = 0., 0., 0.
    elif p0<0:
        # Sch. solution
        f0 = (g0-g0**3)/(2.*r)
        f1 = (n0*(g0**2-1))/(2.*r)
        f2 = 0
    else:
        eps = feps(p0)  # remove the warning
        f0 = (g0 + g0**3*(-1 + eps*r**2))/(2.*r)
        f1 = (n0*(-1 + g0**2*(1 + p0*r**2)))/(2.*r)
        
        n1 = f1
        f2 = -((n1*(eps + p0))/n0)
    return [f0, f1, f2]


def systemBS(r, yV, arg):
    """
    SISTEMA de EDO que describen una estrella de Bosones en el contexto de EKG
    
    L = k R - aB(CD[-a][PhiC]*CD[a][Phi]/2 + m^2 Phi PhiC)

    donde R: escalar de Ricci, k=c^4/(16 G pi), CD indica la derivada covariante y m la masa del campo.

    El sistema es adimensionalizado de tal forma que:
    rF = r/m, 
    PhiF = Mp/Sqrt[2pi] Phi,
    wF = m w

    las cantidades con F son las físicas (no adimensionales)
    
    COMENTARIOS:
    El sistema es dividido en r=0, r>0, ver notebook de mathematica

    IMPLEMENTACION:

    Sistema:
    g' = f0
    N' = f1
    Phi' = f2
    Phi'' = f3

    Variables:
    [g, g', N, N', Phi, Phi'] -> [g0, g1, n0, n1, p0, p1]

    In:
    r  -> ri 
    yV -> [g(ri), n(ri), p(ri), p1(ri)]
    arg -> [w, aB]

    Out:
    [f0, f1, f2, f3]  ->  [g', N', Phi', Phi'']
    """

    g0, n0, p0, p1 = yV
    w, aB = arg

    if p0>80 or p0<-80:
        #sys.exit('El perfil se indeterminó')
        f0, f1, f2, f3 = 0, 0, 0, 0
    elif r > 0:
        f0 = (g0 - g0**3)/(2.*r) + (aB*g0*r*(n0**2*p1**2 + g0**2*p0**2*(2*n0**2 + w**2)))/n0**2
        f1 = ((-1 + g0**2)*n0)/(2.*r) + (aB*r*(n0**2*p1**2 + g0**2*p0**2*(-2*n0**2 + w**2)))/n0
        f2 = p1
        
        g1, n1 = f0, f1
        f3 = (g1*p1)/g0 - (p1*(2*n0 + n1*r))/(n0*r) + g0**2*p0*(2 - w**2/n0**2)
    else:
        f0 = 0
        f1 = 0
        f2 = p1
        f3 = p0*(2-w**2/n0**2)/3
        
    return [f0, f1, f2, f3]