# Estrellas/Utilez

import sys
import numpy as np

def Bis(f, inter, Nit, error, eps):
    """
    The Bisection method
    In:
    f: it is a lambda or def object that represent a function f(x).
    inter: it is the interval [a, b] where we will find the roots.
    Nit: number of iterations. By default it is 1000.
    error: type of error used. The possibilities are {'dist', 'abs', 'rel', 'relMax'}. By default, it is 'dist'.
    eps: error bound, by default it is 1e-05.
    
    Out:
    roots: the found root
    """
    a, b = min(inter), max(inter)
    
    # p0 
    p0 = (a + b)/2
    sig = SigInt(f, [a, p0]) 
    if (sig>0):
        a = p0
    else:
        b = p0

    # pi
    for i in range(Nit):
        pi = (a + b)/2
        sig = SigInt(f, [a, pi])
        
        epsi = errorDic(error)(a, b) if error=='dist' else errorDic(error)(pi, p0)
        if abs(epsi)<=eps or sig==0.:
            root = pi
            break
        
        if (sig>0):
            a, p0 = pi, pi
        else:
            b, p0 = pi, pi
        
    if i==Nit-1:
        print('IMPORTANTE: la raiz encontrada no cumple el criterio de eps')
        root = pi
    return root

def roo_Bis(f, inter, Nit=1000, error='dist', eps=1e-05, Ndiv=100):
    """
    In:
    f: it is a lambda or def object that represent a function f(x).
    inter: it is the interval [a, b] where we will find the roots.
    Nit: number of iterations. By default it is 1000.
    error: type of error used. The possibilities are {'dist', 'abs', 'rel', 'relMax'}. By default, it is 'dist'.
    eps: error bound, by default it is 1e-05.
    Ndiv: number of subdivisions that the interval [a, b] will be divided, by default 100 is used.
    
    Out:
    roots: a list with the found roots
    """
    if not check1(f, inter):
        sys.exit("Revise los parámetros dados")
    
    interF = Findinterv(f, inter, Ndiv=Ndiv)
    if len(interF)==0:
        print('No se encontraron roots en el intervalo')
        roots = None
    else:
        roots = []
        for inter in interF:
            root = Bis(f, inter, Nit, error, eps)
            roots.append(root)
            
    return roots


def secv1(f, approx, Nit, error, eps):
    """
    The secant method
    In:
    f: it is a lambda or def object that represent a function f(x).
    approx: it is the approximations [p0, p1].
    Nit: number of iterations. By default it is 1000.
    error: type of error used. The possibilities are {'dist', 'abs', 'rel', 'relMax'}. By default, it is 'dist'.
    eps: error bound, by default it is 1e-05.
    
    Out:
    roots: the found root
    """
    p0, p1 = approx
    q0, q1 = f(p0), f(p1)
    
    for i in range(Nit):
        pi = p1-q1*(p1-p0)/(q1-q0)
        
        epsi = errorDic(error)(pi, p1) if error=='dist' else errorDic(error)(pi, p1)
        if abs(epsi)<=eps:
            root = pi
            break
        
        # update
        p0, q0 = p1, q1
        p1, q1 = pi, f(pi)
        
    if i==Nit-1:
        print('IMPORTANTE: la raiz encontrada no cumple el criterio de eps')
        root = pi
    return root

def roo_Secv1(f, inter, delt=2., Nit=1000, error='dist', eps=1e-05, Ndiv=100):
    """
    In:
    f: it is a lambda or def object that represent a function f(x).
    inter: it is the interval [a, b] where we will find the roots.
    delt: it is the value that is added to p0.
    Nit: number of iterations. By default it is 1000.
    error: type of error used. The possibilities are {'dist', 'abs', 'rel', 'relMax'}. By default, it is 'dist'.
    eps: error bound, by default it is 1e-05.
    Ndiv: number of subdivisions that the interval [a, b] will be divided, by default 100 is used.
    
    Out:
    roots: a list with the found roots
    """
    if not check1(f, inter):
        sys.exit("Revise los parámetros dados")
    
    interF = Findinterv(f, inter, Ndiv=Ndiv)
    if len(interF)==0:
        print('No se encontraron roots en el intervalo')
        roots = None
    else:
        roots = []
        for inter in interF:
            # first Bis
            p0 = Bis(f, inter, Nit, error, eps)
            p1 = p0 + delt
            root = secv1(f, [p0, p1], Nit, error, eps)
            roots.append(root)         
    return roots


def check1(f, inter):
    """
    Check of f is a object type function and inter have len two
    """
    test1 = True if callable(f) and len(inter)==2 else False
    return test1

def errorDic(met):
    """
    Error options
    """
    metodos = {'abs': lambda pn1, pn: abs(pn1-pn), 
               'rel': lambda pn1, pn: abs(pn1-pn)/abs(pn),
               'relMax': lambda pn1, pn: abs(pn1-pn)/abs(max([pn, pn1])),
               'dist': lambda a, b: (b-a)/2
               }
    return metodos[met]

def SigInt(f, inter):
    """
    Return the sign of the f(a)*f(b). In the case that this product is zero, return 0.
    """
    a, b = min(inter), max(inter)
    ya, yb = f(a), f(b)
    
    if (ya*yb<0):
        sig = -1 
    elif (ya*yb>0):
        sig = 1 
    else:
        sig = 0
    return sig

def Findinterv(f, inter, Ndiv=100):
    """
    encuentra intervalos con raices
    """
    a, b = min(inter), max(inter)
    
    if (a>0) and (b>0) and abs(b-a<1):
        interDisc = np.logspace(np.log10(a), np.log10(b), Ndiv)
    else:    
        interDisc = np.linspace(a, b, Ndiv)       
               
    interFinal = []
    for i in range(1, Ndiv):
        inter2 = [interDisc[i-1], interDisc[i]]
        if (SigInt(f, inter2) < 0):
            interFinal.append(inter2)      
    return interFinal

    