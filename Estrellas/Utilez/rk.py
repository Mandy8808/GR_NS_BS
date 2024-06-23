# Estrellas/Utilez

import numpy as np


# Runge Kutta
def RKMet(func, data, arg=None, rk='rk5'):
    """
    Implementacion del metodo de  Runge Kutta de orden 4 y 5

    In:
    func -> Sistema de EDO, debe tener la estructura:
            (x, yval, a, b, ..) donde a, b, .. son argumentos que se han de pasar como listas a la variable arg
    
    data -> Es una lista [y0v, a, b, npt] que proporciona la información necesaria para la implementación del RK donde
            y0v ->  es una lista con las condiciones de frontera en r=a
            a, b -> rmin, rmax
            npt -> son los números de puntos usados para generar h = (b-a)/(npt-1)
    arg -> argumentos que se le pasarán al sistema de ecuaciones
    rk -> el método que se utilizará puede ser 'rk4' o 'rk5', por defecto está el 'rk5'

    Out:
    xval -> un arreglo discreto de valores de ri, donde se itero  
    yval -> una lista de arreglos discretos correspondientes a las iteraciones con la misma estructura que y0v. Cada fila corresponde a un ri.
    """
    if arg:
        f = lambda x, yv: np.array(func(x, yv, arg))
    else:
        f = lambda x, yv: np.array(func(x, yv))
    
    if rk == 'rk5':
        xval, yval = rk5Gene(f, data)
    else:
        xval, yval = rk4Gene(f, data)
    
    return xval, yval
    

def rk4Gene(f, data):
    """
    Implementacion del metodo de  Runge Kutta de orden 4

    In:
    func -> Sistema de EDO, debe tener la estructura:
            (x, yval, a, b, ..) donde a, b, .. son argumentos que se han de pasar como listas a la variable arg
    
    data -> Es una lista [y0v, a, b, npt] que proporciona la información necesaria para la implementación del RK donde
            y0v ->  es una lista con las condiciones de frontera en r=a
            a, b -> rmin, rmax
            npt -> son los números de puntos usados para generar h = (b-a)/(npt-1)
    
    Out:
    xval -> un arreglo discreto de valores de ri, donde se itero  
    yval -> una lista de arreglos discretos correspondientes a las iteraciones con la misma estructura que y0v. Cada fila corresponde a un ri.
    """
    y0v, a, b, npt = data

    # h val
    h = (b-a)/(npt-1)
    xval = a + np.arange(npt)*h
    nEq = len(y0v)
    yval = np.zeros((npt, nEq))  # crea una matriz (npt, nEq)
    
    yv = np.copy(y0v)
    for j, xi in enumerate(xval):
        yval[j, :] = yv
        k0 = h*f(xi, yv)  # importante como se define la función, tiene que ser de la forma f(x, yv) donde yv=[y0, y1, ...]
        k1 = h*f(xi+h/2, yv+k0/2)
        k2 = h*f(xi+h/2, yv+k1/2)
        k3 = h*f(xi+h, yv+k2)
        yv = yv + (k0 + 2*k1 + 2*k2 + k3)/6
        #yv += (k0 + 2*k1 + 2*k2 + k3)/6
           
    return xval, yval

def rk5Gene(f, data):
    """
    Implementacion del metodo de  Runge Kutta de orden 5

    In:
    func -> Sistema de EDO, debe tener la estructura:
            (x, yval, a, b, ..) donde a, b, .. son argumentos que se han de pasar como listas a la variable arg
    
    data -> Es una lista [y0v, a, b, npt] que proporciona la información necesaria para la implementación del RK donde
            y0v ->  es una lista con las condiciones de frontera en r=a
            a, b -> rmin, rmax
            npt -> son los números de puntos usados para generar h = (b-a)/(npt-1)
    
    Out:
    xval -> un arreglo discreto de valores de ri, donde se itero  
    yval -> una lista de arreglos discretos correspondientes a las iteraciones con la misma estructura que y0v. Cada fila corresponde a un ri.
    """
    y0v, a, b, npt = data

    # h val
    h = (b-a)/(npt-1)
    xval = a + np.arange(npt)*h
    nEq = len(y0v)
    yval = np.zeros((npt, nEq))  # crea una matriz (npt, nEq)
    
    yv = np.copy(y0v)
    for j, xi in enumerate(xval):
        yval[j, :] = yv
        k0 = h*f(xi, yv)  # importante como se define la función, tiene que ser de la forma f(x, yv) donde yv=[y0, y1, ...]
        k1 = h*f(xi+h/4, yv+k0/4)
        k2 = h*f(xi+3*h/8, yv+3*k0/32+9*k1/32)
        k3 = h*f(xi+12/13*h, yv+1932*k0/2197-7200*k1/2197+7296*k2/2197)
        k4 = h*f(xi+h, yv+439*k0/216-8*k1+3680*k2/513-845*k3/4104)
        k5 = h*f(xi+h/2, yv-8*k0/27+2*k1-3544*k2/2565+1859*k3/4104-11*k4/40)
        yv = yv + 16*k0/135 + 6656*k2/12825 + 28561*k3/56430-9*k4/50 + 2*k5/55
    return xval, yval