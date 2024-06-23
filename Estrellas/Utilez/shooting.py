# Estrellas/Utilez

import numpy as np
from scipy.integrate import solve_ivp


def Freq_solve(system, in0, wmin, wmax, rmin=0, rmax=1e03, met='RK45', Rtol=1e-09, Atol=1e-10,
               limw=1e-14, alphB=1, info=False):
    """
    Implementación de un algoritmo de shooting usando un método de bisección para encontrar el valor de la 
    frecuencia w0, dado una amplitud central p0.

    In:
    system -> sitema de ecuaciones: eds.systemBS
    in0 -> condiciones iniciales del problema: [g0, n0, p0, p1] recordar que [g, g', N, N', Phi, Phi'] -> [g0, g1, n0, n1, p0, p1]
    wmin, wmax -> rango de valores en el que se buscará la frecuencia (recomendado añadir un check al código que verifique que wmin<wmax)
    rmax, rmin -> intervalo que se discretizará, por defecto rmin=0, rmax=1e03 (notar que lo q se busca es que alcance el límite en w0)
    met -> metodología usada por solve_ivp, las opciones son: 'DOP853', 'LSODA', 'RK45', por defecto está este último
    Rtol, Atol=1e-10 -> representan la tolerancia relativa y absoluta al usar solve_ivp, por defecto son Rtol=1e-09, Atol=1e-10
    limw -> representa la diferencia limítrofe que aceptaremos: abs((wmax-wmin)/2)<=limw Por defecto es limw=1e-14 
    alphB -> Valor de alphaB, por defecto alphB=1
    info -> Imprime información complementaria, por defecto info=False

    Out:
    w0, rTemp  -> valor de la frecuencia encontrada y el radio máximo de la iteración
    """


    # IMPORTANT: it is not possible to find two event at same time
    # Events
    def Sig(r, yV, arg): return yV[2]
    def dSig(r, yV, arg): return yV[3]

    # establece las direcciones
    Sig.direction = -1  
    dSig.direction = 1
    # If direction is positive, event will only trigger when going from negative to positive, 
    # and vice versa if direction is negative. If 0, then either direction will trigger event. Implicitly 0 if not assigned.
    
    # establece que acción tomar
    Sig.terminal = True
    dSig.terminal = True

    while True:
        w0 = (wmax+wmin)/2  # bisección metodología
        arg = [w0, alphB]

        sol = solve_ivp(system, [rmin, rmax], in0, events=(Sig, dSig),  # notar orden de los eventos
                         args=(arg,), method=met,  rtol=Rtol, atol=Atol)
        
        if sol.t[-1]==rmax:  # caso en que alcanzó el rmax
            print('Found', w0)
            rTemp = sol.t[-1]
            break  
        elif sol.t_events[0].size > 0:  # caso en que hay algún evento tipo Sig, es decir p0 negativa
            wmax = w0
            if info:
                print('Aumentando', 'w0 = ', w0, 'wmax = ', wmax, 'wmin ', wmin)
        elif sol.t_events[1].size > 0:  # caso en que hay algún evento tipo dSig, es decir p' positiva
            wmin = w0
            if info:
                print('Disminuyendo', 'w0 = ', w0, 'wmax = ', wmax, 'wmin ', wmin)

        # checking the lim freq.
        if abs((wmax-wmin)/2)<=limw:
            print('Maxima precisión alcanzada', w0, 'radio', sol.t[-1])
            rTemp = sol.t[-1]
            break

    return w0, rTemp


def Freq_solveNodos(system, in0, wmin, wmax, rmin=0, rmax=1e03, met='RK45', Rtol=1e-09, Atol=1e-10,
               limw=1e-14, alphB=1, nodos=0, info=False):
    """
    Implementación de un algoritmo de shooting usando un método de bisección para encontrar el valor de la 
    frecuencia w0, dado una amplitud central p0.

    In:
    system -> sitema de ecuaciones: eds.systemBS
    in0 -> condiciones iniciales del problema: [g0, n0, p0, p1] recordar que [g, g', N, N', Phi, Phi'] -> [g0, g1, n0, n1, p0, p1]
    wmin, wmax -> rango de valores en el que se buscará la frecuencia (recomendado añadir un check al código que verifique que wmin<wmax)
    rmax, rmin -> intervalo que se discretizará, por defecto rmin=0, rmax=1e03 (notar que lo q se busca es que alcance el límite en w0)
    met -> metodología usada por solve_ivp, las opciones son: 'DOP853', 'LSODA', 'RK45', por defecto está este último
    Rtol, Atol=1e-10 -> representan la tolerancia relativa y absoluta al usar solve_ivp, por defecto son Rtol=1e-09, Atol=1e-10
    limw -> representa la diferencia limítrofe que aceptaremos: abs((wmax-wmin)/2)<=limw Por defecto es limw=1e-14 
    alphB -> Valor de alphaB, por defecto alphB=1
    nodos -> los nodos de la configuracion
    info -> Imprime información complementaria, por defecto info=False

    Out:
    w0, rTemp, nodosPosit  -> valor de la frecuencia encontrada, el radio máximo de la iteración y la posición de los nodos
    """
    
    print('Finding a profile with ', nodos, 'nodes')
    
    def Sig(r, yV, arg): return yV[2]
    def dSig(r, yV, arg): return yV[3]

    # establece las direcciones
    Sig.direction = 0  # como pueden ser varios nodos se pone 0. Notar que no hay acción que tomar, solo almacenamos 
    dSig.direction = 0

    while True:
        w0 = (wmax+wmin)/2
        arg = [w0, alphB]

        sol = solve_ivp(system, [rmin, rmax], in0, events=(Sig, dSig),
                         args=(arg,), method=met,  rtol=Rtol, atol=Atol)
        
        if sol.t_events[1].size == nodos+1 and sol.t_events[0].size == nodos:
            print('Found', w0)
            rTemp = sol.t[-1]
            nodosPosit = sol.t_events[0]
            break
        elif sol.t_events[1].size > nodos+1:  # una vez por nodo
            if sol.t_events[0].size > nodos:  # dos veces por nodo
                wmax = w0
                rTemp = sol.t_events[0][-1]
            else:  # si pasa por cero más veces que 2*nodos se aumenta la w, sino se disminuye
                wmin = w0
                rTemp = sol.t_events[1][-1]
        elif sol.t_events[1].size <= nodos+1:
            if sol.t_events[0].size > nodos:  # dos veces por nodo
                wmax = w0
                rTemp = sol.t_events[0][-1]
            else:
                wmin = w0
                rTemp = sol.t_events[1][-1]

        # checking the lim freq.
        if abs((wmax-wmin)/2)<=limw:
            print('Maxima precisión alcanzada', w0, 'radio', rTemp)
            nodosPosit = sol.t_events[0]
            break
        
        if nodos==0:
            nodosPosit = None

    return w0, rTemp, nodosPosit 