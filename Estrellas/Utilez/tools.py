
# Estrellas/Utilez

def progressbar(current_value, total_value, bar_lengh, progress_char): 
    """
    Barra de progreso
    """
    percentage = int((current_value/total_value)*100)                                                # Percent Completed Calculation 
    progress = int((bar_lengh * current_value ) / total_value)                                       # Progress Done Calculation 
    loadbar = "Progress: [{:{len}}]{}%".format(progress*progress_char,percentage, len = bar_lengh)    # Progress Bar String
    print(loadbar, end='\r')

def find_nearest(array, value):
    """
    Encontrando el valor más cercano
    """
    n = [abs(i-value) for i in array]
    idx = n.index(min(n))
    #print(idx)
    return (array[idx], idx)

def escalamiento(Nmax, gmax):
    """ 
    Si queremos que asintóticamente la configuración sea Sch, se ha de cumplir que para r->inf:

    g_tt*g_rr = 1

    Ahora como tenemos una libertad de escalamiento de g_tt, entonces hemos de encontrar el valor de la constante 'c'
    que garantice la igualdad anterior, es decir:
    
    c*g_tt(rMax)*g_rr(rMax) = 1 

    Una vez encontrada SOLO debemos escalar el perfil métrico como 
    
    g^(sch)_tt = c*g_tt

    y la frecuencia como

    w^(phy) = c*w

    IMPORTANTE: pueden comprobar que si tomo como condición inicial para N0 el valor g^(sch)_tt (r=0), y para w el valor w^(phy) el perfil que se obtiene 
    satisface directamente (sin necesidad de escaleo) el "pegado con Sch".
    """

    c = 1/(Nmax*gmax)
    return c