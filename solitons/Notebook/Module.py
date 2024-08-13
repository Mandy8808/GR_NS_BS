## Modulo v1.0
# 08/8/2024

import numpy as np
import time
import sys
import os
import zipfile
import glob

########
# Main #
########

def solver(I, xmax, tmax, dt, C, t0=0, x0=0, c=1, f=None, V=None, user_action=None,
           simulacion=True, address='Data', filename='.data', ZipName='SimulacionData',
           info=False):
    """
    Resolución de la ecuación diferencial: u_tt = c^2*u_xx + f para una malla (0, L) \times (0,T]

    Variables:
    xi -> Arreglo correspondiente a la discretización del espacio
    ti -> Arreglo correspondiente a la discretización del tiempo
    u -> Arreglo (de puntos espaciales) correspondiente a la solución en la próxima iteración en el tiempo u^{n+1}_{i}
    u_n -> Arreglo correspondiente a la solución en el nivel anterior en el tiempo u^{n}_{i}
    u_nm1 -> Arreglo correspondiente a la solución en dos niveles anteriores en el tiempo u^{n-1}_{i}

    In:

    Out:
    
    """
    if simulacion:
        data = StoreSolution(address, filename, info=info)
    else:
        data = None

    #### Checks
    f, V = funciones(f, V)
    check(dt, tmax, f, V, I)
    
    #### Malla
    ti, xi = malla(t0, tmax, dt, x0, xmax, c, C)

    start = time.time()
    # Cargando dato inicial dentro de u_n
    u_n = I(xi)
    
    ########################################
    utilez = [u_n, ti, xi, -1, data]  # u, ti, xi, n, data
    Operaciones(user_action, simulacion, utilez)
    ########################################

    # Primer paso de iteración temporal
    utilez = [dt, f, V, C]
    u = zerostep(xi, ti, u_n, utilez)

    ########################################
    utilez = [u, ti, xi, 0, data]  # u, ti, xi, n, data
    Operaciones(user_action, simulacion, utilez)
    ########################################

    # Intercambiando variables para el próximo paso de iteración
    u_nm1 = np.copy(u_n)
    u_n = np.copy(u)

    Nt = ti.size - 1
    for n in range(1, Nt):
        utilez = [dt, f, C]
        u = nsteps(n, xi, ti, u_n, u_nm1, utilez)

        ########################################
        utilez = [u, ti, xi, n, data]  # u, ti, xi, n, data
        Operaciones(user_action, simulacion, utilez)
        ########################################

        # Intercambiando variables para el próximo paso de iteración
        u_nm1 = np.copy(u_n)
        u_n = np.copy(u)

    if simulacion:
        data.close_file(ZipName)
    
    end = time.time()
    cpu_time = end - start
    
    if info:
        print('cpu time:', cpu_time)

    return u, xi, ti, cpu_time


def solverO(I, xmax, tmax, dt, C, t0=0, x0=0, c=1, f=None, V=None, user_action=None, info=False):
    """
    Resolución de la ecuación diferencial: u_tt = c^2*u_xx + f para una malla (0, L) \times (0,T]

    Variables:
    xi -> Arreglo correspondiente a la discretización del espacio
    ti -> Arreglo correspondiente a la discretización del tiempo
    u -> Arreglo (de puntos espaciales) correspondiente a la solución en la próxima iteración en el tiempo u^{n+1}_{i}
    u_n -> Arreglo correspondiente a la solución en el nivel anterior en el tiempo u^{n}_{i}
    u_nm1 -> Arreglo correspondiente a la solución en dos niveles anteriores en el tiempo u^{n-1}_{i}

    In:

    Out:
    
    """

    f, V = funciones(f, V)

    #### Checks
    check(dt, tmax, f, V, I)
    
    #### Malla
    ti, xi = malla(t0, tmax, dt, x0, xmax, c, C)

    start = time.time()
    # Cargando dato inicial dentro de u_n
    u_n = I(xi) 

    # Primer paso de iteración temporal
    utilez = [dt, f, V, C]
    u = zerostep(xi, ti, u_n, utilez)

    if user_action:
        user_action(u, xi, ti, n=0)

    # Intercambiando variables para el próximo paso de iteración
    u_nm1 = np.copy(u_n)
    u_n = np.copy(u)

    Nt = ti.size - 1
    for n in range(1, Nt):
        utilez = [dt, f, C]
        u = nsteps(n, xi, ti, u_n, u_nm1, utilez)

        if user_action:
            if user_action(u, xi, ti, n):
                break
        
        # Intercambiando variables para el próximo paso de iteración
        u_nm1 = np.copy(u_n)
        u_n = np.copy(u)

    end = time.time()
    cpu_time = end - start
    if info:
        print('cpu time:', cpu_time)

    return u, xi, ti, cpu_time

###############
# Iteraciones #
###############
def zerostep(xi, ti, u_n, utilez):
    """
    Primera iteración temporal
    
    In:
    xi -> Arreglo correspondiente a la discretización del espacio
    ti -> Arreglo correspondiente a la discretización del tiempo
    u -> Arreglo (de puntos) de la solución en la próxima iteración en el tiempo
    u_n -> Arreglo de la solución en el nivel anterior en el tiempo  Solution at 1 time level back

    Out:
    u
    """

    dt, f, V, C = utilez
    C2 = C**2
    
    Vi = V(xi[1:-1])
    n = 0
    fi = f(xi[1:-1], ti[n])
    
    Nx = u_n.size - 1
    u = np.zeros(Nx+1)
    u[1:-1] = u_n[1:-1] + dt*Vi + 0.5*C2*(u_n[2:] - 2*u_n[1:-1] + u_n[:-2]) + 0.5*dt**2*fi  # revisar signo + dt*Vi 
    
    return u

def nsteps(n, xi, ti, u_n, u_nm1, utilez):
    """
    Primera iteración temporal
    
    In:
    xi -> Arreglo correspondiente a la discretización del espacio
    ti -> Arreglo correspondiente a la discretización del tiempo
    u -> Arreglo (de puntos espaciales) correspondiente a la solución en la próxima iteración en el tiempo u^{n+1}_{i}
    u_n -> Arreglo correspondiente a la solución en el nivel anterior en el tiempo u^{n}_{i}
    u_nm1 -> Arreglo correspondiente a la solución en dos niveles anteriores en el tiempo u^{n-1}_{i}
    
    Out:
    u
    """
    dt, f, C = utilez
    C2 = C**2
    
    fi = f(xi[1:-1], ti[n]*np.ones(len(xi[1:-1])))
    
    Nx = u_n.size
    u = np.zeros(Nx)
    u[1:-1] = -u_nm1[1:-1] + 2*u_n[1:-1] + C2*(u_n[2:] - 2*u_n[1:-1] + u_n[:-2]) + dt**2*fi
    return u


#########################
# Funciones secundarias #
#########################
def funciones(f, V):
    """ 
    """
    if (f is None) or (f==0):
        f=lambda x, t: 0
    
    if (V is None) or (V==0):
        V=lambda x: 0
    return f, V

def check(dt, tmax, f, V, I):
    """ 
    Función que checa que las entradas f, V sean objetos tipo funciones y el paso temporal dado sea menor a tmax
    """
    # check dt
    if dt>tmax:
        sys.exit("El paso temporal dt debe ser mas pequeño que el intervalo temporal tmax")
    # check si es funcion
    if not callable(f):
        sys.exit("La fuente f debe ser una función")
    if not callable(V):
        sys.exit("El potencial V debe ser una función")
    if not callable(I):
        sys.exit("La configuración inicial debe ser una función")

def malla(t0, tmax, dt, x0, xmax, c, C):
    """ 
    Devuelve la malla espacial para un dt y tmax definido, asumiendo una discretizacion espacial de la forma:
    dx = c dt/C
    """
    Nt = int(round(tmax/dt))
    ti = t0 + dt*np.arange(Nt+1)
    dx = dt*c/float(C)  # 
    Nx = int(round((xmax-x0)/dx))
    xi = x0 + dx*np.arange(Nx+1)
    return ti, xi

def Operaciones(user_action, simulacion, utilez):
    u, ti, xi, n, data = utilez
    if user_action:
        user_action(u, xi, ti, n)
    
    if simulacion:
        data.save_file(u, ti, xi, n)
        

def JoinFilesInOneZip(archives, archive_name):
    """
    Uniendo todos los archivos en un .zip

    IN:
    archive_name -> nombre que daremos al .zip
    archives -> lista o tupla de datos, o dirección donde encontraremos los datos

    Out:
    Crea un archivo: archive_name.zip
    """

    # creando el .zip 
    archive = zipfile.ZipFile(archive_name, mode='w', compression=zipfile.ZIP_DEFLATED, allowZip64=True)
    
    # identificando los archivos a comprimir
    if isinstance(archives, (list, tuple)):
        filenames = archives
    elif isinstance(archives, str):
        filenames = glob.glob(archives)

    # Open each archive and write to the common archive
    for filename in filenames:
        f = zipfile.ZipFile(filename, mode='r', compression=zipfile.ZIP_DEFLATED)
        for name in f.namelist():
            data = f.open(name, 'r')
            # Save under name without .npy
            archive.writestr(name[:-4], data.read())
        f.close()
        os.remove(filename)
    archive.close()

class StoreSolution:
    """ 
    Class para salvar soluciones
    """

    # El método __init__ es llamado al crear el objeto
    def __init__(self, address, filename, info=False):
        if info:
            print(f"Usando dirección {address} y nombre {filename}")
        
        # Atributos de instancia
        self.direccion = address
        self.nombre = filename
        self.time = []

    # METODOS
    def save_file(self, dataU, ti, xi, n):

        # verificando si existe o creando la carpeta
        if not os.path.exists(self.direccion):
            os.makedirs(self.direccion)

        # salvando u
        kwargs = {self.nombre + '%d'%(n+1): dataU} 
        fname = self.direccion + '/' + self.nombre + '_u_%d'%(n+1) + '.dat'
        np.savez(fname, **kwargs)
        
        # almacenando el tiempo
        self.time.append(ti[n+1])
        
        # salvando la malla
        if n==-1:
            np.savez(self.direccion + '/' + self.nombre + '_x.dat', x=xi)

    def close_file(self, ZipName):
        """
        Cada np.savez crea un nuevo archivo (zip archivo) con extension .npz (oculto)
        Por conveniencia es mejor unir todos en un solo archivo (zip) con extension .npz
        """
        if self.direccion is not None:
            # Save all the time points where solutions are saved
            np.savez(self.direccion + '/' + self.nombre + '_t.dat', t=np.array(self.time, dtype=float))
            
            # Uniendo todos los datos en un .zip
            archive_name = self.direccion + '/' + ZipName + '.npz'
            #filenames = glob.glob(self.direccion + '/' + self.nombre + '*.dat.npz')
            filenames = self.direccion + '/' + self.nombre + '*.dat.npz'
            JoinFilesInOneZip(filenames, archive_name)



class Visualization:

    def __init__(self, address, figData, info=False):
        if info:
            print(f"Usando dirección {address}")
        
        # Atributos de instancia
        self.direccion = address
        self.figData = figData

    # METODOS
    def loadData(self):
        array_names = np.load(self.direccion)
        return array_names
    
    def frame0(self, data, struc, xlim=None, ylim=None, solEx=None, show=False):
        """ 
        Primer frame
        """
        u0, xi, t0 = data
        xmin, xmax, ls, lw, color = struc
        fig, ax = self.figData

        frame, = ax.plot(xi, u0, ls=ls, c=color, lw=lw)
        tframe = ax.text(xmax-xmax/4, max(u0)-max(u0)/8, s=r'time=$%3.2f$'%t0, fontsize='small')

        if solEx:
            xi = frame.get_xdata()
            xval = np.linspace(xi[0], xi[-1], 1000)
            yi = solEx(xval, t0)
            frameSol, = ax.plot(xval, yi, ls='--', c='k', lw=1)
        else:
            frameSol = None

        if xlim:
            ax.set_xlim(xlim[0], xlim[1])
        else:
            ax.set_xlim(xmin, xmax)
        
        if ylim:
            ax.set_ylim(ylim[0], ylim[1])
        else:    
            ax.set_ylim(min(u0)-min(u0)/8, max(u0)+max(u0)/8)
        
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$u$')

        if show:
            import matplotlib.pyplot as plt
            plt.show()

        return fig, ax, frame, frameSol, tframe
    
    def updateframe(self, ind, frame, frameSol, tframe, 
                    ax, ti, dataname, array_names, solEx, ylim):
        """ 
        """
        t = ti[ind]
        #print(t)
        tframe.set_text(r'time=$%7.6f$'%t)

        ui = array_names[dataname+'%d'%ind]

        # updating axis
        frame.set_ydata(ui)
        if frameSol:
            xi = frameSol.get_xdata()
            uiE = solEx(xi, t)
            frameSol.set_ydata(uiE)

        # updating y-lim
        if ylim:
            ax.set_ylim(ylim[0], ylim[1])
        else:    
            ax.set_ylim(min(ui)-min(ui)/8, max(ui)+max(ui)/8)

        return frame, frameSol, tframe,

    def video(self, n0, dataname, struc, nameV, direc='/', xlim=None, ylim=None,
              vconf=[10, 1000000, ['-vcodec', 'libx264']],
              solEx=None, show=True):
        """ 
        video
        """
        import matplotlib.animation as animation

        array_names = self.loadData()  # cargando datos
        
        u0 = array_names[dataname+'%d'%n0]
        ti = array_names['t']
        xi = array_names['x']
        
        data = [u0, xi, ti[n0]]
        
        fig, ax, frame, frameSol, tframe = self.frame0(data, struc, xlim=xlim, ylim=ylim, solEx=solEx, show=show)
        frames = len(ti[n0:])

        anim = animation.FuncAnimation(fig, self.updateframe,
                                       frames=range(n0, frames),
                                       fargs=(frame, frameSol, tframe, ax, ti, dataname, array_names, solEx, ylim),
                                       interval=50, blit=False)
        
        fps, bitrate, extra_args = vconf
        anim.save(nameV+'.mp4', bitrate=bitrate, fps=fps, extra_args=extra_args)  # direc+nameV+

