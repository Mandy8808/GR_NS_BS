## Modulo v1.0
# 08/8/2024

import numpy as np
import time
import sys
import os
import zipfile
import glob

import matplotlib.pyplot as plt

########
# Main #
########
#   dt, C, m=0,  c=1, NoHomog, V=None, user_action=None, Vp=lambda u: 0

def solver(dt, dx, xmax, tmax, ux, dux=lambda x: 0, Vx=lambda u: 0, dVx=lambda u: 0, xmin=0, tmin=0, NoHomog=[0, 0],
           simulacion=True, address='Data', filename='.data', ZipName='SimulacionData', info=False, density=False,
           user_action=None):
    """
    Resolución de la ecuación diferencial con la estructura: 
        
        u_tt = u_xx - V'(u) 
        
    para una malla (xmin, xmax) \times (tmin, tmax), con dt, dx
    
    In:
        Condiciones a la frontera:
        NoHomog = [u(t, xmin), u(t, xmax)] : una lista con el valor en la frontera, predeterminado NoHomog = [0, 0] 
    
    Funciones complementarias:
        Condiciones iniciales:
        ux -> una función que corresponda a u(t=0, x)
        dux -> una función que corresponda a du(t=0, x)/dt, por defecto, dUx = 0

        Potencial
        Vx => una función que corresponde al potencial (se usa para calcular la densidad de energía)
        
        Derivada del potencial
        dVx -> una función que corresponde a la derivada con respecto al campo del potencial V, por defecto dVx=0
    
    Opciones:
        simulacion -> True para salvar los datos de las iteraciones numéricas, False para no salvarlo
        address -> Folder que se creará para salvar los datos de las iteraciones numéricas
        filename -> Prefijo de los archivos que se salvarán con las iteraciones numéricas
        ZipName -> Nombre del archivo general que salvará TODA la simulación
        info -> True para ver ciertas informaciones, False para no verlas
    
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
    funciones = [ux, dux, dVx]
    check(tmin, tmax, dt, xmin, xmax, dx, funciones)
    
    #### Malla
    ti, xi = malla(tmin, tmax, dt, xmin, xmax, dx)

    start = time.time()
    # Cargando dato inicial dentro de u_n
    u_n = ux(xi)
    
    if density:
        dens = densidad(u=u_n, xi=xi, dux=dux, V=Vx)
    else:
        dens = 0
    
    ########################################
    utilez = [u_n, ti, xi, -1, dens, data]  # u, ti, xi, n, data
    Operaciones(user_action, simulacion, utilez)
    ########################################

    # Primer paso de iteración temporal
    u = zerostep(xi, dx, dt, ux=u_n, dux=dux, dVx=dVx, NoHomog=NoHomog)
    
    if density:
        dens = densidad(u=u, xi=xi, um=u_n, dt=dt, V=Vx)
    else:
        dens = 0
        
    ########################################
    utilez = [u, ti, xi, 0, dens, data]  # u, ti, xi, n, data
    Operaciones(user_action, simulacion, utilez)
    ########################################

    # Intercambiando variables para el próximo paso de iteración
    u_nm1 = np.copy(u_n)
    u_n = np.copy(u)

    Nt = ti.size - 1
    for n in range(1, Nt):
        u = nsteps(dx, dt, ux=u_n, uxm1=u_nm1, dVx=dVx, NoHomog=NoHomog)

        if density:
            dens = densidad(u=u, xi=xi, um=u_n, dt=dt, V=Vx)
        else:
            dens = 0
        
        ########################################
        utilez = [u, ti, xi, n, dens, data]  # u, ti, xi, n, data
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

###############
# Iteraciones #
###############
def zerostep(xi, dx, dt, ux, dux=lambda x: 0, dVx=lambda u: 0, NoHomog=[0, 0]):
    """
    Primera iteración temporal
    
    In:
    xi -> Arreglo correspondiente a la discretización del espacio con dx
    dt -> ancho de la malla temporal
    ux -> Arreglo (de puntos) de la solución en la iteración t=0
    dux -> Función de la derivada de la solución en la iteración t=0
    dVx -> Función del potencial

    Out:
    u
    """
    C2 = (dt/dx)**2
    Nx = ux.size
    duxi = dux(xi[1:-1]) 
    
    u = np.zeros(Nx)
    u[1:-1] = ux[1:-1] + dt*duxi + 0.5*C2*(ux[2:] - 2*ux[1:-1] + ux[:-2]) - 0.5*dt**2*dVx(ux[1:-1])
    
    u[0] = NoHomog[0]
    u[-1] = NoHomog[1]
    return u


def nsteps(dx, dt, ux, uxm1, dVx=lambda u: 0, NoHomog=[0, 0]):
    """
    n-iteración temporal
    
    In:
    dx -> ancho de la malla espacial
    dt -> ancho de la malla temporal
    ux -> Arreglo (de puntos) de la solución en la iteración t=n
    uxm1 -> Arreglo correspondiente a la solución en dos niveles anteriores en el tiempo u^{n-1}_{i}
    dVx -> Función del potencial
    
    Out:
    u
    """
    C2 = (dt/dx)**2
    Nx = ux.size
    
    u = np.zeros(Nx)
    u[1:-1] = -uxm1[1:-1] + 2*ux[1:-1] + C2*(ux[2:] - 2*ux[1:-1] + ux[:-2]) - dt**2*dVx(ux[1:-1])
    
    u[0] = NoHomog[0]
    u[-1] = NoHomog[1]
    return u

def densidad(u, xi, dux=None, um=None, dt=None, V=lambda u:0):
    """
    rho = ((du/dt)^2 + (du/dx)^2)/2 + V(u)
    """ 
    if dux:
        rho = (dux(xi)**2 + np.gradient(u, xi)**2)/2 + V(u)
        #plt.plot(xi, rho)
        #plt.show()
    else:
        du = (u-um)/dt
        rho = (du**2 + np.gradient(u, xi)**2)/2 + V(u)
    return rho

#########################
# Funciones secundarias #
#########################
def check(tmin, tmax, dt, xmin, xmax, dx, funciones) :
    """ 
    Revisando si el dx, y el dt es el adecuado
    Revisando si las funciones auxiliares están bien definidas
    """
    # check dt
    if dt>(tmax-tmin):
        sys.exit("El paso temporal dt debe ser mas pequeño que el intervalo temporal tmax-tmin")
    if dx>(xmax-xmin):
        sys.exit("El paso espacial dx debe ser mas pequeño que el intervalo espacial xmax-xmin")
    
    # check si es funcion
    text = ['Error: u(t=0, x) debe ser una función',
            'Error: du(t=0, x)/dt debe ser una función',
            'Error dV/du debe ser una función'
            ]
    for i in range(len(funciones)):
        if not callable(funciones[i]):
            sys.exit(text[i])

def malla(tmin, tmax, dt, xmin, xmax, dx):
    """ 
    Devuelve la malla espacial para un dt y un intervalo tmin-tmax definido, asumiendo una discretizacion espacial dx entre xmin-xmax
    """
    Nt = int(round((tmax-tmin)/dt))
    ti = tmin + dt*np.arange(Nt+1)
    
    Nx = int(round((xmax-xmin)/dx))
    xi = xmin + dx*np.arange(Nx+1)
    return ti, xi

def Operaciones(user_action, simulacion, utilez):
    u, ti, xi, n, dens, data = utilez
    if user_action:
        user_action(u, xi, ti, n)
    
    if simulacion:
        data.save_file(u, dens, ti, xi, n)
        

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
    def save_file(self, dataU, dataDens, ti, xi, n):

        # verificando si existe o creando la carpeta
        if not os.path.exists(self.direccion):
            os.makedirs(self.direccion)

        # salvando u
        kwargs = {self.nombre + '%d'%(n+1): dataU} 
        fname = self.direccion + '/' + self.nombre + '_u_%d'%(n+1) + '.dat'
        np.savez(fname, **kwargs)
        
        # salvando densidad
        kwargs = {self.nombre + '_dens_' + '%d'%(n+1): dataDens} 
        fname = self.direccion + '/' + self.nombre + '_dens_%d'%(n+1) + '.dat'
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

    def __init__(self, address, figData, solAnalit=None, info=False):
        if info:
            print(f"Usando dirección {address}")
        
        # Atributos de instancia
        self.direccion = address
        self.figData = figData
        self.solAnalit = solAnalit

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
        
        if self.solAnalit:
            vt0 = np.ones(len(xi))*t0
            frame1, = ax.plot(xi, self.solAnalit(vt0, xi), ls=':', color='k', lw=lw)
        else:
            frame1 = None
            
        tframe = ax.text(xlim[-1]-xlim[-1]/2, max(u0)-max(u0)/8, s=r'time=$%3.2f$'%t0, fontsize='small')

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

        return fig, ax, frame, frameSol, tframe, frame1
    
    def updateframe(self, ind, frame, frameSol, tframe, frame1,
                    ax, ti, dataname, array_names, solEx, ylim):
        """ 
        """
        t = ti[ind]
        #print(t)
        tframe.set_text(r'time=$%4.3f$'%t)

        ui = array_names[dataname+'%d'%ind]

        # updating axis
        frame.set_ydata(ui)
        
        if frame1:
            xi = frame1.get_xdata()
            vt = np.ones(len(xi))*t
            frame1.set_ydata(self.solAnalit(vt, xi))
        
        if frameSol:
            xi = frameSol.get_xdata()
            uiE = solEx(xi, t)
            frameSol.set_ydata(uiE)

        # updating y-lim
        if ylim:
            ax.set_ylim(ylim[0], ylim[1])
        else:    
            ax.set_ylim(min(ui)-min(ui)/8, max(ui)+max(ui)/8)

        return frame, frameSol, tframe, frame1

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
        
        fig, ax, frame, frameSol, tframe, frame1 = self.frame0(data, struc, xlim=xlim, ylim=ylim, solEx=solEx, show=show)
        frames = len(ti[n0:])

        anim = animation.FuncAnimation(fig, self.updateframe,
                                       frames=range(n0, frames),
                                       fargs=(frame, frameSol, tframe, frame1, ax, ti, dataname, array_names, solEx, ylim),
                                       interval=50, blit=False)
        
        fps, bitrate, extra_args = vconf
        anim.save(nameV+'.mp4', bitrate=bitrate, fps=fps, extra_args=extra_args)  # direc+nameV+

