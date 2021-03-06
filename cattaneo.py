import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import ini
# therm, _, sources = ini.chelkanov2019()
therm, _, sources = ini.xiong2011()
rho = therm['rho']

k = therm['k']
C_V = therm['C_V']
tau = therm['tau']

grd = ini.Grid('x', 'y', 'z')
dx = grd.dx
dy = grd.dy
dz = grd.dz
dt = grd.dt

source_f, b_source_f = sources

class Model:
    def __init__(self, nx, ny, nz):
        self.cell = np.zeros((nx,ny, nz))
        self.fill_cells()
        self.mat_properties()
        self.fill_arrays()
        
        self.fill_cells()
        self.mat_properties()
        self.fill_arrays()
        
        
    def mat_properties(self):
        nx, ny, nz = self.cell.shape

        self.layers = np.unique(self.cell.flatten())
        self.lay_mask = np.empty((self.layers.shape[0], nx, ny, nz),dtype=np.bool)
        
        for i,lay in enumerate(self.layers):
            self.lay_mask[i,:,:,:] = (self.cell==lay)
        
    def fill_cells(self):
        nx, ny, nz = self.cell.shape
        self.cell[:,ny/3:ny/2,nz/2:] = 1
        self.cell[:,ny/3:ny/2,nz/2:] = 2
    def fill_arrays1(self):
        self.rho[i] = rho

        self.tau_tab[i] = tau
        self.tau[self.lay_mask[i, :, :, :]] = self.tau_tab[i]
        self.k_tab[i] = k
        self.k[self.lay_mask[i, :, :, :]] = self.k_tab[i]
        self.C_V_tab[i] = C_V
        self.C_V[self.lay_mask[i, :, :, :]] = self.C_V_tab[i]

    def fill_arrays(self):

        nx, ny, nz = self.cell.shape
        
        

        self.tau = np.zeros((nx,ny, nz))
        self.k = np.zeros((nx,ny, nz))
        self.C_V = np.zeros((nx,ny, nz))
        
        self.rho = [0]*self.layers.shape[0]
        self.tau_tab = [0]*self.layers.shape[0]
        self.k_tab = [0]*self.layers.shape[0]
        self.C_V_tab = [0]*self.layers.shape[0]
        for i,lay in enumerate(self.layers):
            self.rho[i] = rho

            self.tau_tab[i] = tau
            self.tau[self.lay_mask[i,:,:,:]] = self.tau_tab[i]
            self.k_tab[i] = k
            self.k[self.lay_mask[i,:,:,:]] = self.k_tab[i]
            self.C_V_tab[i] = C_V
            self.C_V[self.lay_mask[i,:,:,:]] = self.C_V_tab[i]
    def update(self, T):
        nx, ny, nz = self.cell.shape
        
        for i,lay in enumerate(self.layers):
            self.tau
            self.tau[self.lay_mask[i,:,:,:]] = data[i][1]
            self.k[self.lay_mask[i,:,:,:]] = data[i][2]
            self.C_V[self.lay_mask[i,:,:,:]] = data[i][3]
   
m = Model(grd.nx, grd.ny, grd.nz)

def show1():
    plt.figure()
    plt.subplot(1, 1, 1)
    plt.contourf([dy*i for i in range(grd.ny)], [dx*i for i in range(grd.nx)],m.tau[:,:,grd.nz/2])
    plt.gca().set_title("tau")
    plt.colorbar()
    plt.show()
#show1()

qx = np.zeros((grd.nx+1,grd.ny, grd.nz))#fluxes - faces
qy = np.zeros((grd.nx, grd.ny+1, grd.nz))
qz = np.zeros((grd.nx, grd.ny, grd.nz+1))

w = np.zeros((grd.nx, grd.ny, grd.nz))
source = np.zeros((grd.nx, grd.ny, grd.nz))
v_source_f = np.vectorize(source_f)

T = np.zeros((grd.nx, grd.ny, grd.nz))
T_old = np.zeros((grd.nx, grd.ny, grd.nz))
T[:, :, :] = 300
T_old[:, :, :] = 300
T_der_ = np.zeros((grd.nx, grd.ny, grd.nz, 6))
for ix in range(grd.nx):
    for iy in range(grd.ny):
        for iz in range(grd.nz):
            T_der_[ix, iy, iz, 0] = 0 - 3*dt
            T_der_[ix, iy, iz, 1] = 0 - 2*dt
            T_der_[ix, iy, iz, 2] = 0 - 1*dt
            T_der_[ix, iy, iz, 3:] = 300
T_der = np.zeros((grd.nx, grd.ny, grd.nz))
T0 = 300

w = T*m.C_V

k0 = 25

def flux_x(qx, T):
    #print np.any(m.tau.flatten() < 1e-14)
    #qx[1:-1,...] = ((dt)/(2*m.tau[1:,:,:]-dt))*qx[1:-1,...] -(k/dx)*((dt*(T[1:,...] - T[:-1,...]))/m.tau[1:,:,:])

    #p*q(n+1)+(1-p)*q(n)
    p = 1
    qx[1:-1, ...] = ((m.tau[1:, :, :] - dt*(1-p))/(p*dt + m.tau[1:, :, :]))*qx[1:-1, ...] \
                    - (k/dx)*((dt)/(dt*p+m.tau[1:, :, :]))*(T[1:, ...] - T[:-1, ...])


def flux_y(qy, T):
    qy[:, 1:-1, :] = qy[:, 1:-1, :]\
                     - (dt*qy[:, 1:-1, :])/m.tau[:, 1:, :] \
                     - (k/dy)*((dt*(T[:, 1:, :] - T[:, :-1, :]))/m.tau[:, 1:, :])


def flux_z(qz, T):
    qz[..., 1:-1] = qz[...,1:-1]-(dt*qz[...,1:-1])/m.tau[:,:,1:] -(k/dz)*((dt*(T[:,:,1:] - T[:,:,:-1]))/m.tau[:,:,1:])


def w_propgtn(qx, qy, qz, w, T,t):
    source[:, grd.ny/2, grd.nz/2] = v_source_f(source[:, grd.ny/2, grd.nz/2])
    source[:, :, :] = source[:, grd.ny / 2:grd.ny / 2+1, grd.nz / 2:grd.nz / 2+1]

    w[:, :, :] = w \
               - (dt/dx)*(qx[1:, :, :] - qx[:-1, :, :]) \
               - (dt/dy)*(qy[:, 1:, :] - qy[:, :-1, :]) \
               - (dt/dz)*(qz[:, :, 1:] - qz[:, :, :-1]) + dt*source
    # T_der[:,:,:] = (w/m.C_V-T)/dt
    tau0 = np.amax(m.tau[0, :, :])
    T[0, :, :] = b_source_f(t)

    w[0, :, :] = T[0, :, :]*m.C_V[0, ...]

    T[0:, :, :] = w[0:, :, :] / m.C_V

    T_der_[:, :, :, 0] = T_der_[:, :, :, 1]
    T_der_[:, :, :, 1] = T_der_[:, :, :, 2]
    T_der_[:, :, :, 2] = t

    T_der_[:, :, :, 3] = T_der_[:, :, :, 4]
    T_der_[:, :, :, 4] = T_der_[:, :, :, 5]
    T_der_[:, :, :, 5] = T


def boundaries(w, qx, qy, qz,T, t):
##    qx[0,:,:] = k0*(T[0,:,:] - 300)
##    qx[-1,:,:] = -k0*(400 - T[-1,:,:])
##    #qx[-1,...] = qx[-1,...]-(dt*qx[-1,...])/m.tau[-1,:,:] -(m.k[-1,...]/dx)*((dt*(400. - T[-1,...]))/m.tau[-1,...])
##
##    qy[:,0,:] = k0*(T[:,0,:] - 300)
##    qy[:,-1,:] = k0*(300 - T[:,-1,:])
##
##    qz[:,:,0] = k0*(T[:,:,0] - 300)
##    qz[:,:,-1] = k0*(300 - T[:,:,-1])
    qx[0, :, :] = 0
    qx[-1, :, :] = 0
    qy[:, 0, :] = 0
    qy[:, -1, :] = 0
    qz[:, :, 0] = 0
    qz[:, :, -1] = 0
    # T[0, :, :] =
    # w[0, :, :] = T[0, :, :]*m.C_V[0, :, :]
    

t = 0
def step(t):
    #print 't'
    flux_x(qx, T)
    flux_y(qy, T)
    flux_z(qz, T)
    w_propgtn(qx, qy, qz, w, T, t)
    boundaries(w, qx, qy, qz, T, t)

def show(step):
    plt.figure()
    plt.subplot(1, 1, 1)
    plt.contourf([dy*i for i in range(grd.ny)], [dx*i for i in range(grd.nx)],w[:,:,grd.nz/2])
    plt.gca().set_title(str(step))
    plt.show()


class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

def show1(step):
    fig,ax = plt.subplots()

    f1, f2 = 297., 303.
    cbarticks=np.arange(f1,f2,(f2-f1)/10., dtype = np.float)
    mesh = ax.pcolormesh([dy*i for i in range(grd.ny)], [dx*i for i in range(grd.nx)],T[:,:,grd.nz/2],cmap='coolwarm', vmin = f1, vmax = f2,
                         norm=MidpointNormalize(vmin = f1, vmax = f2, midpoint=300.))
    #print T[:,:,grd.nz/2].min(), T[:,:,grd.nz/2].max()                         
    cbar = fig.colorbar(mesh, ax=ax , cmap='coolwarm', ticks=cbarticks)
    cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.1lf'))
    cbar.ax.yaxis.set_minor_formatter(FormatStrFormatter('%.1lf'))
    cbar.ax.set_yticklabels(['{:.5f}'.format(x) for x in cbarticks])
    plt.show()

def show_tau():
    plt.figure()
    plt.subplot(1, 1, 1)
    plt.contourf([dy*i for i in range(grd.ny)], [dx*i for i in range(grd.nx)], m.tau[:,:,grd.nz/2],
                 vmin = 0.9*np.amin(m.tau), vmax = 1.1*np.amax(m.tau), extend='both')
    plt.gca().set_title('tau')
    plt.colorbar()
    #plt.show()

def show_cell():
    plt.figure()
    plt.subplot(1, 1, 1)
    plt.contourf([dy*i for i in range(grd.ny)], [dx*i for i in range(grd.nx)], m.cell[:,:,grd.nz/2],
                 vmin = 0.9*np.amin(m.cell), vmax = 1.1*np.amax(m.cell), extend='both')
    
    plt.gca().set_title('cell')
    plt.colorbar()
    #plt.show()

    
if __name__ == '__main__':
    for i in range(int(100*tau/dt)):
        
        step(t)
        t += dt
        if i == int(2*tau/dt) or i == int(10*tau/dt) :
            print i
            show1(t)


