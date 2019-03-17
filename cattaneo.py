import numpy as np
import matplotlib.pyplot as plt

class Model:
    def __init__(self, nx, ny, nz):
        self.cell = np.zeros((nx,ny, nz))
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
        self.cell[:,ny/3:ny/2,nz/2:] = 1
        self.cell[:,ny/3:ny/2,nz/2:] = 2

    
    def fill_arrays(self):
        nx, ny, nz = self.cell.shape
        
        data = np.loadtxt('C:\ks_work\calculations\cattaneo_py_test\data1')

        self.tau = np.zeros((nx,ny, nz))
        self.k = np.zeros((nx,ny, nz))
        self.C_V = np.zeros((nx,ny, nz))
        
        self.rho = [0]*self.layers.shape[0]
        self.tau_tab = [0]*self.layers.shape[0]
        self.k_tab = [0]*self.layers.shape[0]
        self.C_V_tab = [0]*self.layers.shape[0]
        for i,lay in enumerate(self.layers):
            self.rho[i] = data[i][0]

            self.tau_tab[i] = data[i][1]
            self.tau[self.lay_mask[i,:,:,:]] = self.tau_tab[i]
            self.k_tab[i] = data[i][2]
            self.k[self.lay_mask[i,:,:,:]] = self.k_tab[i]
            self.C_V_tab[i] = data[i][3]
            self.C_V[self.lay_mask[i,:,:,:]] = self.C_V_tab[i]
   
class Grid:
    def __init__(self, x, y):
        self.nx = 25
        self.ny = 20
        self.nz = 23
        self.xi = [float(i+1) for i in range(self.nx)]
        self.yi = [0.1*i for i in range(self.ny)]
        #self.zi = [0.1*i for i in range(self.nz)]
      
        self.dxi05 = [i-j for i,j in zip(self.xi[1:], self.xi)]
        self.dxi05.insert(0, self.dxi05[0])
        self.dxi05.append(self.dxi05[-1])
        self.xi05 = [i-0.5*j for i,j in zip(self.xi, self.dxi05)]
        self.xi05.append(self.xi[-1]+0.5*self.dxi05[-1])
      
#read .grd file
grd = Grid('x','y')
m = Model(grd.nx, grd.ny, grd.nz)

qx = np.zeros((grd.nx+1,grd.ny, grd.nz))#fluxes - faces
qy = np.zeros((grd.nx,grd.ny+1, grd.nz))
qz = np.zeros((grd.nx,grd.ny, grd.nz+1))

w = np.zeros((grd.nx,grd.ny, grd.nz))#w
T = np.zeros((grd.nx,grd.ny, grd.nz))#w
T[:,:,:] = 300


k0 = 25 #Wt/m
rho = 2.7

tau = 1e-9
k = 1e6
C_V = 1e5

dt = 1e-8
dx = 1e-3
dy = 1.5e-3
dz = 1.2e-3
a = dx

w = T*C_V
#w[:,:,:] = 300*C_V

def flux_x(qx, w):
    qx[1:-1,...] = (1-dt)*qx[1:-1,...] -(k*dt/(C_V*dx))*(w[1:,...] - w[:-1,...])

def flux_y(qy, w):
    qy[:,1:-1,:] = (1-dt)*qy[:,1:-1,:] -(k*dt/(C_V*dy))*(w[:,1:,:] - w[:,:-1,:])

def flux_z(qz, w):
    qz[...,1:-1] = (1-dt)*qz[...,1:-1] -(k*dt/(C_V*dz))*(w[:,:,1:] - w[:,:,:-1])

def w_propgtn(qx, qy, qz, w):
    w[:,:,:] = w[:,:,:] + (dt/dx)*(qx[1:,:,:] - qx[:-1,:,:])+ (dt/dy)*(qy[:,1:,:] - qy[:,:-1,:])+ (dt/dz)*(qz[:,:,1:] - qz[:,:,:-1])

def boundaries(w, qx, qy, qz):
    qx[0,:,:] = k0*(w[0,:,:]/C_V - 400)
    qx[-1,:,:] = k0*(300 - w[-1,:,:]/C_V)

    qy[:,0,:] = k0*(w[:,0,:]/C_V - 300)
    qy[:,-1,:] = k0*(300 - w[:,-1,:]/C_V)

    qz[:,:,0] = k0*(w[:,:,0]/C_V - 300)
    qz[:,:,-1] = k0*(300 - w[:,:,-1]/C_V)

t = 0
def step(t):
    flux_x(qx, w)
    flux_y(qy, w)
    flux_z(qz, w)
    w_propgtn(qx, qy, qz, w)
    boundaries(w, qx, qy, qz)

def show(step):
    plt.figure()
    plt.subplot(1, 1, 1)
    plt.contourf([dy*i for i in range(grd.ny)], [dx*i for i in range(grd.nx)],w[:,:,grd.nz/2])
    plt.gca().set_title(str(step))
    plt.show()

    
for i in range(20000):
    if i%1000 == 0:
        print i
        show(i)
    step(t)
    t += dt

