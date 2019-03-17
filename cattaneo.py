import numpy as np
import matplotlib.pyplot as plt

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
<<<<<<< HEAD
        
        self.layers = np.unique(self.cell.flatten())
        self.lay_mask = np.empty((self.layers.shape[0], nx, ny, nz),dtype=np.int32)
        
        for i,lay in enumerate(self.layers):
            self.lay_mask[i] = (self.cell==lay)
=======
        self.layers = np.unique(self.cell.flatten())
        self.lay_mask = np.empty((self.layers.shape[0], nx, ny, nz),dtype=np.bool)
        
        for i,lay in enumerate(self.layers):
            self.lay_mask[i,:,:,:] = (self.cell==lay)
>>>>>>> 4fe6b7a84335f1146885b1b3ac9992062f9ca109
        
    def fill_cells(self):
        nx, ny, nz = self.cell.shape
        
        self.cell[nx/2:,ny/3:ny/2,:] = 1
        self.cell[:nx/2,ny/3:ny/2,:] = 2

    
    def fill_arrays(self):
<<<<<<< HEAD
        nx = self.cell.shape[0]
        ny = self.cell.shape[1]
        nz = self.cell.shape[2]
        self.C_V = np.zeros((nx,ny, nz))
        self.tau = np.zeros((nx,ny, nz))
        self.k = np.zeros((nx,ny, nz))
        
        data = np.loadtxt('data')
        print self.layers.shape[0]
        self.rho = np.zeros((self.layers.shape[0]))
        print self.rho
        for i,lay in enumerate(self.layers):
            print data[i,:]
            print self.lay_mask[i,:,:,nz/2]
            self.rho[i] = data[i,0]
            self.tau[self.lay_mask[i]] = data[i,1]
            self.k[self.lay_mask[i]] = data[i,2]
            self.C_V[self.lay_mask[i]] = data[i,3]
=======
        nx, ny, nz = self.cell.shape
        
        data = np.loadtxt('C:\ks_work\calculations\cattaneo_py_test\data1')
>>>>>>> 4fe6b7a84335f1146885b1b3ac9992062f9ca109

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
    def update(self, T):
        nx, ny, nz = self.cell.shape
        
        for i,lay in enumerate(self.layers):
            self.tau
            self.tau[self.lay_mask[i,:,:,:]] = data[i][1]
            self.k[self.lay_mask[i,:,:,:]] = data[i][2]
            self.C_V[self.lay_mask[i,:,:,:]] = data[i][3]
   
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
dt = 1e-8
dx = 1e-3
dy = 1.5e-3
dz = 1.2e-3
a = dx
def show1():
    plt.figure()
    plt.subplot(1, 1, 1)
    plt.contourf([dy*i for i in range(grd.ny)], [dx*i for i in range(grd.nx)],m.tau[:,:,grd.nz/2])
    plt.gca().set_title("tau")
    plt.colorbar()
    plt.show()
show1()

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



w = T*m.C_V

def flux_x(qx, T):
    #print np.any(m.tau.flatten() < 1e-14)
    qx[1:-1,...] = qx[1:-1,...]-(dt*qx[1:-1,...])/m.tau[1:,:,:] -(k/dx)*((dt*(T[1:,...] - T[:-1,...]))/m.tau[1:,:,:])

def flux_y(qy, T):
    qy[:,1:-1,:] = qy[:,1:-1,:]-(dt*qy[:,1:-1,:])/m.tau[:,1:,:] -(k/dy)*((dt*(T[:,1:,:] - T[:,:-1,:]))/m.tau[:,1:,:])

def flux_z(qz, T):
    qz[...,1:-1] = qz[...,1:-1]-(dt*qz[...,1:-1])/m.tau[:,:,1:] -(k/dz)*((dt*(T[:,:,1:] - T[:,:,:-1]))/m.tau[:,:,1:])

def w_propgtn(qx, qy, qz, w, T):
    w[:,:,:] = w + (dt/dx)*(qx[1:,:,:] - qx[:-1,:,:])+ (dt/dy)*(qy[:,1:,:] - qy[:,:-1,:])+ (dt/dz)*(qz[:,:,1:] - qz[:,:,:-1])
    T[:,:,:] = w/m.C_V
    if np.any(T.flatten() > 600) or np.any(T.flatten() < 250):
        print 'ALARM'
        print T[:,:,grd.nz/2] > 600
        raw_input()

def boundaries(w, qx, qy, qz,T):
    qx[0,:,:] = k0*(T[0,:,:] - 400)
    qx[-1,:,:] = k0*(300 - T[-1,:,:])

    qy[:,0,:] = k0*(T[:,0,:] - 300)
    qy[:,-1,:] = k0*(300 - T[:,-1,:])

    qz[:,:,0] = k0*(T[:,:,0] - 300)
    qz[:,:,-1] = k0*(300 - T[:,:,-1])

t = 0
def step(t):
    flux_x(qx, T)
    flux_y(qy, T)
    flux_z(qz, T)
    w_propgtn(qx, qy, qz, w, T)
    boundaries(w, qx, qy, qz, T)

def show(step):
    plt.figure()
    plt.subplot(1, 1, 1)
    plt.contourf([dy*i for i in range(grd.ny)], [dx*i for i in range(grd.nx)],w[:,:,grd.nz/2])
    plt.gca().set_title(str(step))
    plt.show()


<<<<<<< HEAD


show1()
for i in range(20000):
    
    if i%1000 == 0:
        print i
        #show(i)
    step(t)
    t += dt
=======
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

    f1, f2 = 295., 305.
    cbarticks=np.arange(f1,f2,(f2-f1)/10., dtype = np.float)
    mesh = ax.pcolormesh([dy*i for i in range(grd.ny)], [dx*i for i in range(grd.nx)],T[:,:,grd.nz/2],cmap='coolwarm', vmin = f1, vmax = f2,
                         norm=MidpointNormalize(vmin = f1, vmax = f2, midpoint=300.))
    print T[:,:,grd.nz/2].min(), T[:,:,grd.nz/2].max()                         
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
    for i in range(100):
        print i
        step(t)
        t += dt
        if i%20 == 0:
            show1(t)
>>>>>>> 4fe6b7a84335f1146885b1b3ac9992062f9ca109

