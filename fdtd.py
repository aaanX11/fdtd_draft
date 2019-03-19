import numpy as np
import matplotlib.pyplot as plt

class Space:
    def __init__(self, nx, ny):
        self.cell = np.zeros((nx,ny))

class Grid:
    def __init__(self, nx, ny):
        self.nx = nx
        self.ny = ny
        #self.nz = 10
        self.xi = [float(i+1) for i in range(self.nx)]
        self.yi = [0.1*i for i in range(self.ny)]
        #self.zi = [0.1*i for i in range(self.nz)]
      
        self.dxi05 = [i-j for i,j in zip(self.xi[1:], self.xi)]
        self.dxi05.insert(0, self.dxi05[0])
        self.dxi05.append(self.dxi05[-1])
        self.xi05 = [i-0.5*j for i,j in zip(self.xi, self.dxi05)]
        self.xi05.append(self.xi[-1]+0.5*self.dxi05[-1])
      
#read .grd file
import cattaneo as therm
grd = Grid(therm.grd.nx, therm.grd.ny)
s = Space(grd.nx, grd.ny)

vx = np.zeros((grd.nx+1,grd.ny))#velocities - edges
vy = np.zeros((grd.nx,grd.ny+1))

sigxx = np.zeros((grd.nx,grd.ny))#p, lambda, mu - centers
sigyy = np.zeros((grd.nx,grd.ny))

sigxy = np.zeros((grd.nx+1,grd.ny+1))#shear stress, density - corners

rho = 1740
llmbda = 9.4e10
mu = 4.0e10
dt = 1e-10
dx = 6e-8
dy = 6e-8
a = dx

def get_source():
    import math
    vec = np.asarray([dy*(i-grd.ny/2) for i in range(grd.ny)])
    A = 2/np.sqrt(3*a)* math.pow(math.pi,1./3)
    res = A*(1 - vec**2/a**2)*np.exp(-vec**2/a**2)
    plt.plot( res)
    plt.show()
    return res#A*(1 - vec**2/a**2) 
    
source = get_source()

def velocity(vx, vy, sigxx, sigyy, sigxy):
    #nx+1 ny
    #inner edges X cells
    vx[1:-1,:] = vx[1:-1,:] + (dt/(rho*dx))*(sigxx[1:,:]-sigxx[:-1,:]) + (dt/(rho*dy))*(sigxy[1:-1,1:]-sigxy[1:-1,:-1])
    #nx ny+1
    #cells X inner edges
    vy[:,1:-1] = vy[:,1:-1] + (dt/(rho*dy))*(sigyy[:,1:]-sigyy[:,:-1]) + (dt/(rho*dx))*(sigxy[1:,1:-1]-sigxy[:-1,1:-1])
    
alpha = 2.05
K = (llmbda+mu)/3

nz = therm.grd.nz
def sigma(vx, vy, sigxx, sigyy, sigxy):
    #cells
    sigxx[:,:] = sigxx[:,:] + (dt*(llmbda+2*mu)/dx)*(vx[1:,:] - vx[:-1,:]) + (dt*llmbda/dy)*(vy[:,1:] - vy[:,:-1])-alpha*K*(therm.T[:,:,nz/2]-therm.T0)
    sigyy[:,:] = sigyy[:,:] + (dt*llmbda/dx)*(vx[1:,:] - vx[:-1,:]) + (dt*(llmbda+2*mu)/dy)*(vy[:,1:] - vy[:,:-1])-alpha*K*(therm.T[:,:,nz/2]-therm.T0)

    #inner corners
    sigxy[1:-1,1:-1] = sigxy[1:-1,1:-1] + (dt*mu/dy)*(vx[1:-1,1:] - vx[1:-1,:-1]) + (dt*mu/dx)*(vy[1:,1:-1] - vy[:-1,1:-1])

def boundaries(vx, vy, sigxx, sigyy, sigxy, t):
    #x = 0
    #-------------fixed-------------------------------
    #vx, sigxy = 0
    
    vx[0,:] = 0
    #vy[-0.5,:] = 0 = 0.5*(vy[-1,:] + vy[0,:])

    #sigxy requires fake point (vy[0,:] - vy[-1,:])/dx
    sigxy[0,:] = sigxy[0,:] + (dt*mu/dx)*2.0*vy[0,:]

    #x = -1
    #-------------free-------------------------------
    #sigxx, sigxy = 0
    
    sigxy[-1,:] = 0
    #sigxx[-1+1,:] = 2*alpha*K*(therm.T[:,:,nz/2]-therm.T0)- sigxx[-1,:] - sigxx[-1,:]
    #sigxx[-0.5,:] = 0 = 0.5*(sigxx[-1,:] + sigxx[-1+1,:])-alpha*K*(T-T0)

    
    #vx requires fake points 
    vx[-1,:] = vx[-1,:] + (dt/(rho*dx))*(2*alpha*K*(therm.T[-1,:,nz/2]-therm.T0)- sigxx[-1,:] - sigxx[-1,:])
    
    #FORCING SOURCE
     
    ##vx[-1,:] = vx[-1,:]+ source*np.exp(-t**2/a**2)

    #y = 0
    #-------------free--------------------------------
    #sigyy, sigxy = 0
    
    sigxy[:,0] = 0
    #sigyy[:,-0.5] = 0

    vy[:,0] = vy[:,0] + (dt/(rho*dx))*2.0*sigyy[:,0]

    #y = -1
    #-------------free-------------------------------
    #sigyy, sigxy = 0

    sigxy[:,0] = 0
    
    vy[:,-1] = vy[:,-1] + (dt/(rho*dx))*(-2.0)*sigyy[:,-1]

t = 0
def step(t):
    velocity(vx, vy, sigxx, sigyy, sigxy)
    sigma(vx, vy, sigxx, sigyy, sigxy)
    boundaries(vx, vy, sigxx, sigyy, sigxy, t)

def show(step):
    plt.figure()
    plt.subplot(1, 1, 1)
    plt.contourf([dy*i for i in range(grd.ny)], [dx*i for i in range(grd.nx+1)],vx)
    plt.gca().set_title(str(step))
    plt.show()

    
for i in range(1000):
    if i%50 == 0:
        print i
        show(i)
    therm.step(t)
    step(t)
    t += dt

