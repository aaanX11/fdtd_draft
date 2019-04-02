import numpy as np
import matplotlib.pyplot as plt

import ini
_, mech = ini.xiong2011()
rho = mech['rho']

C11 = mech['C11']
C12 = mech['C12']
C44 = mech['C44']
alpha = mech['alpha']

grd = ini.Grid('x','y','z')
dx = grd.dx
dy = grd.dy
dz = grd.dz
dt = grd.dt

K = (C11+2*C12)/3

T0 = 300


class Space:
    def __init__(self, nx, ny):
        self.cell = np.zeros((nx,ny))

s = Space(grd.nx, grd.ny)

vx = np.zeros((grd.nx+1,grd.ny))#velocities - edges
vy = np.zeros((grd.nx,grd.ny+1))

sigxx = np.zeros((grd.nx,grd.ny))#p, lambda, mu - centers
sigyy = np.zeros((grd.nx,grd.ny))

sigxy = np.zeros((grd.nx+1,grd.ny+1))#shear stress, density - corners


def get_source():
    a = dx
    import math
    vec = np.asarray([dy*(i-grd.ny/2) for i in range(grd.ny)])
    A = 2/np.sqrt(3*a)* math.pow(math.pi,1./3)
    res = A*(1 - vec**2/a**2)*np.exp(-vec**2/a**2)
    plt.plot( res)
    plt.show()
    return res#A*(1 - vec**2/a**2) 
    
#source = get_source()

def velocity(vx, vy, sigxx, sigyy, sigxy):
    #nx+1 ny
    #inner edges X cells
    vx[1:-1,:] = vx[1:-1,:] + (dt/(rho*dx))*(sigxx[1:,:]-sigxx[:-1,:]) + (dt/(rho*dy))*(sigxy[1:-1,1:]-sigxy[1:-1,:-1])
    #nx ny+1
    #cells X inner edges
    vy[:,1:-1] = vy[:,1:-1] + (dt/(rho*dy))*(sigyy[:,1:]-sigyy[:,:-1]) + (dt/(rho*dx))*(sigxy[1:,1:-1]-sigxy[:-1,1:-1])
    

def sigma(vx, vy, sigxx, sigyy, sigxy, T):
    #cells
    sigxx[:,:] = sigxx[:,:] + (dt*(C11)/dx)*(vx[1:,:] - vx[:-1,:]) + (dt*C12/dy)*(vy[:,1:] - vy[:,:-1])-alpha*K*(T-T0)
    sigyy[:,:] = sigyy[:,:] + (dt*C12/dx)*(vx[1:,:] - vx[:-1,:]) + (dt*(C11)/dy)*(vy[:,1:] - vy[:,:-1])-alpha*K*(T-T0)

    #inner corners
    sigxy[1:-1,1:-1] = sigxy[1:-1,1:-1] + (dt*C44/dy)*(vx[1:-1,1:] - vx[1:-1,:-1]) + (dt*C44/dx)*(vy[1:,1:-1] - vy[:-1,1:-1])

def boundaries(vx, vy, sigxx, sigyy, sigxy, t, T):
    #x = 0
    #-------------fixed-------------------------------
    #vx, sigxy = 0
    
    vx[0,:] = 0
    #vy[-0.5,:] = 0 = 0.5*(vy[-1,:] + vy[0,:])

    #sigxy requires fake point (vy[0,:] - vy[-1,:])/dx
    sigxy[0,:] = sigxy[0,:] + (dt*C44/dx)*2.0*vy[0,:]

    #x = -1
    #-------------free-------------------------------
    #sigxx, sigxy = 0
    
    sigxy[-1,:] = 0
    #sigxx[-1+1,:] = 2*alpha*K*(therm.T[:,:,nz/2]-therm.T0)- sigxx[-1,:] - sigxx[-1,:]
    #sigxx[-0.5,:] = 0 = 0.5*(sigxx[-1,:] + sigxx[-1+1,:])-alpha*K*(T-T0)

    
    #vx requires fake points 
    vx[-1,:] = vx[-1,:] + (dt/(rho*dx))*(2*alpha*K*(T[-1,:]-T0)- sigxx[-1,:] - sigxx[-1,:])
    
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

#t = 0
def step(t, T):
    velocity(vx, vy, sigxx, sigyy, sigxy)
    sigma(vx, vy, sigxx, sigyy, sigxy, T)
    boundaries(vx, vy, sigxx, sigyy, sigxy, t, T)

def show(step):
    plt.figure()
    plt.subplot(1, 1, 1)
    plt.contourf([dy*i for i in range(grd.ny)], [dx*i for i in range(grd.nx+1)],vx)
    plt.gca().set_title('Vx'+str(step))
    plt.show()

def showSig(step):
    plt.figure()
    plt.subplot(1, 1, 1)
    plt.contourf([dy*i for i in range(grd.ny)], [dx*i for i in range(grd.nx)],sigxx)
    plt.gca().set_title('SIGMA_xx'+str(step))
    plt.show()

def showSig_1d(step):
    plt.figure()
    plt.subplot(1, 1, 1)
    plt.plot(sigxx[:,grd.ny/2])
    #plt.contourf([dy*i for i in range(grd.ny)], [dx*i for i in range(grd.nx)],sigxx)
    plt.gca().set_title('SIGMA_xx'+str(step))
    plt.show()

def showSigT_1d(step, T):
    plt.figure()
    plt.subplot(1, 1, 1)
    plt.plot(sigxx[:,grd.ny/2])
    #plt.contourf([dy*i for i in range(grd.ny)], [dx*i for i in range(grd.nx)],sigxx)
    plt.gca().set_title('SIGMA_xx'+str(step))
    plt.show()


