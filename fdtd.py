import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import UnivariateSpline

# -*- coding: cp1251 -*-
import ini
# _, mech, _ = ini.chelkanov2019()
_, _, mech, _ = ini.xiong2011()
rho = mech['rho']

C11 = mech['C11']
C12 = mech['C12']
C44 = mech['C44']
alpha = mech['alpha']

grd = ini.Grid('x', 'y', 'z')
dx = grd.dx
dy = grd.dy
dz = grd.dz
dt = grd.dt
dt_loc = 0.0001*min(dx, dy, dz)/np.sqrt(C11/rho)

K = (C11+2*C12)/3

T0 = 300


class Space:
    def __init__(self, nx, ny):
        self.cell = np.zeros((nx,ny))

s = Space(grd.nx, grd.ny)

vx = np.zeros((grd.nx+1, grd.ny))#velocities - edges
vy = np.zeros((grd.nx, grd.ny+1))

sigxx = np.zeros((grd.nx, grd.ny))#p, lambda, mu - centers
sigyy = np.zeros((grd.nx, grd.ny))

T_der_0 = np.zeros((grd.nx, grd.ny))
T_der2_0 = np.zeros((grd.nx, grd.ny))
T_cur = np.zeros((grd.nx, grd.ny))
T_old = np.zeros((grd.nx, grd.ny))

sigxy = np.zeros((grd.nx+1, grd.ny+1))#shear stress, density - corners


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
    vx[1:-1,:] = vx[1:-1,:] + (dt_loc/(rho*dx))*(sigxx[1:,:]-sigxx[:-1,:]) + (dt_loc/(rho*dy))*(sigxy[1:-1,1:]-sigxy[1:-1,:-1])
    #nx ny+1
    #cells X inner edges
    vy[:,1:-1] = vy[:,1:-1] + (dt_loc/(rho*dy))*(sigyy[:,1:]-sigyy[:,:-1]) + (dt_loc/(rho*dx))*(sigxy[1:,1:-1]-sigxy[:-1,1:-1])

def derivative(arr, *args):
    x0 = args[0]
    x = arr[:3]
    f = arr[3:]

    # http://fourier.eng.hmc.edu/e176/lectures/ch6/node3.html

    xx1 = np.asarray([2*x0 - x[1] - x[2], 2*x0 - x[0] - x[2], 2*x0 - x[0] - x[1]])
    xx2 = np.asarray([(x[0] - x[1])*(x[0] - x[2]), (x[1] - x[0])*(x[1] - x[2]), (x[2] - x[0])*(x[2] - x[1])])

    #return np.sum((xx1/xx2)*f)
    return f[0]*(2*x0 - x[1] - x[2])/((x[0] - x[1])*(x[0] - x[2])) \
           + f[1]*(2*x0 - x[0] - x[2])/((x[1] - x[0])*(x[1] - x[2])) \
           + f[2]*(2*x0 - x[0] - x[1])/((x[2] - x[0])*(x[2] - x[1]))


def derivative_spline(arr, *args):
    t = args[0]
    spl = UnivariateSpline(arr[:3], arr[3:], k=2)
    der = spl.derivative()
    return der([t])[0]


def sigma(vx, vy, sigxx, sigyy, sigxy, t, T):
    vx_d = vx[1:, :] - vx[:-1, :]
    
    # print 'sigma_terms\nalpha *K *derivative = ', alpha*K*T[-3:-1, grd.ny/2]
    # print 'C11*vx_deriv_x = ', C11*vx_d[-3:-1,grd.ny/2]/dx

    # T_der[:, :] = np.apply_along_axis(derivative_spline, 2, T, t)

    #cells
    sigxx[:, :] = sigxx[:, :] + (dt_loc*(C11)/dx)*(vx[1:, :] - vx[:-1, :]) + (dt_loc*C12/dy)*(vy[:, 1:] - vy[:, :-1])-dt_loc*alpha*K*T
    sigyy[:, :] = sigyy[:, :] + (dt_loc*C12/dx)*(vx[1:, :] - vx[:-1, :]) + (dt_loc*(C11)/dy)*(vy[:, 1:] - vy[:, :-1])-dt_loc*alpha*K*T

    #inner corners
    sigxy[1:-1,1:-1] = sigxy[1:-1,1:-1] + (dt_loc*C44/dy)*(vx[1:-1,1:] - vx[1:-1,:-1]) + (dt_loc*C44/dx)*(vy[1:,1:-1] - vy[:-1,1:-1])

def boundaries(vx, vy, sigxx, sigyy, sigxy, t, T):
    #x = 0
    #-------------free-------------------------------
    #sigxx, sigxy = 0
    
    sigxy[0,:] = 0
    vx[0,:] = vx[0,:] + dt_loc*(sigxx[0,:] + sigxx[0,:])/(dx*rho)

    #x = -1
    #-------------fixed-------------------------------
    #vx, vy = 0
    
    sigxy[-1,:] = sigxy[-1,:] - (dt_loc*C44/dx)*2.0*vy[-1,:]
    vx[-1,:] = 0
    
    #y = 0
    #-------------free--------------------------------
    #sigyy, sigxy = 0
    
    sigxy[:,0] = 0
    #sigyy[:,-0.5] = 0

    vy[:,0] = vy[:,0] + (dt_loc/(rho*dx))*2.0*sigyy[:,0]

    #y = -1
    #-------------free-------------------------------
    #sigyy, sigxy = 0

    sigxy[:, 0] = 0
    
    vy[:, -1] = vy[:,-1] + (dt_loc/(rho*dx))*(-2.0)*sigyy[:,-1]

t_loc = 0
def step_loc_old(t, T):
    #print 'm'
    velocity(vx, vy, sigxx, sigyy, sigxy)
    sigma(vx, vy, sigxx, sigyy, sigxy, t, T)
    boundaries(vx, vy, sigxx, sigyy, sigxy, t, T)

def step_old(t, T):
    global t_loc
    global dt
    #print dt_loc, t - t_loc
    while t_loc < t:
        dt = min(dt_loc, t - t_loc)
        #print "mech dt = ", dt
        #raw_input()
        step_loc_old(t_loc+dt, T)
        t_loc += dt

def derivative(dt, dT):
    return

def step_loc(t, T):
    #print 'm'
    velocity(vx, vy, sigxx, sigyy, sigxy)
    sigma(vx, vy, sigxx, sigyy, sigxy, t, T)
    boundaries(vx, vy, sigxx, sigyy, sigxy, t, T)

def step(t, T):
    global t_loc, dt

    global T_cur, T_der_0

    T_old = T_cur
    T_cur = T - T_cur
    T_der_i = T_der_0

    dt_outer = t - t_loc
    while t_loc < t:
        dt = min(dt_loc, t - t_loc)
        step_loc(t_loc+dt, T_der_i)
        print dt_loc, dt_outer, np.amax(T_cur), np.amax(T_old), np.amax(T_der_0)
        T_der_i += (2*dt_loc/(dt_outer*dt_outer))*T_cur - (2*dt_loc/dt_outer)*T_der_0

        t_loc += dt

    T_der_0 = 2 * T_cur / dt_outer - T_der_0
    T_cur += T_old


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


