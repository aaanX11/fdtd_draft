import math
import numpy as np


def write_to_file(grd, therm, mech):
    import os
    path = raw_input("Path to tm? ")

    with open(os.path.join(path, "therm_ini"), 'w') as f:
        f.write(str(therm["rho"]))
        f.write('\n')
        f.write(str(therm["C_V"]))
        f.write('\n')
        f.write(str(therm["k"]))
        f.write('\n')
        f.write(str(therm["tau"]))
        f.write('\n')

    with open(os.path.join(path, "mech_ini"), 'w') as f:
        f.write(str(mech["rho"]))
        f.write('\n')
        f.write(str(mech["C11"]))
        f.write('\n')
        f.write(str(mech["C12"]))
        f.write('\n')
        f.write(str(mech["C44"]))
        f.write('\n')
        f.write(str(mech["alpha"]))
        f.write('\n')

    with open(os.path.join(path, "grid_ini"), 'w') as f:
        f.write(str(grd.nx))
        f.write('\n')
        f.write(str(grd.ny))
        f.write('\n')
        f.write(str(grd.nz))
        f.write('\n')
        f.write(str(grd.dx))
        f.write('\n')
        f.write(str(grd.dy))
        f.write('\n')
        f.write(str(grd.dz))
        f.write('\n')
        f.write(str(1e5))  # nt
        f.write('\n')
        f.write(str(grd.dt))
        f.write('\n')


def chelkanov2019():
    global x_ph, dt
    # ------------------------------------------
    alpha_x = 2.05e-6  # K^-1

    C11 = 1.66e11  # N m^-2
    C12 = 0.639e11
    C44 = 0.5 * (C11 - C12)  # isotropic

    alpha = alpha_x

    # ------------------------------------------
    C_V = 1.8e6  # J m^-3 K^-1
    rho = 2330.  # kg/m^3
    k = 1.5e2  # W m^-1 K^-1

    c_str = 8440  # m s-1

    tau = 1.17e-12  # s

    dt = tau*0.1
    x_ph = 0.0001*1e-2  # m

    # thermal shock-----------------
    # deltaT(x) : 1e-7, 5e-7, 1e-6 s
    # sigma(x)  : 5e-7, 1e-6
    # q(x)      : 1e-7, 5e-7, 1e-6 s

    # laser-----------------
    # T(x)      : 1e-6, 3e-6
    # sig(x)    : 1e-6, 3e-6
    # T(t) left boundary
    # T(t) middle

    def source_laser(x):
        q0 = 4.5e5  # Wt m^-2
        tau_ = 3.e-6

        # http://www.miigaik.ru/upload/iblock/2fb/2fb9d53e46a7aae7eff545ce1ae56bfc.pdf :
        alpha_ = 1e3  # sm^-1
        return q0*alpha_*np.exp(-alpha_*x)/tau_

    def boundary_source_laser(t):
        return 0

    def source_thermal_shock(x):
        return 0

    def boundary_source_thermal_shock(t):
        return 300+300*(1-np.exp(-t/(1.528402175606507e-06)))

    sources = [source_laser, boundary_source_laser]
    sources = [source_thermal_shock, boundary_source_thermal_shock]

    return ('chelkanov',
            {'rho': rho, 'C_V': C_V, 'k': k, 'tau': tau},
            {'rho': rho, 'C11': C11, 'C12': C12, 'C44': C44, 'alpha': alpha},
            sources)


def xiong2011():
    global x_ph, dt
    #------------------------------------------
    lmbda_x = 9.4e10 #N m^-2
    mu_x = 4.0e10 #N m^-2
    alpha_x = 2.05e-5 #K^-1

    C11_ph = lmbda_x+2*mu_x
    C12_ph = lmbda_x
    C44_ph = 0.5*(C11_ph - C12_ph)#isotropic 
    alpha_ph = 3*alpha_x

    #------------------------------------------
    c_E = 1.04e3 #J kg^-1 K^-1
    rho_x = 1740. #kg/m^3
    k_x = 1.7e2 #W m^-1 K^-1
    
    c_str2 = (lmbda_x+2*mu_x)/rho_x
    print math.sqrt(c_str2)
    eta = rho_x*c_E/k_x
    tau_star  = 0.05

    tau_x = tau_star/(c_str2*eta)
    #print tau
    #4.7e-14

    x_star = 3
    x = x_star/(math.sqrt(c_str2)*eta)
    #print x 
    #x_ph = 2.8e-8
    x_ph = 2.8e-7
    nx = 20
    #------------------------------------------
    

    C_V_ph = c_E*rho_x
    #------------------------------------------
    dt = tau_x*0.1

    def source_thermal_shock(x):
        return 0

    def boundary_source_thermal_shock(t):
        tau = 1e-12
        return 300 + 300 * (1 - np.exp(-t / tau))

    sources = [source_thermal_shock, boundary_source_thermal_shock]

    return ('xi',
            {'rho': rho_x, 'C_V': C_V_ph, 'k': k_x, 'tau': tau_x},
            {'rho': rho_x, 'C11': C11_ph, 'C12': C12_ph, 'C44': C44_ph, 'alpha': alpha_ph},
            sources)


def reutov():
    global x_ph, dt

    C_V = 1.8e6  # J m^-3 K^-1
    rho = 2330.  # kg/m^3
    k = 1.5e2  # W m^-1 K^-1
    alpha = 3 * alpha_x



    # ------------------------------------------

    C11 = 1.66e11  # N m^-2
    C12 = 0.639e11
    C44 = 0.5 * (C11 - C12)  # isotropic


    # ------------------------------------------



    c_str2 = 8440  # m s-1

    tau = 2e-12  # s

    dt = tau*0.1
    x_ph = 0.0001*1e-2  # m

    # thermal shock-----------------
    # deltaT(x) : 1e-7, 5e-7, 1e-6 s
    # sigma(x)  : 5e-7, 1e-6
    # q(x)      : 1e-7, 5e-7, 1e-6 s

    # laser-----------------
    # T(x)      : 1e-6, 3e-6
    # sig(x)    : 1e-6, 3e-6
    # T(t) left boundary
    # T(t) middle

    def source_laser(x):
        q0 = 4.5e5  # Wt m^-2
        tau_ = 3.e-6

        # http://www.miigaik.ru/upload/iblock/2fb/2fb9d53e46a7aae7eff545ce1ae56bfc.pdf :
        alpha_ = 1e3  # sm^-1
        return q0*alpha_*np.exp(-alpha_*x)/tau_

    def boundary_source_laser(t):
        return 0

    def source_thermal_shock(x):
        return 0

    def boundary_source_thermal_shock(t):
        return 300+300*(1-np.exp(-t/tau))

    sources = [source_laser, boundary_source_laser]
    sources = [source_thermal_shock, boundary_source_thermal_shock]

    return ({'rho': rho, 'C_V': C_V, 'k': k, 'tau': tau},
            {'rho': rho, 'C11': C11, 'C12': C12, 'C44': C44, 'alpha': alpha},
            sources)

class Grid:
    def __init__(self, x, y, z):
        print x_ph, dt
        self.nx = 100
        self.ny = 20
        self.nz = 23
        self.xi = np.linspace(0, x_ph, self.nx)
        self.yi = np.linspace(0, x_ph, self.ny)

        self.dxi05 = [i-j for i, j in zip(self.xi[1:], self.xi)]
        self.dxi05.insert(0, self.dxi05[0])
        self.dxi05.append(self.dxi05[-1])
        self.xi05 = [i-0.5*j for i,j in zip(self.xi, self.dxi05)]
        self.xi05.append(self.xi[-1]+0.5*self.dxi05[-1])
        self.dx = min(self.dxi05)
        self.dy = self.dx
        self.dz = self.dx

        self.dt = dt


def from_file():
    # data = np.loadtxt('C:\ks_work\calculations\cattaneo_py_test\data1')
    data = np.loadtxt('data_thermo', ndmin=2)

    rho = data[:, 0]
    C_V = data[:, 1]
    k = data[:, 2]
    tau = data[:, 3]

    data = np.loadtxt('data_elastic', ndmin=2)

    rho2 = data[:, 0]
    C11_ = data[:, 1]
    C12_ = data[:, 2]
    C44_ = data[:, 3]
    alpha1 = data[:, 4]

x_ph = 1
dt = 1
if __name__ == '__main__':

    if raw_input('calculation? x or ch') == 'x':
        r = xiong2011()
    else:
        r = chelkanov2019()

    print r[0]
    print r[1]
    grd = Grid('x', 'y', 'z')

    write_to_file(grd, r[0], r[1])

