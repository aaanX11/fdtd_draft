import math
import numpy as np

def xiong2011():
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
    eta = rho_x*c_E/k_x
    tau_star  = 0.05

    tau_x = tau_star/(c_str2*eta)
    #print tau
    #4.7e-14

    x_star = 3
    x = x_star/(math.sqrt(c_str2)*eta)
    #print x 
    x_ph = 2.8e-8
    nx = 20
    #------------------------------------------
    

    C_V_ph = c_E*rho_x
    #------------------------------------------
    dt = 3*tau_x
    

    
    return ({'rho' : rho_x, 'C_V' : C_V_ph, 'k' : k_x, 'tau' : tau_x}, {'C11' : C11_ph, 'C12' : C12_ph, 'C44' : C44_ph, 'alpha' : alpha_ph}, [i*x_ph for i in range(nx)], dt)


def from_file():
    #data = np.loadtxt('C:\ks_work\calculations\cattaneo_py_test\data1')
    data = np.loadtxt('data_thermo', ndmin = 2)   

    rho = data[:,0]
    C_V = data[:,1]
    k = data[:,2]
    tau = data[:,3]

    data = np.loadtxt('data_elastic', ndmin = 2) 

    rho2 = data[:,0]
    C11_ = data[:,1]
    C12_ = data[:,2]
    C44_ = data[:,3]
    alpha1 = data[:,4]

if __name__ == '__main__':
    r = xiong2011()
    print r[0]
    print r[1]
    print r[2]
    print r[3]
 
