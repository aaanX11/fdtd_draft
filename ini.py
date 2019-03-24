import math

lmbda = 9.4e10 #N m^-2
mu = 4.0e10 #N m^-2
alpha = 0.5e10 #N m^-2
rho = 1740. #kg/m^3
k = 1.7e2 #W m^-1 K^-1
c_E = 1.04e3 #J kg^-1 K^-1
c_str2 = (lmbda+2*mu)/rho
eta = rho*c_E/k
tau_star  = 0.05

tau = tau_star/(c_str2*eta)
#print tau
#4.7e-14

x_star = 3
x = x_star/(math.sqrt(c_str2)*eta)
#print x 
x = 20*2.8e-8

#??????????????????????????
k0 = 25 #Wt/m
