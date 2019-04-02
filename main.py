import fdtd as mech

import cattaneo as therm

import ini
grd = ini.Grid('x','y','z')

t = 0
for i in range(100):
    #if i == 23 or i == 117:
    if i%10 == 0:
        print i
        mech.showSigT_1d(i,therm.T[:,:,grd.nz/2] )
    therm.step(t)
    mech.step(t, therm.T[:,:,grd.nz/2])
    t += grd.dt
