import fdtd as mech

import cattaneo as therm

t = 0
for i in range(1000):
    if i%50 == 0:
        print i
        mech.show(i)
    therm.step(t)
    mech.step(t)
    t += ini.dt
