import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import fdtd as mech

import cattaneo as therm

import ini
grd = ini.Grid('x','y','z')


def addSubplot(fig, ax, arr):
    
    
    ax.plot(arr)
    #fig.canvas.draw_idle()
    

fig = plt.figure()
ax1 = fig.add_subplot(2, 2, 1)

ax2 = fig.add_subplot(2, 2, 2)

ax3 = fig.add_subplot(2, 2, 3)
plt.ion()
plt.show()


if __name__ == "__main__":
    t = 0
    for i in range(100):
        #if i == 23 or i == 117:
        if i%10 == 0:
            print i
            #mech.showSigT_1d(i,therm.T[:,:,grd.nz/2] )
            ax1.clear()
            addSubplot(fig,ax1, mech.sigxx[:,grd.ny/2])
            ax1.set_title('SIGMA_xx st={} t={:.2e}'.format(i, t))
            ax1.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
            ax1.set_ylim(min(mech.sigxx[:,grd.ny/2]), max(mech.sigxx[:,grd.ny/2]))

            ax2.clear()
            addSubplot(fig,ax2, mech.vx[:,grd.ny/2])
            ax2.set_title('v_x st={} t={:.2e}'.format(i, t))
            ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
            ax2.set_ylim(min(mech.vx[:,grd.ny/2]), max(mech.vx[:,grd.ny/2]))

            ax3.clear()
            addSubplot(fig,ax3, therm.T[:,grd.ny/2,grd.nz/2])
            ax3.set_title('T st={} t={:.2e}'.format(i, t))
            ax3.yaxis.set_major_formatter(FormatStrFormatter('%.5e'))
            ax3.set_ylim(min(therm.T[:,grd.ny/2,grd.nz/2]), max(therm.T[:,grd.ny/2,grd.nz/2]))
            fig.canvas.draw()
            #plt.draw()
            raw_input()
        therm.step(t)
        mech.step(t, therm.T[:,:,grd.nz/2])
        t += grd.dt
