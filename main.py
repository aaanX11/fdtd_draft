import numpy as np
import matplotlib
matplotlib.use("TkAgg")

import time
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import Tkinter as tk
import sys

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
if matplotlib.__version__ < "2.2.0":
    from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg as NavigationToolbar2Tk
else:
    from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
from matplotlib.backend_bases import key_press_handler

import fdtd as mech
import cattaneo as therm
import ini

grd = ini.Grid('x', 'y', 'z')


def addSubplot(fig, ax, arr, title):
    ax.clear()
    ax.plot(arr)
    

    Marr = max(arr)
    marr = min(arr)
    #print Marr, marr, title
    ampl = Marr-marr
    if ampl == 0:
        ampl = 1e-7
    ax.set_ylim(marr-0.01*ampl, Marr+0.01*ampl)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax.set_title(title)


def show4(fig, ax1, ax2, ax3, ax4):

    fig.suptitle(' st={} t={:.2e}'.format(i, t))
    
    addSubplot(fig, ax1, mech.sigxx[:,grd.ny/2], 'SIGMA_xx')
    addSubplot(fig, ax2, mech.vx[:,grd.ny/2], 'v_x')
    addSubplot(fig, ax3, therm.T[:,grd.ny/2,grd.nz/2], 'T')
    #addSubplot(fig, ax3, therm.T_der[:, grd.ny/2, grd.nz/2], 'T\'')
    addSubplot(fig, ax4, therm.qx[:,grd.ny/2,grd.nz/2], 'qx')

    fig.canvas.draw()


def show(fig, ax1, ax2, ax3):

    fig.suptitle(' st={} t={:.2e}'.format(i, t))
    
    addSubplot(fig, ax1, mech.sigxx[:,grd.ny/2], 'SIGMA_xx')
    addSubplot(fig, ax2, mech.vx[:,grd.ny/2], 'v_x')
    addSubplot(fig, ax3, therm.T[:,grd.ny/2,grd.nz/2], 'T')
    #addSubplot(fig, ax3, therm.T_der[:, grd.ny/2, grd.nz/2], 'T\'')

    #addSubplot(fig, ax4, therm.T[:,grd.ny/2,grd.nz/2], 'T')

    fig.canvas.draw()


i = 0
t = 0


def step_series(steps):
    print 'starting step series'
    global i
    j = i + steps
    while i <= j:
        print i
        i += 1
        step()
    
    if i > 11000:
        window.quit()
        window.destroy()
        sys.exit(0)
    show4(fig, ax1, ax2, ax3, ax4)


def step():
    global t
    # if i == 23 or i == 117:
    start_time = time.time()
    therm.step(t)
    print("--- therm %s seconds ---" % (time.time() - start_time))

    start_time = time.time()
    t += grd.dt
    if False:  # builtin spline
        mech.step(t, np.resize(therm.T_der_[:, :, grd.nz/2, :],(grd.nx, grd.ny, 6)))
    else:  # quad spline
        mech.step(t, therm.T[:, :, grd.nz / 2])

    print("--- mech  %s seconds ---" % (time.time() - start_time))
    # print therm.T_der_[0, 0, grd.nz/2, :]



def onclick(event):
    start_time = time.time()
    step_series(10)
    print("--- %s seconds ---" % (time.time() - start_time))


def formula(ax):
    ax.axis('off')
    ax.xticks=()
    ax.yticks=()
    
    strings = [r"$\frac{\partial \sigma_{11}}{\partial t} = C_{11}\frac{\partial v_1}{\partial x_1} + C_{12}\frac{\partial v_2}{\partial x_2} - \alpha K\dot{T}$",
               r"$\frac{\partial \sigma_{22}}{\partial t} = C_{12}\frac{\partial v_1}{\partial x_1} + C_{11}\frac{\partial v_2}{\partial x_2} - \alpha K\dot{T}$",
               r"$\frac{\partial \sigma_{12}}{\partial t} =\frac{C_{44}}{2}\left(\frac{\partial v_1}{\partial x_2} +\frac{\partial v_2}{\partial x_1}\right)$",
               r"$\rho\frac{\partial v_1}{\partial t} = \frac{\partial \sigma_{11}}{\partial x_1} + \frac{\partial \sigma_{12}}{\partial x_2}$",
               r"$\rho\frac{\partial v_2}{\partial t} = \frac{\partial \sigma_{12}}{\partial x_1} + \frac{\partial \sigma_{22}}{\partial x_2}$"
               ]
    nn = len(strings)
    for y, string in zip([0.95-(1.0/(nn+1))*i for i in range(nn)], strings):
        ax.text(0.05, y, string, fontsize=16)
        #print string
    strings = [r"$\frac{\partial w}{\partial t} +div \mathbb{q} = 0$",
               r"$\tau \frac{\partial \mathbb{q}}{\partial t} +\lambda\nabla T = -\mathbb{q}$"
               ]
    for y, string in zip([0.95-(1.0/(nn+1))*i for i in range(nn)], strings):
        ax.text(0.6, y, string, fontsize=16)


if __name__ == "__main__":
    window = tk.Tk()
    
    fig = plt.figure()
    ax1 = fig.add_subplot(2, 2, 1)
    ax2 = fig.add_subplot(2, 2, 2)
    ax3 = fig.add_subplot(2, 2, 3)
    ax4 = fig.add_subplot(2, 2, 4)

    formula(ax4)
    # show4(fig, ax1, ax2, ax3)
    show4(fig, ax1, ax2, ax3, ax4)

    fig_canvas_agg = FigureCanvasTkAgg(fig, master=window)
    fig_canvas_agg.draw()
    fig_canvas_agg.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    toolbar = NavigationToolbar2Tk(fig_canvas_agg, window)
    toolbar.update()
    fig_canvas_agg._tkcanvas.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=1)

    fig_canvas_agg.mpl_connect('button_press_event', onclick)

    window.protocol("WM_DELETE_WINDOW", lambda : window.quit())
    window.mainloop()
    window.destroy()
