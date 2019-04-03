import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

import Tkinter as tk
import sys

matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2TkAgg)
from matplotlib.backend_bases import key_press_handler

import fdtd as mech
import cattaneo as therm
import ini

grd = ini.Grid('x','y','z')


def addSubplot(fig, ax, arr, title):
    #ax.clear()
    ax.plot(arr)
    

    Marr = max(arr)
    marr = min(arr)
    ampl = Marr-marr
    if ampl == 0:
        ampl = 1e-7
    ax.set_ylim(marr-0.01*ampl, Marr+0.01*ampl)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax.set_title(title)
    

    
def show(fig, ax1, ax2, ax3):

    fig.suptitle(' st={} t={:.2e}'.format(i, t))
    
    addSubplot(fig,ax1, mech.sigxx[:,grd.ny/2], 'SIGMA_xx')
    addSubplot(fig,ax2, mech.vx[:,grd.ny/2], 'v_x')
    addSubplot(fig,ax3, therm.T[:,grd.ny/2,grd.nz/2], 'T')

    fig.canvas.draw()
    
i = 0
t = 0

def step_series(steps):
    global i
    j = i + steps
    while i <= j:
        i += 1
        step()
    
    if i > 117:
        window.quit()
        window.destroy()
        sys.exit(0)
    show(fig, ax1, ax2, ax3)

def step():
    global t
    #if i == 23 or i == 117:
    therm.step(t)
    mech.step(t, therm.T[:,:,grd.nz/2])
    t += grd.dt
    
def onclick(event):
    step_series(10)
    
def formula(ax):
    pass


if __name__ == "__main__":
    window = tk.Tk()
    
    fig = plt.figure()
    ax1 = fig.add_subplot(2, 2, 1)
    ax2 = fig.add_subplot(2, 2, 2)
    ax3 = fig.add_subplot(2, 2, 3)
    ax4 = fig.add_subplot(2, 2, 4)

    formula(ax4)
    show(fig, ax1, ax2, ax3)

    fig_canvas_agg = FigureCanvasTkAgg(fig, master = window)
    fig_canvas_agg.draw()
    fig_canvas_agg.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    toolbar = NavigationToolbar2TkAgg(fig_canvas_agg, window)
    toolbar.update()
    fig_canvas_agg._tkcanvas.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=1)

    fig_canvas_agg.mpl_connect('button_press_event', onclick)


    window.protocol("WM_DELETE_WINDOW", lambda : window.quit())
    window.mainloop()
    window.destroy()  
    
        

    
        
