import ini

import numpy as np
import matplotlib

matplotlib.use("TkAgg")

import time
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

if matplotlib.__version__ < "2.2.0":
    from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg as NavigationToolbar2Tk
else:
    from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
from matplotlib.backend_bases import key_press_handler

import Tkinter as tk
import ttk

import json
import sys
import os


grd = ini.Grid('x', 'y', 'z')

def addSubplot(fig, ax, arr, title):
    ax.clear()
    ax.plot(arr)

    Marr = max(arr)
    marr = min(arr)
    # print Marr, marr, title
    ampl = Marr - marr
    if ampl == 0:
        ampl = 1e-7
    ax.set_ylim(marr - 0.01 * ampl, Marr + 0.01 * ampl)
    if ampl > 1e3:
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax.set_title(title)


def show4(fig, ax1, ax2, ax3, ax4):
    fig.suptitle(' st={} t={:.2e}'.format(i, t))

    addSubplot(fig, ax1, mech.sigxx[:, grd.ny / 2], 'SIGMA_xx')
    addSubplot(fig, ax2, mech.vx[:, grd.ny / 2], 'v_x')
    addSubplot(fig, ax3, therm.T[:, grd.ny / 2, grd.nz / 2], 'T')
    # addSubplot(fig, ax3, therm.T_der[:, grd.ny/2, grd.nz/2], 'T\'')
    addSubplot(fig, ax4, therm.qx[:, grd.ny / 2, grd.nz / 2], 'qx')

    fig.canvas.draw()


def show(fig, ax1, ax2, ax3):
    fig.suptitle(' st={} t={:.2e}'.format(i, t))

    addSubplot(fig, ax1, mech.sigxx[:, grd.ny / 2], 'SIGMA_xx')
    addSubplot(fig, ax2, mech.vx[:, grd.ny / 2], 'v_x')
    addSubplot(fig, ax3, therm.T[:, grd.ny / 2, grd.nz / 2], 'T')
    # addSubplot(fig, ax3, therm.T_der[:, grd.ny/2, grd.nz/2], 'T\'')

    # addSubplot(fig, ax4, therm.T[:,grd.ny/2,grd.nz/2], 'T')

    fig.canvas.draw()


def formula(ax):
    ax.axis('off')
    ax.xticks = ()
    ax.yticks = ()

    strings = [
        r"$\frac{\partial \sigma_{11}}{\partial t} = C_{11}\frac{\partial v_1}{\partial x_1} + C_{12}\frac{\partial v_2}{\partial x_2} - \alpha K\dot{T}$",
        r"$\frac{\partial \sigma_{22}}{\partial t} = C_{12}\frac{\partial v_1}{\partial x_1} + C_{11}\frac{\partial v_2}{\partial x_2} - \alpha K\dot{T}$",
        r"$\frac{\partial \sigma_{12}}{\partial t} =\frac{C_{44}}{2}\left(\frac{\partial v_1}{\partial x_2} +\frac{\partial v_2}{\partial x_1}\right)$",
        r"$\rho\frac{\partial v_1}{\partial t} = \frac{\partial \sigma_{11}}{\partial x_1} + \frac{\partial \sigma_{12}}{\partial x_2}$",
        r"$\rho\frac{\partial v_2}{\partial t} = \frac{\partial \sigma_{12}}{\partial x_1} + \frac{\partial \sigma_{22}}{\partial x_2}$"
        ]
    nn = len(strings)
    for y, string in zip([0.95 - (1.0 / (nn + 1)) * i for i in range(nn)], strings):
        ax.text(0.05, y, string, fontsize=16)
        # print string
    strings = [r"$\frac{\partial w}{\partial t} +div \mathbb{q} = 0$",
               r"$\tau \frac{\partial \mathbb{q}}{\partial t} +\lambda\nabla T = -\mathbb{q}$"
               ]
    for y, string in zip([0.95 - (1.0 / (nn + 1)) * i for i in range(nn)], strings):
        ax.text(0.6, y, string, fontsize=16)


def draw_ini_vals(ax):
    ax.axis('off')
    ax.xticks = ()
    ax.yticks = ()

    strings = [
        r"$\alpha = $",
        r"$T = $",
        r"$C_V$",
        r"$\lambda =$",
        r"$\tau =$",
        r"$C_{11}$",
        r"$C_{44}$",
        r"$C_{12}$"
        ]
    nn = len(strings)
    for y, string in zip([0.95 - (1.0 / (nn + 1)) * i for i in range(nn)], strings):
        ax.text(0.05, y, string, fontsize=16)
        # print string
    strings = [r"$\frac{\partial w}{\partial t} +div \mathbb{q} = 0$",
               r"$\tau \frac{\partial \mathbb{q}}{\partial t} +\lambda\nabla T = -\mathbb{q}$"
               ]
    for y, string in zip([0.95 - (1.0 / (nn + 1)) * i for i in range(nn)], strings):
        ax.text(0.6, y, string, fontsize=16)


class FigFrame(ttk.Frame):
    def __init__(self, master, data):
        self.master = master
        ttk.Frame.__init__(self, master)

        #frame_main = tk.Frame(self)
        #frame_main.pack(side=tk.TOP)

        self.data = data

        self.fig = plt.figure(figsize=(12, 10))

        fig_canvas_agg = FigureCanvasTkAgg(self.fig, master=self)
        fig_canvas_agg.draw()
        fig_canvas_agg.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        toolbar = NavigationToolbar2Tk(fig_canvas_agg, self)
        toolbar.update()
        # fig_canvas_agg.get_tk_widget().grid(row=5, column =0, columnspan=4)
        fig_canvas_agg._tkcanvas.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=1)

        self.axes = self.fig.subplots(2, 2)

        self.arrays = dict()
        self.read(data)


        if data['1d']:
            data['prj_name']
            data['n']
            data['x']
            data['nt']
            data['t']
            data['field_names'] = ('cells', 'w', 'q', 'T', 'T\'', 'v', 'sig')

        elif data['2d']:
            data['field_names'] = ('cells', 'w', 'q_x', 'q_y', 'T', 'T\'', 'v_x', 'v_y', 'sig_xx', 'sig_yy', 'sig_xy')

    def draw_pics(self, fields, vars):
        for ax, field_name, var in zip(self.axes, fields, vars):
            field_idx = self.data['field_names'].index(field_name)
            if var == 'x':
                ydata = self.arrays[field_idx][-1, :]
                xdata = self.data['x']
            elif var == 't':
                ydata = self.arrays[field_idx][:, len(self.data['x'])/2]
                xdata = self.data['t']
            ax.plot(xdata, ydata)

    def read(self):
        for i, name in enumerate(self.data['field_names']):
            fp = open(os.path.join(self.data['prj_name'], '//res//{}'.format(name)))
            self.arrays[i] = np.fromfile(fp)
            fp.close()


class Gui(tk.Tk):
    def __init__(self, f):
        tk.Tk.__init__(self)

        data = self.read_meta(f)
        self.data = data

        self.frame1 = tk.Frame(self)

        self.figures = FigFrame(self, data)

        self.frame3 = tk.Frame(self)

        self.frame1.pack(side=tk.LEFT)
        self.figures.pack(side=tk.LEFT)
        self.frame3.pack(side=tk.LEFT)

        self.box1_value = tk.StringVar()
        self.box1 = ttk.Combobox(self.frame1, textvariable=self.box1_value,
                                state='readonly')
        self.box2_value = tk.StringVar()
        self.box2 = ttk.Combobox(self.frame1, textvariable=self.box2_value,
                                 state='readonly')
        self.box3_value = tk.StringVar()
        self.box3 = ttk.Combobox(self.frame1, textvariable=self.box3_value,
                                 state='readonly')
        self.box4_value = tk.StringVar()
        self.box4 = ttk.Combobox(self.frame1, textvariable=self.box4_value,
                                 state='readonly')
        self.box1.grid(column=1, row=0)
        self.box2.grid(column=1, row=1)
        self.box3.grid(column=1, row=2)
        self.box4.grid(column=1, row=3)

        self.fig2 = plt.figure(figsize=(4, 10))

        fig_canvas_agg2 = FigureCanvasTkAgg(self.fig, master=self.frame2)
        fig_canvas_agg2.draw()
        fig_canvas_agg2.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        self.fig_canvas_agg2.mpl_connect('button_press_event', self.onclick)

        ax_ini_vals = self.fig2.add_subplot(111)
        draw_ini_vals(ax_ini_vals)

    def read_meta(self, f):
        d = dict()
        with open(f) as fp:
            d = json.load(fp)

        return d

    def onclick(self):
        #change material
        pass

    def update_figures(self):
        self.figures.draw_pics((self.box1.get(), self.box2.get(), self.box3.get(), self.box4.get()),
                               ('x', 'x', 'x', 'x'))

    def manage_draw_controls(self):
        values = self.data['field_names']
        if self.data["dim"] == 1:
            self.box1['values'] = values
            self.box1.current(0)
            self.box2['values'] = values
            self.box2.current(1)
            self.box3['values'] = values
            self.box3.current(2)
            self.box4['values'] = values
            self.box4.current(3)
        elif self.data["dim"] == 2:
            self.box1['values'] = values
            self.box1.current(0)
            self.box2['values'] = values
            self.box2.current(1)
            self.box3['values'] = values
            self.box3.current(2)
            self.box4['values'] = values
            self.box4.current(3)

        self.box1.bind("<<ComboboxSelected>>", self.update_figures)
        self.box3.bind("<<ComboboxSelected>>", self.update_figures)
        self.box4.bind("<<ComboboxSelected>>", self.update_figures)
        self.box2.bind("<<ComboboxSelected>>", self.update_figures)


if __name__ == "__main__":
    addr = raw_input("Project?")
    window = Gui(addr)

    window.protocol("WM_DELETE_WINDOW", lambda: window.quit())
    window.mainloop()
    window.destroy()
