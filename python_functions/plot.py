import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
#matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from tkinter import *
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from pyevtk.hl import gridToVTK
from pyevtk.vtk import VtkGroup
import numpy as np


def plot(axes, save_data, time_step):
  if time_step < 0 or time_step >= len(save_data):
    print("Time step for plotting does not exist!")
    return

  if np.size(axes) == 1:  axes.append(1)

  root = Tk()
  animated_cam(root, axes, save_data, time_step)
  root.mainloop()

def plot_to_file(axes, save_data, file_name, text = 'boxed italics text in data coords', ax=plt ):
  if np.size(axes) == 1:  axes.append(1)
  fig, ax = config_plot(axes)
  # ax.set_title(text)
  plot_update(axes, save_data, ax)
  fig.savefig(file_name, dpi=400)
# save to vtk file for analysing in paraview 
def plot_to_vtk(filename, data, shape):
  n_steps = np.shape(data)[0]
  g = VtkGroup("./group")
  # Dimensions
  if len(shape) == 3:
    nx, ny, nz = shape[0], shape[1], shape[2]
  if len(shape) == 2:
    nx, ny, nz = shape[0], shape[1], 1
  lx, ly, lz = 1.0, 1.0, 1.0
  dx, dy, dz = lx / nx, ly / ny, lz / nz

  # Coordinates
  X = np.arange(0, lx + 0.1 * dx, dx, dtype="float64")
  Y = np.arange(0, ly + 0.1 * dy, dy, dtype="float64")
  Z = np.arange(0, lz + 0.1 * dz, dz, dtype="float64")
  
  
  for step in range(n_steps):
    file = filename + str(step)
    g.addFile(filepath=file, sim_time=step)
    gridToVTK(
    file,
    X,
    Y,
    Z,
    cellData={"cells": data[step]},
    )
  g.save()


def plot_update(axes, data, ax=plt):
  data = np.reshape(data, axes)
  if np.size(axes) == 1:  axes.append(1)

  if np.size(axes) == 2:
    data = (data != 0)
    # cmap = plt.get_cmap('hot')
    cmap = colors.ListedColormap(['white', 'k'])
    # cmap = colors.ListedColormap(['white' ,'sandybrown','darkgray'])
    ax.pcolor(data[::-1], cmap=cmap,edgecolors='b', linewidths=0)#
  elif np.size(axes) == 3:
    # data = (data != 0)
    Colors = np.empty(axes + [4], dtype=np.float32)
    # Control Transparency
    alpha = .9
    Colors[data == 1] = [0.957, 0.643, 0.376, alpha]
    Colors[data == 2] = [0.66, 0.66, 0.66, alpha]
    # Colors[data] = [0, 0, 0, alpha]
    ax.voxels(data, facecolors=Colors, edgecolors='black')
  return ax

def config_plot(axes):
  dim = np.size(axes)
  if dim == 1:
    dim = 2

  if dim == 2:    
    fig, ax = plt.subplots()

  elif dim == 3:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
  return (fig, ax)

class animated_cam:
  def __init__(self, master, axes, save_data, time_step):
    self.master = master
    self.frame = Frame(self.master)
    self.fig, self.ax = config_plot(axes)
    self.canvas = FigureCanvasTkAgg(self.fig, self.master)  
    self.config_window()
    self.axes      = axes
    self.save_data = save_data
    self.time_step = time_step
    self.n_steps   = len(save_data)
    self.draw_cam()
    self.frame.pack(expand=YES, fill=BOTH)

  def config_window(self):
    toolbar = NavigationToolbar2Tk(self.canvas, self.master)
    toolbar.update()
    self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
    self.button = Button(self.master, text="Quit", command=self._quit)
    self.button.pack(side=BOTTOM)
    self.button_next = Button(self.master, text="Next", command=self.plot_next)
    self.button_next.pack(side=RIGHT)
    self.button_prev = Button(self.master, text="Previous", command=self.plot_previous)
    self.button_prev.pack(side=LEFT)
    self.fig.canvas.mpl_connect('close_event', self.on_close)

  def draw_cam(self):
    self.ax.clear() # clear current axes
    data = self.save_data[self.time_step]
    self.ax = plot_update(self.axes, data, self.ax)
    self.ax.set(title='Timestep: ' + str(self.time_step))
    self.ax.set_aspect('equal', 'box')
    self.ax.axis('off')
    self.canvas.draw()

  def _quit(self):
    self.master.quit()  # stops mainloop

  def plot_next(self):
    self.time_step = (self.time_step + 1 ) % self.n_steps
    self.draw_cam()

  def plot_previous(self):
    self.time_step = (self.time_step - 1 + self.n_steps ) % self.n_steps
    self.draw_cam()

  def on_close(self, event):
    self.master.quit()
