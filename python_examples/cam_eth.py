from __future__ import print_function
from datetime import datetime

import os, sys
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
import csv


# --------------------------------------------------------------------------------------------------
# Function bilaplacian_test.
# --------------------------------------------------------------------------------------------------
def cam_test(n_steps, debug_mode=False):
  start_time = datetime.now()
  print("Starting time is", start_time)

  try:
    import CAM
  except (ImportError, ModuleNotFoundError) as error:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep  + ".." + os.sep + "import")
    import CAM

  try:
    from plot import plot, plot_to_file, plot_to_vtk
  except (ImportError, ModuleNotFoundError) as error:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep  + ".." + os.sep + 
      "python_functions")
    from plot import plot, plot_to_file, plot_to_vtk

  
  const                 = CAM.config()
  const.nx              = [500,500]
  const.debug_mode      = debug_mode
  jump_parameter_composites  = 5
  jump_parameter = 5
  porosity = 0.75
  faces = [1] * 4
  PyCAM = CAM.include(const)
  Domain = PyCAM(jump_parameter_composites)
  #Domain.place_single_cell_bu_randomly(jump_parameter, porosity , 0)
  
  Nx =500
  mat = scipy.io.loadmat("particleShapeLibrary/particleShapes" + str(Nx) + "rotations.mat")
  particleList = mat['fullParticleList']
  minFeretDiam = mat['fullParticleMinFeretDiameters']
  print(minFeretDiam[0])
  particle = particleList[0][0]

  Domain.place_particle(jump_parameter, particle, -1, faces)
  
  intLowBound = []
  intUpBound = []  
  intSize = []
  with open('particleSizeDistribution/clay34_geoderma_waterStable_200.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        line_count += 1
        intLowBound.append(float(row[0]))
        intUpBound.append(float(row[1]))
        intSize.append(float(row[2]))

  interval = len(intSize) - 1   
  print(intLowBound[interval] )
  a = [i for i,x in enumerate(minFeretDiam) if x >= intLowBound[interval] and x < intUpBound[interval]]#



  save_data = np.zeros( (n_steps + 1, np.prod(const.nx)) ) 
  save_data[0] = Domain.fields()

  for step in range(n_steps):
    Domain.do_cam()
    save_data[step+1] = Domain.fields()


  end_time = datetime.now() 
  print("Program ended at", end_time, "after", end_time-start_time)
  save_data[(save_data < 2500) & (save_data > 0)] = 1
  save_data[save_data >= 2500] = 2
  if not os.path.exists('output'):  os.makedirs('output')
  #plot_to_vtk("output/cam", save_data, const.nx)
  plot_to_file(const.nx, save_data[-1], 'output/cam.png')
  #plot(const.nx, save_data, 0)

  os.system("cp output/ /mnt/c/users/maxro/Downloads/ -r")
# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  n_steps =10
  cam_test(n_steps, debug_mode)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")
