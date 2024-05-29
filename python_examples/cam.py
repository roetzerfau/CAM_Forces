from __future__ import print_function
from datetime import datetime

import os, sys
import numpy as np


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
  PyCAM = CAM.include(const)
  Domain = PyCAM(jump_parameter_composites)

  Domain.place_single_cell_bu_randomly(jump_parameter, porosity , 0)

  # success = 0
  # while success < 10:
  #   success = success + Domain.place_sphere(jump_parameter, 5, -1)
  # print("Nr of spheres " + str(success))
  # success = 0
  # while success < 2500:
  #   success = success + Domain.place_plane(jump_parameter, [2,2,30], -1, 1)
  # print("Nr of planes " + str(success))
  # success = 0
  # while success < 5000:
  #   success = success + Domain.place_plane(jump_parameter, [2,2,6], -1, -1)
  # print("Nr of planes " + str(success))
  
  # Domain.print_array()

  # Domain.place_particles()
  
  save_data = np.zeros( (n_steps + 1, np.prod(const.nx)) ) 
  save_data[0] = Domain.fields()

  for step in range(n_steps):
    Domain.do_cam()
    save_data[step+1] = Domain.fields()
  
  # print(Domain.average_particle_size_d())
  # print(Domain.average_particle_size(save_data[-1]))
  # print(Domain.bulk_distance(save_data[0],save_data[-1]))
  end_time = datetime.now() 
  print("Program ended at", end_time, "after", end_time-start_time)
  save_data[(save_data < 2500) & (save_data > 0)] = 1
  save_data[save_data >= 2500] = 2
  if not os.path.exists('output'):  os.makedirs('output')
  #plot_to_vtk("output/cam", save_data, const.nx)
  plot_to_file(const.nx, save_data[-1], 'output/cam.png')
  #plot(const.nx, save_data, 0)
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
