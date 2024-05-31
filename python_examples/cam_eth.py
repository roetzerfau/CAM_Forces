from __future__ import print_function
from datetime import datetime

import os, sys
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
import csv
import random

def stencil_size(jump_parameter, area, nx):
  return (jump_parameter/(area ** (1.0/len(nx))))

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
  print("Compiled")
  Nx = 500

  const                 = CAM.config()
  const.nx              = [Nx,Nx]
  const.debug_mode      = debug_mode
  #macro_names = ["DFACE_ATTRACTIVITY", "DROTATION", "DROTATION_COMPOSITES", "DSTENCIL_4_ALL_BUS", "DSUB_COMPOSITES"]
  #const.ca_settings     = [False, False, False, False, False]
  jump_parameter_composites  = 100
  jump_parameter = 10
  aimPor = 0.35
  faces = [1] * 4
  PyCAM = CAM.include(const)
  Domain = PyCAM(jump_parameter_composites)
  print("Domain folded")
  properties = [1] 
  #Domain.place_single_cell_bu_randomly(jump_parameter, aimPor , 0)
  
  
  mat = scipy.io.loadmat("particleShapeLibrary/particleShapes" + str(Nx) + "rotations.mat")
  particleList = mat['fullParticleList'][0]
  minFeretDiam = mat['fullParticleMinFeretDiameters']
  particleAreas = mat['fullParticleAreas'][0]

  particleInd_position = []
  intLowBound = []
  intUpBound = []  
  intSize = []
  with open('particleSizeDistribution/clay19_geoderma_waterStable_200.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        line_count += 1
        intLowBound.append(float(row[0]))
        intUpBound.append(float(row[1]))
        intSize.append(float(row[2]))
 

  newDomain = False
  if newDomain:
    numIntervals = len(intSize)

    interval = numIntervals - 1   

    particleCandidatesInt = [i for i in range(len(minFeretDiam)) if minFeretDiam[i] >= intLowBound[interval] and minFeretDiam[i] < intUpBound[interval]]

    data = Domain.fields()

    currPor = Domain.porosity_d()

    maxUndershoot = 0.005
    numCells = np.prod(const.nx)
    while (currPor > aimPor):
      print("Porosity ", currPor)
    
      # ind max. particle Size s.t. interval is only undershot by maxUndershoot

      maxParticleSize = (np.sum(intSize[interval:numIntervals]) + maxUndershoot) * (1-aimPor) * numCells -  (1-currPor) * numCells
  
    # helper to find candidates that don't undershoot too much
      particleCandidatesHelper = [i for i in range(len(particleAreas)) if particleAreas[i] <= maxParticleSize]
      
      particleCandidatesInds = list(set(particleCandidatesHelper).intersection(set(particleCandidatesInt)))

      numCandidates = len(particleCandidatesInds)

      #continue if no candidate found
      if(numCandidates == 0 or ( np.sum(intSize[interval:numIntervals]) * (1-aimPor) < (1-currPor) )):
          interval = interval - 1

          print('Start placing particles of size >= ',intLowBound[interval], 'and <', intUpBound[interval])
          particleCandidatesInt = [i for i in range(len(minFeretDiam)) if minFeretDiam[i] >= intLowBound[interval] and minFeretDiam[i] < intUpBound[interval]]
      else:
        
        randInd =  random.randrange(numCandidates)
        candInd = particleCandidatesInds[randInd]
        particle = particleList[candInd]

        position = random.randrange(numCells)
        stencil = stencil_size(jump_parameter, len(particle) ,const.nx)
        if( Domain.place_particle(stencil, particle, position, faces, properties)):
          particleInd_position.append([candInd, position]) 
          
      currPor = Domain.porosity_d() 

    with open('particleInd_domain.txt', 'w') as file:
      for item in particleInd_position:
        file.write (f"{item[0]}\t{item[1]}\n")
  else:
    with open('particleInd_domain.txt', 'r') as file:
      for line in file:
        # Split the line into two items using the tab character as delimiter
        item = line.strip().split('\t')
        particleInd_position.append(item)
    for [particleInd, position] in particleInd_position:
      particle = particleList[int(particleInd)]
      stencil = stencil_size(jump_parameter, len(particle) ,const.nx)
      Domain.place_particle(stencil, particle, int(position), faces, properties)


  #--------------POM Particle----------------------
  mat = scipy.io.loadmat("particleShapeLibrary/POMshapes" + str(Nx) + ".mat")
  POMparticleShapesList = mat['POMparticleShapesList'][0]
  POMparticleAreas = mat['POMparticleAreas']
  POMminFeret = mat['POMminFeret'][0]
  amount = 0.005
  numSolidCells = np.sum([1 for x in Domain.fields() if x > 0] )
  numPOMCells = round(amount*numSolidCells)
 # numPOMparticles = zeros(size(POMsizeDistr,1),1);
  #for i = 1 : length(numPOMparticles)
  #      numPOMparticles(i) = round(numPOMCells * POMsizeDistr(i,2) / POMsizeDistr(i,1));
 #   end  















  save_data = np.zeros( (n_steps + 1, np.prod(const.nx)) ) 
  save_data[0] = Domain.fields()
  print("hat")
  for step in range(n_steps):
    Domain.do_cam()
    print("tes1")
    save_data[step+1] = Domain.fields()


  end_time = datetime.now() 
  print("Program ended at", end_time, "after", end_time-start_time)
  save_data[(save_data < 2500) & (save_data > 0)] = 1
  save_data[save_data >= 2500] = 2
  if not os.path.exists('output'):  os.makedirs('output')
  plot_to_vtk("output/cam", save_data, const.nx)
  plot_to_file(const.nx, save_data[-1], 'output/camlast.png')
  plot_to_file(const.nx, save_data[0], 'output/cam0.png')
  #plot(const.nx, save_data, 0)

  os.system("cp output/ /mnt/c/users/maxro/Downloads/ -r")
# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  n_steps =100
  cam_test(n_steps, debug_mode)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")
