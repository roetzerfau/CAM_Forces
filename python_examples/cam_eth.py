from __future__ import print_function
from datetime import datetime

import os, sys
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
import csv
import random

n_steps =100
def stencil_size(jump_parameter, area, nx):
  return min(5,np.ceil(jump_parameter/(area ** (1.0/len(nx)))))

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
  scaling = 2.5
  Nx = int(2500/scaling)
  
  print(Nx)

  Nx_file = 250
  nx_base =  [Nx_file, Nx_file]

  const                 = CAM.config()
  const.nx              = [Nx,Nx]
  const.debug_mode      = debug_mode
  #macro_names = ["DFACE_ATTRACTIVITY", "DROTATION", "DROTATION_COMPOSITES", "DSTENCIL_4_ALL_BUS", "DSUB_COMPOSITES"]
  const.ca_settings     = [True, False, False, False, False]
  
  jump_parameter_composites  = 20/scaling
  jump_parameter = 20/scaling
  aimPor = 0.35
  numCells = np.prod(const.nx)
  
  faces = [1] * 4
  PyCAM = CAM.include(const)
  Domain = PyCAM(jump_parameter_composites)

  texture = 'loam_bayreuth'
  print("Domain folded")
  
  #Domain.place_single_cell_bu_randomly(jump_parameter, aimPor , 0)
  
  
  mat = scipy.io.loadmat("particleShapeLibrary/particleShapes" + str(Nx_file) + "rotations.mat")
  particleList = mat['fullParticleList'][0]
  minFeretDiam = mat['fullParticleMinFeretDiameters']
  particleAreas = mat['fullParticleAreas'][0]

  particleInd_position = []
  intLowBound = []
  intUpBound = []  
  intSize = []
  with open('particleSizeDistribution/' + texture + '.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        line_count += 1
        intLowBound.append(float(row[0]))
        intUpBound.append(float(row[1]))
        intSize.append(float(row[2]))


  mat = scipy.io.loadmat("particleShapeLibrary/POMshapes" + str(Nx_file) + ".mat")
  POMparticleShapesList = mat['POMparticleShapesList'][0]
  POMparticleAreas = mat['POMparticleAreas'][0]
  POMminFeret = mat['POMminFeret'][0]
  POMsizeDistr = [[10, 0.2],[ 15, 0.3],[ 20 ,0.5]]
  amount = 0.005
  file_properties = 'domain_' + str(Nx)  + '_por_' + str(aimPor) + '_texture_' + texture + '_pxsize_' + str(scaling * 2)
  file_particleInd_position = 'MineralParticleInd_position_domain_' + file_properties + '.txt'
  file_POMParticleInd_position = 'POMParticleInd_position_domain_' + file_properties + '.txt'
  fileDomain = 'file_domain_' + file_properties+  '.txt'
  folder_saveData = 'saveData/'+ file_properties + '/'

  if not os.path.exists(folder_saveData):  os.makedirs(folder_saveData)



  if False:
    folder_saveData = 'saveData_test/'+ file_properties + '/'
    properties = [1, 1, scaling] 
    particle = particleList[10000]
    #particle = particleList[-100]
    #particle = particleList[0]
    #print(particle)
    #particle = [0,100,1,101,2,102]
    for i in range(100):
      particle = particleList[i]
      position = random.randrange(numCells)
      #position = 0
      stencil = stencil_size(jump_parameter, len(particle) ,const.nx)
      Domain.place_particle(stencil, particle,nx_base, position, faces, properties)#, 
      #Domain.place_single_cell_bu_randomly(jump_parameter, aimPor , 0)
  else: 
    newDomain = True
    if newDomain:
      numIntervals = len(intSize)

      interval = numIntervals - 1   

      particleCandidatesInt = [i for i in range(len(minFeretDiam)) if minFeretDiam[i] >= intLowBound[interval] and minFeretDiam[i] < intUpBound[interval]]
      fields = Domain.fields()
      positionCandidates = [i for i in range(len(fields)) if fields[i] == 0]

      currPor = Domain.porosity_d()

      maxUndershoot = 0.005
      
      while (currPor > aimPor):
        print("Porosity ", currPor)
      
        # ind max. particle Size s.t. interval is only undershot by maxUndershoot

        maxParticleSize = ((np.sum(intSize[interval:numIntervals]) + maxUndershoot) * (1-aimPor) * numCells -  (1-currPor) * numCells) * scaling**len(const.nx)
    
      # helper to find candidates that don't undershoot too much
        particleCandidatesHelper = [i for i in range(len(particleAreas)) if particleAreas[i] <= maxParticleSize]
        
        particleCandidatesInds = list(set(particleCandidatesHelper).intersection(set(particleCandidatesInt)))

        numCandidates = len(particleCandidatesInds)

        #continue if no candidate found
        if(numCandidates == 0 or ( np.sum(intSize[interval:numIntervals]) * (1-aimPor) < (1-currPor) )):
            interval = interval - 1

            print('Start placing particles of size >= ',intLowBound[interval], 'and <', intUpBound[interval])
            particleCandidatesInt = [i for i in range(len(minFeretDiam)) if minFeretDiam[i] >= intLowBound[interval] and minFeretDiam[i] < intUpBound[interval]]
            #fields = Domain.fields()
            #positionCandidates = [i for i in range(len(fields)) if fields[i] == 0]
        else:
          
          randInd =  random.randrange(numCandidates)
          candInd = particleCandidatesInds[randInd]
          particle = particleList[candInd]

          properties = [1, 1, scaling] 
          if (minFeretDiam[candInd] < 6.3):
            properties[1] = 1
          elif (minFeretDiam[candInd] < 20):
            properties[1] = 0.5
          elif (minFeretDiam[candInd] < 63):
            properties[1] = 0.25
          else:
            properties[1] = 0.1
          





          #position = random.randrange(numCells)
          randIndex = random.randrange(len(positionCandidates))
          position = positionCandidates[randIndex]
          stencil = stencil_size(jump_parameter, len(particle) ,const.nx)
          if( Domain.place_particle(stencil, particle,nx_base, position, faces, properties)):
            particleInd_position.append([candInd, position]) 
            fields = Domain.fields()
            positionCandidates = [i for i in range(len(fields)) if fields[i] == 0]
          #else:
            #positionCandidates[randIndex] = []
            
        currPor = Domain.porosity_d() 

      with open(folder_saveData +file_particleInd_position, 'w') as file:
        for item in particleInd_position:
          file.write (f"{item[0]}\t{item[1]}\n")
    else:
      with open(folder_saveData +file_particleInd_position, 'r') as file:
        for line in file:
          # Split the line into two items using the tab character as delimiter
          item = line.strip().split('\t')
          particleInd_position.append(item)
      for [particleInd, position] in particleInd_position:
        particleInd = int(particleInd)
        particle = particleList[particleInd]
        stencil = stencil_size(jump_parameter, len(particle) ,const.nx)
        properties = [1, 1, scaling] 
        if (minFeretDiam[particleInd] < 6.3):
          properties[1] = 1
        elif (minFeretDiam[particleInd] < 20):
          properties[1] = 0.5
        elif (minFeretDiam[particleInd] < 63):
          properties[1] = 0.25
        else:
            properties[1] = 0.1

        Domain.place_particle(stencil, particle,nx_base,int(position), faces, properties)



    POMParticleInd_position = []
    properties = [2, 0, scaling] 
    print("POMParticles")
    if newDomain:
      #--------------POM Particle----------------------
      numSolidCells = np.sum([1 for x in Domain.fields() if x > 0] )
      numPOMCells = round(amount*numSolidCells * scaling**len(const.nx))
      #print(np.shape(POMsizeDistr)[0])
      numPOMparticles = [0] * len(POMsizeDistr)
      for i  in range(len(numPOMparticles)):
          numPOMparticles[i] = ( round(numPOMCells * POMsizeDistr[i][1] / POMsizeDistr[i][0]))
      #print(numPOMparticles)
      for interval in range(len(numPOMparticles)-1, -1, -1):
        particleCandidatesInds = [i for i in range(len(POMparticleAreas)) if POMparticleAreas[i] == POMsizeDistr[interval][0]]
        for j in range(numPOMparticles[interval]):
          numCandidates = len(particleCandidatesInds)
          # choose candidate
          randInd =  random.randrange(numCandidates)
          candInd = particleCandidatesInds[randInd]
          particle = POMparticleShapesList[candInd]

          
          stencil = stencil_size(jump_parameter, len(particle) ,const.nx)

          newPositionFound = False
          while(not newPositionFound):
            position =  random.randrange(numCells)
            #candPos = [i for i in range(numCells) if Domain.fields()[i] == 0]
            #position =candPos[ random.randrange(len(position))]
            #print(position)
            if(Domain.place_particle(stencil, particle,nx_base,position, faces, properties)):
              newPositionFound = True
              POMParticleInd_position.append([candInd, position]) 
              currPor = Domain.porosity_d() 
              print("Porosity ", currPor)
          
            
        # print("Len ", len(POMParticleInd_position))
      with open(folder_saveData + file_POMParticleInd_position, 'w') as file:
        for item in POMParticleInd_position:
          file.write (f"{item[0]}\t{item[1]}\n")
    else:
      with open(folder_saveData +file_POMParticleInd_position, 'r') as file:
        for line in file:
          # Split the line into two items using the tab character as delimiter
          item = line.strip().split('\t')
          POMParticleInd_position.append(item)
        for [particleInd, position] in POMParticleInd_position:
          particle = POMparticleShapesList[int(particleInd)]
          stencil = stencil_size(jump_parameter, len(particle) ,const.nx)
          Domain.place_particle(stencil, particle,nx_base, int(position), faces, properties)

  save_data = np.zeros( (n_steps + 1, np.prod(const.nx)) ) 
  save_data_single = np.zeros( (2, np.prod(const.nx)) ) 
  save_data[0] = Domain.fields()
  for step in range(n_steps):
    Domain.do_cam()
    print("step ", step)
    save_data[step+1] = Domain.fields()


  end_time = datetime.now() 
  print("Program ended at", end_time, "after", end_time-start_time)
  #save_data[(save_data < 2500) & (save_data > 0)] = 1
  #save_data[save_data >= 2500] = 2
  

  with open(folder_saveData + fileDomain, 'w') as file:
    for data in save_data[-1]:
      file.write('%d \n' % data)
      #file.write(f"{data}\n")
  print("Plot: ", folder_saveData)

  save_data_single[0] = save_data[0]
  save_data_single[1] = save_data[-1]
  
  plot_to_vtk(folder_saveData+ "cam", save_data_single, const.nx)
  
  #plot_to_file(const.nx, save_data[-1], folder_saveData + 'camlast.png')
  #plot_to_file(const.nx, save_data[0], folder_saveData+ 'cam0.png')

  #plot(const.nx, save_data, 0)

  os.system("cp " + folder_saveData + " /mnt/c/users/maxro/Downloads/ -r")
# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):

  cam_test(n_steps, debug_mode)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")
