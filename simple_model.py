#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 15:11:54 2020

@author: ryanstrickland
"""

from landlab.components import FlowAccumulator, FastscapeEroder, PerronNLDiffuse, LinearDiffuser
from landlab.plot import imshow_grid
from landlab import RasterModelGrid
from matplotlib.pyplot import figure
import matplotlib
import numpy as np

cell_dim = .5 #cell width/length in km
grid_width = 300 #number of nodes width
grid_height = 300 #number of nodes height

#Diffusion and stream power input parameters
uplift_rate = 0.01 #Units of km/yr
k = 0.0001 #Constant 0 < k < 1 generally
t_trans = 1000 #Number of years for transition from k_initial to k_final
trans_start_t = 25000 #Year the k transition begins
k_streampower = 0.3 
m_streampower = 0.5
total_t = 20000 #Total time 
dt =  100 #Timestep size
n_steps = total_t // dt #The number of steps in the model run


'''MODEL INSTANTIATION-DO NOT EDIT'''
#Instantiate Landlab RasterModelGrid
mg = RasterModelGrid((grid_width, grid_height), cell_dim) #Instantiate the model space


#Create the topographic__elevation array, fill with zeros
z = mg.add_zeros('node', 'topographic__elevation') #Base topographic elevation is zeros
initial_roughness = np.random.rand(z.size) / 100000.
z += initial_roughness

#Create the diffusivity array, fill with zeros
k = mg.zeros('node', dtype = float)
k += [k_initial for _ in range(mg.number_of_nodes)]
k_field = mg.add_field('node', 'diffusivity', k, noclobber = False)

            
#Establish boundary conditions
mg.set_fixed_value_boundaries_at_grid_edges(False, True, True, True)

'''SETUP THE COMPONENTS'''
#Setup the modules
linear_diffuser = LinearDiffuser(mg, linear_diffusivity = k)
flowRouter = FlowAccumulator(mg)
streamPower = FastscapeEroder(mg, K_sp = k_streampower, m_sp = m_streampower)


'''RUN THE MODEL'''
#Establish a steady state landscape under the initial conditions
for i in range(n_steps):
    z[mg.core_nodes] += uplift_rate * dt #add uplift
    linear_diffuser.run_one_step(dt)
    flowRouter.run_one_step() #run the flow router
    streamPower.run_one_step(dt) #run stream power incision
    
    if i % 100 == 0:
        print(i * dt, "years have elapsed")


figure()
imshow_grid(mg, 'topographic__elevation')



