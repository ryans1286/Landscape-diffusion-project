# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 15:07:48 2020

@author: rstrickland
"""

from landlab.components import FlowAccumulator, FastscapeEroder, PerronNLDiffuse
from landlab.plot import imshow_grid
from landlab import RasterModelGrid, CLOSED_BOUNDARY, FIXED_VALUE_BOUNDARY, load_params
from matplotlib.pyplot import figure, show, plot, xlabel, ylabel, title
import numpy as np
import random
from landlab.io.netcdf import write_netcdf


#Model geometry
cell_dim = 2. #cell width/length in km
grid_width = 100 #number of nodes width
grid_height = 100 #number of nodes height

#User input parameters
uplift_rate = 0.001 #Units of m/yr
k_initial = 0.0002 #Constant 0 < k < 1 generally
k_final = 0.002
total_t = 200000 #Total time 
dt =  100 #Timestep size
n_steps = total_t // dt #The number of steps in the model run

#Instantiate Landlab RasterModelGrid
mg = RasterModelGrid((grid_width, grid_height), cell_dim) #Instantiate the model space

#Establish the boundary conditions
for edge in (mg.nodes_at_left_edge, mg.nodes_at_right_edge):
    mg.status_at_node[edge] = FIXED_VALUE_BOUNDARY
for edge in (mg.nodes_at_top_edge, mg.nodes_at_bottom_edge):
    mg.status_at_node[edge] = FIXED_VALUE_BOUNDARY

#Create the topographic__elevation array, fill with zeros
z = mg.add_zeros('node', 'topographic__elevation') #Base topographic elevation is zeros
initial_roughness = np.random.rand(z.size) / 100000.
z += initial_roughness

#Create the diffusivity array, fill with zeros
k = mg.zeros('node', dtype = float)
k += [k_initial for _ in range(mg.number_of_nodes)]
k_field = mg.add_field('node', 'diffusivity', k, noclobber = False)

nonlinear_diffuser = PerronNLDiffuse(mg, nonlinear_diffusivity = k)
flowRouter = FlowAccumulator(mg)
streamPower = FastscapeEroder(mg, K_sp = 0.3, m_sp = 0.5)

for i in range(n_steps):
    z[mg.core_nodes] += uplift_rate * dt #add uplift
    nonlinear_diffuser.run_one_step(dt) #run diffusion one step
    flowRouter.run_one_step() #run the flow router
    streamPower.run_one_step(dt) #run stream power incision
    
    if i % 50 == 0:
        print(i * dt, "years have elapsed")

print(total_t, "years have elapsed. First round complete.")

figure('topo with diffusion')
write_netcdf('projectTest/run1_hi-lo-1.nc', mg, format = 'NETCDF3_64BIT', names = 'topographic__elevation')


#Redefine diffusion parameter
k = mg.zeros('node', dtype = float)
k += [k_final for _ in range(mg.number_of_nodes)]

nonlinear_diffuser = PerronNLDiffuse(mg, nonlinear_diffusivity = k)


for i in range(n_steps):
    z[mg.core_nodes] += uplift_rate * dt #add uplift
    nonlinear_diffuser.run_one_step(dt) #run diffusion one step
    flowRouter.run_one_step() #run the flow router
    streamPower.run_one_step(dt) #run stream power incision
    
    if i % 50 == 0:
        print(total_t + i * dt, "years have elapsed")
print(total_t, "years have elapsed. Model run complete.")

figure('topo with diffusion')
write_netcdf('projectTest/run1_hi-lo-2.nc', mg, format = 'NETCDF3_64BIT', names = 'topographic__elevation')

