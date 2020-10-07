# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 15:07:48 2020

@author: rstrickland
"""

from landlab.components import FlowAccumulator, FastscapeEroder, PerronNLDiffuse
from landlab.plot import imshow_grid
from landlab import RasterModelGrid#, load_params
from landlab.ca.celllab_cts import Transition, CAPlotter
from landlab.ca.raster_cts import RasterCTS
from matplotlib.pyplot import figure#, show, plot, xlabel, ylabel, title
import matplotlib
#import matplotlib.pyplot as plt
#from landlab.io.netcdf import write_netcdf
import numpy as np




'''FUNCTIONS-DO NOT EDIT'''
#Define the cellular automaton transition function
def transitions():
    xn_list = []
    
    #Transition((initial-state, initial-state, orientation), (final-state, final-state, orientation), rate, 'name')
    #Note the rate is transitions per second. 
    xn_list.append(Transition((0, 1, 0), (1, 1, 0), .5, 'left-right transition'))
    xn_list.append(Transition((1, 0, 0), (1, 1, 0), .5, 'left-right transition'))
    #These transitions are the most basic and symmetrical.
    #Once a node is "seeded" with state "1", it will gradually change all the  
    #bordering nodes to state "1"
    
    return xn_list




'''INPUT PARAMETERS'''
#Model geometry parameters
cell_dim = 1. #cell width/length in km
grid_width = 100 #number of nodes width
grid_height = 100 #number of nodes height

#Diffusion and stream power input parameters
uplift_rate = 0.001 #Units of m/yr
k_initial = 0.0001 #Constant 0 < k < 1 generally
k_final = 0.001
k_streampower = 0.3
m_streampower = 0.5
total_t = 60000 #Total time 
dt =  100 #Timestep size
n_steps = total_t // dt #The number of steps in the model run

#Cellular automaton parameters
plot_interval = 10. #needed if plotting cellular automaton transitions
run_duration = 100. #The number of seconds the cellular automaton model should run per timestep dt
#report_interval = plot_interval #needed if plotting node state transitions

#Create a csv file with the model parameters, if desired
#probably turn this into a dictionary later on
# parameter_names = ['cell_dim', 'grid_width', 'grid_height', 'uplift_rate', 
#                    'k_initial', 'k_final', 'k_streampower', 'm_streampower',
#                    'total_t', 'dt', 'n_steps', 'plot_interval', 'run_duration']
# parameter_values = [cell_dim, grid_width, grid_height, uplift_rate, 
#                    k_initial, k_final, k_streampower, m_streampower,
#                    total_t, dt, n_steps, plot_interval, run_duration]

# import pandas as pd
# fname = #csv filepath and name 
# df = pd.DataFrame([parameter_names, parameter_values], columns = ['Parameter', 'Value'])
# df.to_csv(fname, index = False)



'''MODEL INSTANTIATION-DO NOT EDIT'''
#Instantiate Landlab RasterModelGrid
mg = RasterModelGrid((grid_width, grid_height), cell_dim) #Instantiate the model space


# for edge in (mg.nodes_at_left_edge, mg.nodes_at_right_edge):
#     mg.status_at_node[edge] = FIXED_VALUE_BOUNDARY
# for edge in (mg.nodes_at_top_edge, mg.nodes_at_bottom_edge):
#     mg.status_at_node[edge] = FIXED_VALUE_BOUNDARY


#Create the topographic__elevation array, fill with zeros
z = mg.add_zeros('node', 'topographic__elevation') #Base topographic elevation is zeros
initial_roughness = np.random.rand(z.size) / 100000.
z += initial_roughness

#Create the diffusivity array, fill with zeros
k = mg.zeros('node', dtype = float)
k += [k_initial for _ in range(mg.number_of_nodes)]
k_field = mg.add_field('node', 'diffusivity', k, noclobber = False)

#Create the node state array for the cellular automaton component
#Create the node-state dictionary
ns_dict = {0 : 'k_i', 1 : 'k_f'}

#Setup cellular automaton plot colors
#This is only needed if you intend to plot the cellular automaton grids.
k_i = 'skyblue'
k_f = 'red'
colorList = [k_i, k_f]
my_cmap = matplotlib.colors.ListedColormap(colorList)

xn_list = transitions() #Initiate the transitions
node_state_grid = mg.zeros('node', dtype = int)
node_state_field = mg.add_field('node', 'k_diff_transitions', node_state_grid)

#Create initial conditions for cellular automaton component
#This is not working as intended, but does the job
for x in range(100):
    xRand = np.random.randint(1, grid_width)
    yRand = np.random.randint(1, grid_height)
    randSeed = np.where((mg.x_of_node == xRand) & (mg.y_of_node == yRand))[0]
    node_state_grid[randSeed] = 1

#Needed for plotting cellular automaton transitions
#current_real_time = time.time()
#next_report = current_real_time + report_interval
            
#Establish boundary conditions
mg.set_fixed_value_boundaries_at_grid_edges(True, True, True, True)




'''SETUP THE COMPONENTS'''
#Setup the modules
nonlinear_diffuser = PerronNLDiffuse(mg, nonlinear_diffusivity = k)
flowRouter = FlowAccumulator(mg)
streamPower = FastscapeEroder(mg, K_sp = k_streampower, m_sp = m_streampower)
ca_diffusion_transition = RasterCTS(mg, ns_dict, xn_list, node_state_grid)

#Run the cellular automaton model
ca_plotter = CAPlotter(ca_diffusion_transition, cmap = my_cmap)




'''RUN THE MODEL'''
#Establish a steady state landscape under the initial conditions
for i in range(n_steps):
    z[mg.core_nodes] += uplift_rate * dt #add uplift
    nonlinear_diffuser.run_one_step(dt) #run diffusion one step
    flowRouter.run_one_step() #run the flow router
    streamPower.run_one_step(dt) #run stream power incision
    
    if i % 50 == 0:
        print(i * dt, "years have elapsed")

print(total_t, "years have elapsed. First round complete.")

figure()
imshow_grid(mg, 'topographic__elevation')


#Initiate the cellular automaton model
#Note- when running in Spyder, both figures appear at the end
# for i in range(n_steps):
#     z[mg.core_nodes] += uplift_rate * dt #add uplift
#     nonlinear_diffuser.run_one_step(dt) #run diffusion one step
#     flowRouter.run_one_step() #run the flow router
#     streamPower.run_one_step(dt) #run stream power incision
    
#     current_time = 0. #run the cellular automaton for run_duration each timestep i
#     while current_time < run_duration:
#         ca_diffusion_transition.run(current_time + plot_interval, ca_diffusion_transition.node_state, plot_each_transition = False)
#         current_time += plot_interval
#     #This block is really slowing me down!     
#     for j in range(mg.number_of_nodes): #change diffusion values
#         if node_state_grid[j] == 1:
#             k[j] = k_final
#     if i % 50 == 0:
#         print(i * dt, "years have elapsed")

# figure()
# imshow_grid(mg, 'topographic__elevation')
#write_netcdf('run1_hi-lo-1.nc', mg, format = 'NETCDF3_64BIT', names = 'topographic__elevation')




'''
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
#write_netcdf('projectTest/run1_hi-lo-2.nc', mg, format = 'NETCDF3_64BIT', names = 'topographic__elevation')
'''
