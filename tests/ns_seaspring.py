from __future__ import division
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'model')))
import numpy as np
from saving import *
from make_boot_geometry import *
from solver import *

#Constant to use in later viscosity coefficient definition
seconds_per_year = 3.154E7
#Starting directory
main_directory = os.getcwd()
#Base name for the set of simulations
base_name = "ns_seaspring"
#Number of simulations
num_sims = 10
#Mesh resolution (meters)
hires = 10
lores = 40
#Generate geometry for use in meshing
xz_boundary = make_boot_geometry(10000, 400, 50, 20, hires, lores)
#Define the array of time steps to use for the set of simulations
sim_dts = [np.power(10.0,j)/(365*24*60*60) for j in range(0, int(num_sims))]
#Generate the set of folder names
sim_names = []
for j in range(len(sim_dts)):
    sim_names.append(base_name + "_dt_" + str(j))
#Make an empty list to use for solvers
Solvers = [None]*len(sim_dts)

for j in range(len(Solvers)):
    #Move to main directory and create each simulation
    os.chdir(main_directory)
    setup_folder(sim_names[j])
    #Define solver and run model for one time step.
    Solvers[j] = SimulationSolver(xz_boundary, 1E-9, 3.5E-25*seconds_per_year, 9800, lores, hires, 4000, True) #geometry, picard solver relative tolerance, glen's flow law parameter, high resolution start distance, low resolution, high resolution, inflow velocity (m/s), include_inertial_term
    Solvers[j].run(0.0, sim_dts[j])