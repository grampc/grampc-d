# This file is part of GRAMPC-D - (https://github.com/grampc-d/grampc-d.git)
#
# GRAMPC-D -- A software framework for distributed model predictive control (DMPC)
# 
#
# Copyright 2023 by Daniel Burk, Maximilian Pierer von Esch, Andreas Voelz, Knut Graichen
# All rights reserved.
#
# GRAMPC-D is distributed under the BSD-3-Clause license, see LICENSE.txt

from pathlib import Path
import matplotlib.pyplot as plt
import sys, os
import numpy as np

# generate path to module
path = os.getcwd()
sys.path.append(os.path.join(path, 'bin'))

path = str(Path(path).parents[1])
path = os.path.join(path, 'bin')

# append path to list of folders where python is searching for modules
sys.path.append(path)

# import python interface
import grampcd_interface

# create interface
interface = grampcd_interface.interface()

# initialize communication interface
interface.initialize_local_communicationInterface_as_coordinator(7777)

# set optimization info
optimization_info = grampcd_interface.optimization_info()
optimization_info.COMMON_Nhor_ = 30
optimization_info.COMMON_Thor_ = 1.0
optimization_info.COMMON_dt_ = 0.01
optimization_info.GRAMPC_MaxGradIter_ = 5
optimization_info.GRAMPC_MaxMultIter_ = 1

# ADMM does not converge monotonically if too few gradient iterations are used
optimization_info.COMMON_Solver_ = "ADMM"
optimization_info.ADMM_maxIterations_ = 5
optimization_info.ADMM_AdaptPenaltyParameter_ = False

interface.set_optimizationInfo(optimization_info)

Tsim = 6.0
num_agents = 3
num_couplings = 2

# wait for agents to connect
interface.wait_for_connections(num_agents, num_couplings)

# run distributed controller
interface.run_DMPC(0, Tsim)

# get solutions
solution = interface.get_solution('all')

# terminate agents
interface.send_flag_to_agents('all')

# plot solution
Nx = [4, 3, 3]
Nu = [3, 2, 2]

fig, axs = plt.subplots(num_agents, 3)
t = solution[0].agentState_.t_
Jsum = np.zeros((len(t)))
for i in range(0, num_agents):
    t = solution[i].agentState_.t_
    x = np.reshape(solution[i].agentState_.x_, (len(solution[i].cost_), Nx[i]))
    u = np.reshape(solution[i].agentState_.u_, (len(solution[i].cost_), Nu[i]))
    J = np.array(solution[i].cost_)
    Jsum = Jsum + J

    axs[i,0].plot(t, x)
    axs[i,0].set_title('Agent ' + str(solution[i].agentState_.i_))
    axs[i,0].set_xlabel('Simulation time')
    axs[i,0].set_ylabel('States')

    axs[i,1].plot(t, u)
    axs[i,1].set_title('Agent ' + str(solution[i].agentState_.i_))
    axs[i,1].set_xlabel('Simulation time')
    axs[i,1].set_ylabel('Controls')

    axs[i,2].plot(t, J)
    axs[i,2].set_title('Agent ' + str(solution[i].agentState_.i_))
    axs[i,2].set_xlabel('Simulation time')
    axs[i,2].set_ylabel('Cost')

plt.show()

fig, axs = plt.subplots(1, 2)
for i in range(0, num_agents):
    x = np.reshape(solution[i].agentState_.x_, (len(solution[i].cost_), Nx[i]))
    axs[0].plot(x[:, 0], x[:, 1])
axs[0].add_patch(plt.Circle((0, 0), 1, fill=False, linestyle='--'))
axs[0].axis('equal')
axs[0].set_xlabel('Position x')
axs[0].set_ylabel('Position y')

axs[1].plot(t, Jsum)
axs[1].set_xlabel('Simulation time')
axs[1].set_ylabel('Total cost')

plt.show()