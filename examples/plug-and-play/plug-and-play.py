# This file is part of GRAMPC-D - (https://github.com/grampc-d/grampc-d.git)
#
# GRAMPC-D -- A software framework for distributed model predictive control (DMPC)
# based on the alternating direction method of multipliers (ADMM).
#
# Copyright 2020 by Daniel Burk, Andreas Voelz, Knut Graichen
# All rights reserved.
#
# GRAMPC-D is distributed under the BSD-3-Clause license, see LICENSE.txt

from pathlib import Path
import matplotlib.pyplot as plt
import sys, os, math, numpy
from operator import add

# generate path to module
path = os.getcwd()
sys.path.append(os.path.join(path, 'bin'))

path = str(Path(path).parents[1])
path = os.path.join(path, 'bin')

# append path to list of folders where python is searching for modules
sys.path.append(path)

#import python interface
import grampcd_interface

# create interface
interface = grampcd_interface.interface()

# initialize communication interface
interface.initialize_central_communicationInterface()

# set optimization info
optimization_info = grampcd_interface.OptimizationInfo()
optimization_info.COMMON_Nhor_ = 30
optimization_info.COMMON_Thor_ = 12
optimization_info.COMMON_dt_ = 0.1
optimization_info.GRAMPC_MaxGradIter_ = 40
optimization_info.GRAMPC_MaxMultIter_ = 1
optimization_info.ADMM_innerIterations_ = 2
optimization_info.ADMM_maxIterations_ = 20
optimization_info.ADMM_ConvergenceTolerance_ = 0.0005
interface.set_optimizationInfo(optimization_info)

Tsim = 20

# parameters for model
Omega = 1;
I = 1; 
P_max = 0.1;
kappa = 1e-3; 

# parameters for cost function
P = 0.1;
Q = 1;
R = 0.01;

# initial and desired states and controls
uinit = [0]
xinit = [0, 0] 
xdes = [0, 0] 
udes = [0];

# register agents
n_agents = 3
agent = grampcd_interface.AgentInfo()
agent.model_name_ = "smartGrid_agentModel"

agent_id = 0
P0 = 0
p = 1
agent.id_ = agent_id
agent.model_parameters_ = [I, Omega, kappa, P0, p]
agent.cost_parameters_ = [0, P, 0, Q, R]
interface.register_agent(agent, xinit, uinit, xdes, udes)

agent_id = 1
P0 = -0.01
p = 0
agent.id_ = agent_id
agent.model_parameters_ = [I, Omega, kappa, P0, p]
agent.cost_parameters_ = [0, P, 0, Q, R]
interface.register_agent(agent, xinit, uinit, xdes, udes)

# register couplings
coupling_info = grampcd_interface.CouplingInfo()

coupling_info.model_name_ = 'smartGrid_couplingModel'
coupling_info.model_parameters_ = [I, Omega, P_max]

coupling_info.agent_id_ = 0
coupling_info.neighbor_id_ = 1
interface.register_coupling(coupling_info)

coupling_info.agent_id_ = 1
coupling_info.neighbor_id_ = 0
interface.register_coupling(coupling_info)

# run distributed controller
interface.run_DMPC(0, Tsim)

# register agent
agent_id = 2
agent.id_ = agent_id
interface.register_agent(agent, xinit, uinit, xdes, udes)

# register couplings
coupling_info.agent_id_ = 1
coupling_info.neighbor_id_ = 2
interface.register_coupling(coupling_info)

coupling_info.agent_id_ = 2
coupling_info.neighbor_id_ = 1
interface.register_coupling(coupling_info)

# continue simulation
interface.run_DMPC(Tsim, Tsim)

solution = interface.get_solution('all')

idx = len(solution[0].cost_) - len(solution[2].cost_)
global_cost = list(map(add, solution[0].cost_, solution[1].cost_));
for i in range(0, len(solution[2].cost_)):
    global_cost[idx+i] = global_cost[idx+i] + solution[2].cost_[i]

diff_01 = [0]*len(solution[0].agentState_.x_)
for i in range(0, len(diff_01)):
    diff_01[i] =  solution[0].agentState_.x_[i] - solution[1].agentState_.x_[i] 

diff = len(solution[1].agentState_.x_) - len(solution[2].agentState_.x_)
diff_12 = [0]*len(solution[2].agentState_.x_)
for i in range(0, len(solution[2].agentState_.x_)):
    diff_12[i] =  solution[1].agentState_.x_[diff+i] - solution[2].agentState_.x_[i] 
   
fig, axes = plt.subplots(1, 2)
fig.suptitle('Simulation example \n smart grid', y=1)

axes[0].plot(solution[0].agentState_.t_, global_cost)
axes[0].set(ylabel = 'Global cost')
axes[0].set(xlabel = 'Simulation time')

axes[1].plot(solution[0].agentState_.t_, diff_01[0::2], label='Plant - Sink 1')
axes[1].plot(solution[2].agentState_.t_, diff_12[0::2], label='Sink 1 - Sink 2')
plt.legend(loc='upper right')
axes[1].set(ylabel = 'Phase shift')
axes[1].set(xlabel = 'Simulation time')

plt.subplots_adjust(wspace = 0.6)

plt.show()