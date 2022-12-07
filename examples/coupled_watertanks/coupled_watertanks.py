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
interface.initialize_central_communicationInterface()

# set optimization info
optimization_info = grampcd_interface.optimization_info()
optimization_info.COMMON_Nhor_ = 21
optimization_info.COMMON_Thor_ = 5
optimization_info.COMMON_dt_ = 0.1
optimization_info.GRAMPC_MaxGradIter_ = 10
optimization_info.GRAMPC_MaxMultIter_ = 2
optimization_info.ADMM_maxIterations_ = 10
optimization_info.ADMM_ConvergenceTolerance_ = 0.02

approx = 1
optimization_info.APPROX_ApproximateCost_ = approx
optimization_info.APPROX_ApproximateConstraints_ = approx
optimization_info.APPROX_ApproximateDynamics_ = approx

interface.set_optimizationInfo(optimization_info)

Tsim = 25

# parameters for cost function
P = 1
Q = 1
R = 0.1;

# parameters for model
A = 0.1
a = 0.005
d = 0.01

# inital and desired states and controls
xinit = [0.5]
uinit = [0]
xdes = [2] 
udes = [0];

# register agents
agent = grampcd_interface.agent_info()
agent.model_name_ = "water_tank_agentModel"

agent_id = 1
agent.id_ = agent_id
agent.model_parameters_ = [A, 1, 0]
agent.cost_parameters_ = [0, 0, R]
interface.register_agent(agent, xinit, uinit, xdes, udes)

agent_id = 2
agent.id_ = agent_id
agent.model_parameters_ = [A, 0, 0]
agent.cost_parameters_ = [0, 0, 0]
interface.register_agent(agent, xinit, uinit, xdes, udes)

agent_id = 3
agent.id_ = agent_id
agent.model_parameters_ = [A, 0, 0]
agent.cost_parameters_ = [0, 0, 0]
interface.register_agent(agent, xinit, uinit, xdes, udes)

agent_id = 4
agent.id_ = agent_id
agent.model_parameters_ = [A, 0, d]
agent.cost_parameters_ = [P, Q, 0]
interface.register_agent(agent, xinit, uinit, xdes, udes)

# register couplings
coupling_info = grampcd_interface.coupling_info()
coupling_info.model_name_ = 'water_tank_couplingModel'
coupling_info.model_parameters_ = [A, a]

coupling_info.agent_id_ = 1
coupling_info.neighbor_id_ = 2
interface.register_coupling(coupling_info)

coupling_info.agent_id_ = 1
coupling_info.neighbor_id_ = 3
interface.register_coupling(coupling_info)

coupling_info.agent_id_ = 2
coupling_info.neighbor_id_ = 1
interface.register_coupling(coupling_info)

coupling_info.agent_id_ = 2
coupling_info.neighbor_id_ = 4
interface.register_coupling(coupling_info)

coupling_info.agent_id_ = 3
coupling_info.neighbor_id_ = 1
interface.register_coupling(coupling_info)

coupling_info.agent_id_ = 3
coupling_info.neighbor_id_ = 4
interface.register_coupling(coupling_info)

coupling_info.agent_id_ = 4
coupling_info.neighbor_id_ = 2
interface.register_coupling(coupling_info)

coupling_info.agent_id_ = 4
coupling_info.neighbor_id_ = 3
interface.register_coupling(coupling_info)

# activate progressbar
interface.set_print_progressbar(True)

# run distributed controller
interface.run_DMPC(0, Tsim)

# get solutions
solution = interface.get_solution('all')

# plot solution
fig = plt.figure()
fig.suptitle('Simulation example \n Coupled water tanks', y=1)

for i in range(0, 4):
    plt.subplot(2, 2, i+1)
    plt.title('Agent ' + str(solution[i].agentState_.i_))
    plt.plot(solution[i].agentState_.t_, solution[i].agentState_.x_)
    plt.ylabel('Water height')
    plt.xlabel('Simulation time in s')

plt.subplots_adjust(hspace = 0.75, wspace = 0.5)
plt.show()