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
optimization_info.COMMON_Thor_ = 1
optimization_info.COMMON_dt_ = 0.1
optimization_info.GRAMPC_MaxGradIter_ = 10
optimization_info.GRAMPC_MaxMultIter_ = 1
optimization_info.ADMM_maxIterations_ = 10
optimization_info.ADMM_ConvergenceTolerance_ = 0

optimization_info.COMMON_Debug_cost_ = 0
optimization_info.ASYNC_Delay_ = 0


interface.set_optimizationInfo(optimization_info)

Tsim = 3

# parameters for cost function
P = 1; Q = 1; R = 0.1;

# initial and desired states and controls
uinit = [0];
xdes = [0, 0]; 
udes = [0];
model_parameters = [1, 1, 1];
cost_parameters = [0, P, 0, Q, R];

# register agent
agent_id = 0;
xinit = [0.5, 0];
agent = grampcd_interface.agent_info()
agent.id_ = agent_id;
agent.model_name_ = "vdp_agentModel";
agent.model_parameters_ = model_parameters
agent.cost_parameters_ = cost_parameters
interface.register_agent(agent, xinit, uinit, xdes, udes);

agent_id = 1;
xinit = [0.25, 0];
agent = grampcd_interface.agent_info()
agent.id_ = agent_id;
agent.model_name_ = "vdp_agentModel";
agent.model_parameters_ = model_parameters
agent.cost_parameters_ = cost_parameters
interface.register_agent(agent, xinit, uinit, xdes, udes);

agent_id = 2;
xinit = [0.75, 0];
agent = grampcd_interface.agent_info()
agent.id_ = agent_id;
agent.model_name_ = "vdp_agentModel";
agent.model_parameters_ = model_parameters
agent.cost_parameters_ = cost_parameters
interface.register_agent(agent, xinit, uinit, xdes, udes);

# register coupling
coupling_info = grampcd_interface.coupling_info();
coupling_info.model_name_ = 'vdp_synchronize_couplingModel';
coupling_info.model_parameters_ = [1];
coupling_info.cost_parameters_ = [1];

coupling_info.agent_id_ = 0;
coupling_info.neighbor_id_ = 1;
interface.register_coupling(coupling_info);

coupling_info.agent_id_ = 1;
coupling_info.neighbor_id_ = 0;
interface.register_coupling(coupling_info);

coupling_info.agent_id_ = 1;
coupling_info.neighbor_id_ = 2;
interface.register_coupling(coupling_info);

coupling_info.agent_id_ = 2;
coupling_info.neighbor_id_ = 1;
interface.register_coupling(coupling_info);

# run distributed controller
interface.run_DMPC(0, Tsim)

# get solutions
solution = interface.get_solution('all')

# plot solution
fig = plt.figure()
fig.suptitle('Simulation example \n Coupled cost functions', y=1)

for i in range(0, 3):
    plt.plot(solution[i].agentState_.t_, solution[i].agentState_.x_[::2])
    plt.ylabel('Phase angle')
    plt.xlabel('Simulation time [s]')

plt.show()