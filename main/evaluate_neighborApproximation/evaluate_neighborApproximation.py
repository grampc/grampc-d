# This file is part of GRAMPC-D - (https://github.com/DanielBurk/GRAMPC-D.git)
#
# GRAMPC-D -- A software framework for distributed model predictive control (DMPC)
# based on the alternating direction method of multipliers (ADMM).
#
# Copyright 2020 by Daniel Burk, Andreas Voelz, Knut Graichen
# All rights reserved.
#
# GRAMPC-D is distributed under the BSD-3-Clause license, see LICENSE.txt

from operator import add
from pathlib import Path
import matplotlib.pyplot as plt
import sys, os, numpy

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

optimization_info.COMMON_Nhor_ = 16
optimization_info.COMMON_Thor_ = 4
optimization_info.COMMON_dt_ = 0.5
optimization_info.GRAMPC_MaxGradIter_ = 40
optimization_info.GRAMPC_MaxMultIter_ = 2
optimization_info.ADMM_innerIterations_ = 5
optimization_info.ADMM_maxIterations_ = 100
optimization_info.ADMM_ConvergenceTolerance_ = 0.009
optimization_info.ADMM_DebugCost_ = 1

approx = 1
optimization_info.APPROX_ApproximateCost_ = approx
optimization_info.APPROX_ApproximateConstraints_ = approx
optimization_info.APPROX_ApproximateDynamics_ = approx

interface.set_optimizationInfo(optimization_info)

Tsim = 0.01

# paremeters for cost function
P = 1
Q = 1
R = 0.1;

# parameters for model
A = 0.1
a = 0.005
d = 0.01

# initial and desired states and controls
xinit = [0.5]
uinit = [0]
xdes = [2] 
udes = [0];

# register agent
agent = grampcd_interface.AgentInfo()
agent.model_name_ = "water_tank_agentModel"

agent_id = 0
agent.id_ = agent_id
agent.model_parameters_ = [A, 1, 0]
agent.cost_parameters_ = [0, 0, R]

interface.register_agent(agent, xinit, uinit)
interface.set_desiredAgentState(agent_id, xdes, udes)

n_agents = 5

for i in range(1, n_agents-1):
    agent.id_ = i
    agent.model_parameters_ = [A, 0, 0]
    agent.cost_parameters_ = [0, 0, 0]
    interface.register_agent(agent, xinit, uinit)
    interface.set_desiredAgentState(i, xdes, udes)

agent_id = n_agents-1
agent.id_ = agent_id
agent.model_parameters_ = [A, 0, d]
agent.cost_parameters_ = [P, Q, 0]
interface.register_agent(agent, xinit, uinit)
interface.set_desiredAgentState(agent_id, xdes, udes)

# register couplings
coupling_info = grampcd_interface.CouplingInfo()
coupling_info.model_name_ = 'water_tank_couplingModel'
coupling_info.model_parameters_ = [A, a]

for i in range(0, n_agents):
    coupling_info.agent_id_ = i
    if(i > 0):
        coupling_info.neighbor_id_ = i-1
        interface.register_coupling(coupling_info)
    if(i < n_agents-1):
        coupling_info.neighbor_id_ = i+1
        interface.register_coupling(coupling_info)

# run distributed controller
interface.run_DMPC(0, Tsim)

interface.print_solution_to_file('all', 'with')
solution = interface.get_solution('all')

cost_with = [0]*len(solution[0].debug_cost_)
for i in range(0, n_agents):
    cost_with = list(map(add, cost_with, solution[i].debug_cost_));

numpy.savetxt("cost_with.csv", cost_with, delimiter=",")

# deregister agents
for i in range(0, n_agents):
    agent.id_ = i
    interface.deregister_agent(agent)

approx = 0
optimization_info.APPROX_ApproximateCost_ = approx
optimization_info.APPROX_ApproximateConstraints_ = approx
optimization_info.APPROX_ApproximateDynamics_ = approx

interface.set_optimizationInfo(optimization_info)

# register agents
agent_id = 0
agent.id_ = agent_id
agent.model_parameters_ = [A, 1, 0]
agent.cost_parameters_ = [0, 0, R]

interface.register_agent(agent, xinit, uinit)
interface.set_desiredAgentState(agent_id, xdes, udes)

for i in range(1, n_agents-1):
    agent.id_ = i
    agent.model_parameters_ = [A, 0, 0]
    agent.cost_parameters_ = [0, 0, 0]
    interface.register_agent(agent, xinit, uinit)
    interface.set_desiredAgentState(i, xdes, udes)

agent_id = n_agents-1
agent.id_ = agent_id
agent.model_parameters_ = [A, 0, d]
agent.cost_parameters_ = [P, Q, 0]
interface.register_agent(agent, xinit, uinit)
interface.set_desiredAgentState(agent_id, xdes, udes)

# register couplings
for i in range(0, n_agents):
    coupling_info.agent_id_ = i
    if(i > 0):
        coupling_info.neighbor_id_ = i-1
        interface.register_coupling(coupling_info)
    if(i < n_agents-1):
        coupling_info.neighbor_id_ = i+1
        interface.register_coupling(coupling_info)

# run distributed controller
interface.run_DMPC(0, Tsim)

solution = interface.get_solution('all')

fig=plt.figure()

cost_without = [0]*len(solution[0].debug_cost_)
for i in range(0, n_agents):
    cost_without = list(map(add, cost_without, solution[i].debug_cost_));
numpy.savetxt("cost_without.csv", cost_without, delimiter=",")
plt.plot(cost_with, 'b', label='With neighbor approximation')
plt.plot(cost_without, 'r--', label='Without neighbor approximation')
plt.legend(loc='lower right')
plt.title('Simulation example \n Evaluate neighbor approximation', y=1)
plt.ylabel('Global cost')
plt.xlabel('ADMM iterations')

plt.subplots_adjust(hspace=0.5)

plt.show()