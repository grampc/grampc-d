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
import sys, os
import time, numpy
from random import seed, random

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
optimization_info.COMMON_Nhor_ = 21
optimization_info.COMMON_Thor_ = 2
optimization_info.COMMON_dt_ = 0.02
optimization_info.GRAMPC_MaxGradIter_ = 15
optimization_info.GRAMPC_MaxMultIter_ = 1
optimization_info.ADMM_maxIterations_ = 10;
optimization_info.ADMM_ConvergenceTolerance_ = 0.002;

interface.set_optimizationInfo(optimization_info)

Tsim = 4

# parameters for cost function
P_x = 1; P_vx = 1; P_y = 1; P_vy = 1;
Q_x = 5; Q_vx = 2; Q_y = 5; Q_vy = 2;
R_ux = 0.01; R_uy = 0.01;

# model parameters
m_agent = 7.5; c = 0.5;

# number of agents
n_agents_x = 5;
n_agents_y = 5;
n_agents = n_agents_x * n_agents_y;

# register agents
agentInfo = grampcd_interface.AgentInfo()
agentInfo.model_name_ = "ssms2d_agentModel"
agentInfo.model_parameters_ = [m_agent, 1]
agentInfo.cost_parameters_ = [P_x, P_vx, P_y, P_vy, Q_x, Q_vx, Q_y, Q_vy, R_ux, R_uy]

x_init = [[0 for x in range(4)] for y in range(n_agents)];
for i in range(0, n_agents_x) :
    for j in range(0, n_agents_y) :
        id = i*n_agents_x + j;

        seed(time.time())
        initial = [i + random(), 0, j + 0.4*random(), 0];

        for k in range(0, 4) :
            x_init[id][k] = initial[k];

for i in range(0, n_agents_x) :
    for j in range(0, n_agents_y) :
        agentInfo.id_ = i*n_agents_x + j;

        interface.register_agent(agentInfo, x_init[agentInfo.id_], [0, 0]);

        x_des = [i, 0, j, 0];
        interface.set_desiredAgentState(agentInfo.id_, x_des, [0, 0]);

# register couplings
coupling_info = grampcd_interface.CouplingInfo()
coupling_info.model_name_ = 'ssms2d_couplingModel'
coupling_info.model_parameters_ = [m_agent, c]

idx = 0
for i in range(0, n_agents_y) :
    for j in range (0, n_agents_x) :

        coupling_info.agent_id_ = idx;

        # coupling with neighbor on the left
        if j > 0 :
            coupling_info.neighbor_id_ = idx - 1;
            interface.register_coupling(coupling_info);

        # coupling with neighbor on the right
        if j < n_agents_x - 1 :
            coupling_info.neighbor_id_ = idx + 1;
            interface.register_coupling(coupling_info);

        # coupling with neighbor above
        if i > 0 :
            coupling_info.neighbor_id_ = idx - n_agents_x;
            interface.register_coupling(coupling_info);

        # coupling with neighbor below
        if i < n_agents_y - 1 :
            coupling_info.neighbor_id_ = idx + n_agents_x;
            interface.register_coupling(coupling_info);

        idx = idx + 1;

# run MPC
interface.run_MPC(0, Tsim)

# get solution
solution_MPC = interface.get_solution('all')

# calculate global cost
global_cost_MPC = [0] * len(solution_MPC[0].agentState_.t_);
for i in range(0, n_agents) :
    for j in range(0, len(solution_MPC[0].agentState_.t_)) :
        global_cost_MPC[j] = global_cost_MPC[j] + solution_MPC[i].cost_[j]

# deregister agents
for i in range(0, n_agents_x) :
    for j in range(0, n_agents_y) :
        agentInfo.id_ = i*n_agents_x + j;
        interface.deregister_agent(agentInfo);

for i in range(0, n_agents_x) :
    for j in range(0, n_agents_y) :
        agentInfo.id_ = i*n_agents_x + j;

        interface.register_agent(agentInfo, x_init[agentInfo.id_], [0, 0]);

        x_des = [i, 0, j, 0];
        interface.set_desiredAgentState(agentInfo.id_, x_des, [0, 0]);

idx = 0
for i in range(0, n_agents_y) :
    for j in range (0, n_agents_x) :

        coupling_info.agent_id_ = idx;

        # coupling with neighbor on the left
        if j > 0 :
            coupling_info.neighbor_id_ = idx - 1;
            interface.register_coupling(coupling_info);

        # coupling with neighbor on the right
        if j < n_agents_x - 1 :
            coupling_info.neighbor_id_ = idx + 1;
            interface.register_coupling(coupling_info);

        # coupling with neighbor above
        if i > 0 :
            coupling_info.neighbor_id_ = idx - n_agents_x;
            interface.register_coupling(coupling_info);

        # coupling with neighbor below
        if i < n_agents_y - 1 :
            coupling_info.neighbor_id_ = idx + n_agents_x;
            interface.register_coupling(coupling_info);

        idx = idx + 1;

#run DMCP
interface.run_DMPC(0, Tsim)

# get solutions
solution_DMPC = interface.get_solution('all')

#calculate global cost
global_cost_DMPC = [0] * len(solution_DMPC[0].agentState_.t_);
for i in range(0, n_agents) :
    for j in range(0, len(solution_DMPC[0].agentState_.t_)) :
        global_cost_DMPC[j] = global_cost_DMPC[j] + solution_DMPC[i].cost_[j]
        

numpy.savetxt("t.csv", solution_MPC[0].agentState_.t_, delimiter=",")
numpy.savetxt("MPC.csv", global_cost_MPC, delimiter=",")
numpy.savetxt("DMPC.csv", global_cost_DMPC, delimiter=",")

fig = plt.figure()
fig.suptitle('Simulation example \n High-scaled system', y=1)

plt.plot(solution_MPC[0].agentState_.t_, global_cost_MPC, solution_MPC[0].agentState_.t_, global_cost_DMPC, 'r--')
plt.ylabel('Global cost')
plt.xlabel('Simulation time in s')

plt.show()