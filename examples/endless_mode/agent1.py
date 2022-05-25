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

# show logging
interface.set_print_message(1)
interface.set_print_warning(1)
interface.set_print_error(1)

agent_id = 1

# set communication info
comm_info_coordinator = grampcd_interface.communication_info()
comm_info_coordinator.ip_ = '127.0.0.1'
comm_info_coordinator.port_ = '7777'

# initialize communication interface
interface.initialize_local_communicationInterface_as_agent(comm_info_coordinator);

# parameters for cost function
P = 1
Q = 1
R = 0.1;

# parameters for model
A = 0.1
a = 0.005
d = 0.01

# initial and desired states and controls
xinit = [1]
uinit = [0]
xdes = [2] 
udes = [0];

# register agent
agent = grampcd_interface.agent_info()
agent.id_ = agent_id;
agent.model_name_ = "water_tank_agentModel";
agent.model_parameters_ = [A, 0, d];
agent.cost_parameters_ = [P, Q, 0];

interface.register_agent(agent, xinit, uinit, xdes, udes);

# register coupling
coupling_info = grampcd_interface.coupling_info();
coupling_info.agent_id_ = agent_id;
coupling_info.model_name_ = 'water_tank_couplingModel';
coupling_info.model_parameters_ = [A, a];

coupling_info.neighbor_id_ = 0;
interface.register_coupling(coupling_info);

interface.cap_stored_data(100)

fig = plt.figure()

max_cost = 0.01
while 1:
	solution = interface.get_solution(1)
	
	fig.suptitle('Simulation example \n Coupled water tanks \n Agent 1', y=1)
	
	plt.subplot(1, 3, 1)
	plt.plot(solution.agentState_.t_, solution.agentState_.x_)
	plt.ylabel('Waterheight')
	axes = plt.gca()
	axes.set_ylim([0, 3.5])
	
	plt.subplot(1, 3, 2)
	plt.plot(solution.agentState_.t_, solution.agentState_.u_)
	plt.ylabel('Control trajectory')
	axes = plt.gca()
	axes.set_ylim([0, 0.25])

	plt.subplot(1, 3, 3)
	plt.plot(solution.agentState_.t_, solution.cost_)
	plt.ylabel('Cost')
	axes = plt.gca()
	if len(solution.cost_) > 0:
		max_cost = max(max_cost, max(solution.cost_))

	axes.set_ylim([0, max_cost])
	
	plt.subplots_adjust(wspace = 0.8)
	plt.draw()
	plt.pause(0.5)
	plt.clf()