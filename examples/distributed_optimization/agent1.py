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
P = 1; Q = 1; R = 0.1;

# inital and desired states and controls
xinit = [0.5, 0];
uinit = [0];
xdes = [0, 0]; 
udes = [0];

# register agent
agent = grampcd_interface.agent_info()
agent.id_ = agent_id;
agent.model_name_ = "vdp_agentModel";
agent.model_parameters_ = [1, 1, 1]
agent.cost_parameters_ = [P, P, Q, Q, R]

interface.register_agent(agent, xinit, uinit, xdes, udes);

# register couplings
coupling_info = grampcd_interface.coupling_info();
coupling_info.agent_id_ = agent_id;
coupling_info.model_name_ = 'vdp_linear_couplingModel';
coupling_info.model_parameters_ = [1];

coupling_info.neighbor_id_ = 0;
interface.register_coupling(coupling_info);

coupling_info.neighbor_id_ = 2;
interface.register_coupling(coupling_info);

# wait for flag
interface.waitFor_flag_from_coordinator()