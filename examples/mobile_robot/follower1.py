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

# show logging
interface.set_print_message(1)
interface.set_print_warning(1)
interface.set_print_error(1)

# initialize communication interface
comm_info_coordinator = grampcd_interface.communication_info()
comm_info_coordinator.ip_ = '127.0.0.1'
comm_info_coordinator.port_ = '7777'
interface.initialize_local_communicationInterface_as_agent(comm_info_coordinator);

# register agent
agent = grampcd_interface.agent_info()
agent.model_name_ = "mobile_robot_follower_model"
agent.id_ = 2
interface.register_agent(agent, [0.25, 0, 0], [0, 0], [0, 0, 0], [0, 0])

# register coupling
coupling_info = grampcd_interface.coupling_info()
coupling_info.model_name_ = 'mobile_robot_coupling_model'
coupling_info.cost_parameters_ = [1, 5, 0.1, 0.125, 0.0125]
coupling_info.agent_id_ = 2
coupling_info.neighbor_id_ = 1
coupling_info.model_parameters_ = [-0.2, -0.2]
interface.register_coupling(coupling_info)

# wait for flag
interface.waitFor_flag_from_coordinator()