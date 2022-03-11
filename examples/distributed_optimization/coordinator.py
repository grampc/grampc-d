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

# show logging
interface.set_print_message(1)
interface.set_print_warning(1)
interface.set_print_error(1)

# initialize communication interface
interface.initialize_local_communicationInterface_as_coordinator(7777)

# set optimization info
optimization_info = grampcd_interface.OptimizationInfo()
optimization_info.COMMON_Nhor_ = 21
optimization_info.COMMON_Thor_ = 1
optimization_info.COMMON_dt_ = 0.1
optimization_info.GRAMPC_MaxGradIter_ = 10
optimization_info.GRAMPC_MaxMultIter_ = 1
optimization_info.ADMM_maxIterations_ = 5;
optimization_info.ADMM_ConvergenceTolerance_ = 0;

optimization_info.ADMM_DebugCost_ = 0
optimization_info.ASYNC_Delay_ = 0

interface.set_optimizationInfo(optimization_info)

Tsim = 5;

# wait for agents to connect
interface.wait_for_connections(3, 4)

# run distributed controller
interface.run_DMPC(0, Tsim)

# get solutions
solution_DMPC = interface.get_solution('all')

interface.send_flag_to_agents('all')

#calculate global cost
global_cost_DMPC = [0] * len(solution_DMPC[0].agentState_.t_);
for i in range(0, 3) :
    for j in range(0, len(solution_DMPC[0].agentState_.t_)) :
        global_cost_DMPC[j] = global_cost_DMPC[j] + solution_DMPC[i].cost_[j]

fig = plt.figure()
fig.suptitle('Simulation example \n Distributed optimization', y=1)

plt.plot(solution_DMPC[0].agentState_.t_, global_cost_DMPC, 'r--')
plt.ylabel('Global cost')
plt.xlabel('Simulation time in s')

plt.show()
