# This file is part of GRAMPC-D - (https://github.com/grampc-d/grampc-d.git)
#
# GRAMPC-D -- A software framework for distributed model predictive control (DMPC)
# 
#
# Copyright 2023 by Daniel Burk, Maximilian Pierer von Esch, Andreas Voelz, Knut Graichen
# All rights reserved.
#
# GRAMPC-D is distributed under the BSD-3-Clause license, see LICENSE.txt

add_executable(distributed_optimization_coordinator coordinator.cpp)
target_link_libraries(distributed_optimization_coordinator grampcd)

add_executable(distributed_optimization_agent0 agent0.cpp)
target_link_libraries(distributed_optimization_agent0 grampcd)

add_executable(distributed_optimization_agent1 agent1.cpp)
target_link_libraries(distributed_optimization_agent1 grampcd)

add_executable(distributed_optimization_agent2 agent2.cpp)
target_link_libraries(distributed_optimization_agent2 grampcd)