# This file is part of GRAMPC-D - (https://github.com/grampc-d/grampc-d.git)
#
# GRAMPC-D -- A software framework for distributed model predictive control (DMPC)
# based on the alternating direction method of multipliers (ADMM).
#
# Copyright 2020 by Daniel Burk, Andreas Voelz, Knut Graichen
# All rights reserved.
#
# GRAMPC-D is distributed under the BSD-3-Clause license, see LICENSE.txt

add_executable(endless_mode_coordinator coordinator.cpp)
target_link_libraries(endless_mode_coordinator dmpc)

add_executable(endless_mode_agent0 agent0.cpp)
target_link_libraries(endless_mode_agent0 dmpc)

add_executable(endless_mode_agent1 agent1.cpp)
target_link_libraries(endless_mode_agent1 dmpc)