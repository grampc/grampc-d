# This file is part of GRAMPC-D - (https://github.com/DanielBurk/GRAMPC-D.git)
#
# GRAMPC-D -- A software framework for distributed model predictive control (DMPC)
# based on the alternating direction method of multipliers (ADMM).
#
# Copyright 2020 by Daniel Burk, Andreas Voelz, Knut Graichen
# All rights reserved.
#
# GRAMPC-D is distributed under the BSD-3-Clause license, see LICENSE.txt

add_executable(evaluate_neighborApproximation evaluate_neighborApproximation.cpp)
target_link_libraries(evaluate_neighborApproximation dmpc)