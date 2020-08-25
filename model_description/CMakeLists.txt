# This file is part of GRAMPC-D - (https://github.com/DanielBurk/GRAMPC-D.git)
#
# GRAMPC-D -- A software framework for distributed model predictive control (DMPC)
# based on the alternating direction method of multipliers (ADMM).
#
# Copyright 2020 by Daniel Burk, Andreas Voelz, Knut Graichen
# All rights reserved.
#
# GRAMPC-D is distributed under the BSD-3-Clause license, see LICENSE.txt

add_library(factory
	"general_model_factory.cpp"

	water_tank/src/water_tank_agent_model.cpp
	water_tank/src/water_tank_coupling_model.cpp

	van_der_pol_oscillator/src/vdp_agent_model.cpp
	van_der_pol_oscillator/src/vdp_linear_coupling_model.cpp
	van_der_pol_oscillator/src/vdp_nonlinear_coupling_model.cpp

	scalable_spring_mass_system/src/ssms_agent_model.cpp
	scalable_spring_mass_system/src/ssms_coupling_model.cpp

	scalable_spring_mass_system_2D/src/ssms2d_agent_model.cpp
	scalable_spring_mass_system_2D/src/ssms2d_coupling_model.cpp

	scalable_spring_mass_system_3D/src/ssms3d_agent_model.cpp
	scalable_spring_mass_system_3D/src/ssms3d_coupling_model.cpp

	oscillators/src/oscillators_agent_model.cpp
	oscillators/src/oscillators_coupling_model.cpp

	smart_grid/src/smartGrid_agentModel.cpp
	smart_grid/src/smartGrid_couplingModel.cpp
)

target_link_libraries(factory dmpc)

include_directories(water_tank)

include_directories(van_der_pol_oscillator)

include_directories(copter)