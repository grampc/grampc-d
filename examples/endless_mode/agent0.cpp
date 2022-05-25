/* This file is part of GRAMPC-D - (https://github.com/grampc-d/grampc-d.git)
 *
 * GRAMPC-D -- A software framework for distributed model predictive control (DMPC)
 * based on the alternating direction method of multipliers (ADMM).
 *
 * Copyright 2020 by Daniel Burk, Andreas Voelz, Knut Graichen
 * All rights reserved.
 *
 * GRAMPC-D is distributed under the BSD-3-Clause license, see LICENSE.txt
 *
 */

#include "grampcd/interface/dmpc_interface.hpp"

int main(int argc, char** argv)
{
	// create interface
	const auto interface = std::make_shared<grampcd::DmpcInterface>();

	// show logging
	interface->set_print_message(true);
	interface->set_print_warning(true);
	interface->set_print_error(true);

	const int agent_id = 0;

	// set communication info
	auto comm_info_coordinator = interface->communication_info();
    comm_info_coordinator.ip_ = "127.0.0.1";
	comm_info_coordinator.port_ = "7777";

	// initialize communication interface
	interface->initialize_local_communicationInterface_as_agent(comm_info_coordinator);

	// parameters for cost function
	typeRNum P = 1; typeRNum Q = 1; typeRNum R = 0.1;

	// parameters for model
	typeRNum A = 0.1; typeRNum a = 0.005; typeRNum d = 0.01;

	// initial and desired states and controls
	std::vector<typeRNum> xinit(1, 1);
	std::vector<typeRNum> uinit(1, 0.0);
	std::vector<typeRNum> xdes(1, 2.0);
	std::vector<typeRNum> udes(1, 0.0);

	// register agent
	auto agent = interface->agent_info();
	agent.id_ = agent_id;
	agent.model_name_ = "water_tank_agentModel";
	agent.model_parameters_ = { A, 1, 0 };
	agent.cost_parameters_ = { 0, 0, R };

	interface->register_agent(agent, xinit, uinit, xdes, udes);

	// register coupling
	auto coupling_info = interface->coupling_info();
	coupling_info.agent_id_ = agent_id;
	coupling_info.model_name_ = "water_tank_couplingModel";
	coupling_info.model_parameters_ = { A, a };

	coupling_info.neighbor_id_ = 1;
    interface->register_coupling(coupling_info);

	interface->cap_stored_data(100);

	while (true)
		interface->wait_blocking_s(1);

	return 0;
}
