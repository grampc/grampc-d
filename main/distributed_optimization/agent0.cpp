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

#include "dmpc/interface/dmpc_interface.hpp"

int main(int argc, char** argv)
{
	// create interface
	const auto interface = std::make_shared<dmpc::DmpcInterface>();

	// show logging
	interface->set_print_message(true);
	interface->set_print_warning(true);
	interface->set_print_error(true);

	const int agent_id = 0;

	// set communication info
	auto comm_info_coordinator = interface->communicationInfo();
    comm_info_coordinator.ip_ = "127.0.0.1";
	comm_info_coordinator.port_ = "7777";

	// initialize communication interface
	interface->initialize_local_communicationInterface_as_agent(comm_info_coordinator);

	// parameters for cost function
	typeRNum P = 1; typeRNum Q = 1; typeRNum R = 0.1;

	// initial and desired states and controls
	std::vector<typeRNum> xinit = {0.5, 0.0};
	std::vector<typeRNum> uinit = {0.0};
	std::vector<typeRNum> xdes = {0.0, 0.0};
	std::vector<typeRNum> udes = {0.0};

	// register agent
	auto agent = interface->agentInfo();
	agent.id_ = agent_id;
	agent.model_name_ = "vdp_agentModel";
	agent.model_parameters_ = { 1, 1, 1 };
	agent.cost_parameters_ = { P, P, Q, Q, R };

	interface->register_agent(agent, xinit, uinit);
	interface->set_desiredAgentState(agent_id, xdes, udes);

	// register coupling
	auto coupling_info = interface->couplingInfo();
	coupling_info.agent_id_ = agent_id;
	coupling_info.model_name_ = "vdp_linear_couplingModel";
	coupling_info.model_parameters_ = { 1 };

	coupling_info.neighbor_id_ = 1;
    interface->register_coupling(coupling_info);

	// wait for flag
	interface->waitFor_flag_from_coordinator();

	return 0;
}
