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

	// initialize communication interface
	interface->initialize_central_communicationInterface();

	// set optimization info
	auto optimization_info = interface->optimizationInfo();
	optimization_info.COMMON_Nhor_ = 21;
	optimization_info.COMMON_Thor_ = 1;
	optimization_info.COMMON_dt_ = 0.1;
	optimization_info.GRAMPC_MaxGradIter_ = 10;
	optimization_info.GRAMPC_MaxMultIter_ = 1;
	optimization_info.ADMM_maxIterations_ = 10;
	optimization_info.ADMM_ConvergenceTolerance_ = 0;

	interface->set_optimizationInfo(optimization_info);

	const typeRNum Tsim = 3;

	// parameters for cost function
	ctypeRNum P = 1; ctypeRNum Q = 1; ctypeRNum R = 0.1;
	const std::vector<typeRNum> model_parameters = { 1, 1, 1 };
	const std::vector<typeRNum> cost_parameters = { 0, P, 0, Q, R };

	// initial and desired states and controls
	std::vector<typeRNum> uinit = { 0.0 };
	std::vector<typeRNum> xdes = { 0.0, 0.0 };
	std::vector<typeRNum> udes = { 0.0 };

	// register agent
	unsigned int agent_id = 0;
	std::vector<typeRNum> xinit = { 0.5, 0.0 };
	auto agent = interface->agentInfo();
	agent.id_ = agent_id;
	agent.model_name_ = "vdp_agentModel";
	agent.model_parameters_ = model_parameters;
	agent.cost_parameters_ = cost_parameters;
	interface->register_agent(agent, xinit, uinit);
	interface->set_desiredAgentState(agent_id, xdes, udes);

	agent_id = 1;
	xinit = { 0.25, 0.0 };
	agent = interface->agentInfo();
	agent.id_ = agent_id;
	agent.model_name_ = "vdp_agentModel";
	agent.model_parameters_ = model_parameters;
	agent.cost_parameters_ = cost_parameters;
	interface->register_agent(agent, xinit, uinit);
	interface->set_desiredAgentState(agent_id, xdes, udes);

	agent_id = 2;
	xinit = { 0.75, 0.0 };
	agent = interface->agentInfo();
	agent.id_ = agent_id;
	agent.model_name_ = "vdp_agentModel";
	agent.model_parameters_ = model_parameters;
	agent.cost_parameters_ = cost_parameters;
	interface->register_agent(agent, xinit, uinit);
	interface->set_desiredAgentState(agent_id, xdes, udes);

	// register coupling
	auto coupling_info = interface->couplingInfo();
	coupling_info.model_name_ = "vdp_synchronize_couplingModel";
	coupling_info.model_parameters_ = { 1 };
	coupling_info.cost_parameters_ = { 1 };

	coupling_info.agent_id_ = 0;
	coupling_info.neighbor_id_ = 1;
	interface->register_coupling(coupling_info);

	coupling_info.agent_id_ = 1;
	coupling_info.neighbor_id_ = 0;
	interface->register_coupling(coupling_info);

	coupling_info.agent_id_ = 1;
	coupling_info.neighbor_id_ = 2;
	interface->register_coupling(coupling_info);

	coupling_info.agent_id_ = 2;
	coupling_info.neighbor_id_ = 1;
	interface->register_coupling(coupling_info);

	// run distributed controller
	interface->run_DMPC(0, Tsim);

	// print solution to file
	interface->print_solution_to_file("all");

	return 0;
}
