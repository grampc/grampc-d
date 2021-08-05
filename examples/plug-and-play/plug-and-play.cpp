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

	// initialize communication interface
	interface->initialize_central_communicationInterface();

	// set optimization info
	auto optimization_info = interface->optimizationInfo();
	optimization_info.COMMON_Nhor_ = 30;
	optimization_info.COMMON_Thor_ = 12;
	optimization_info.COMMON_dt_ = 0.1;
	optimization_info.GRAMPC_MaxGradIter_ = 40;
	optimization_info.GRAMPC_MaxMultIter_ = 1;
	optimization_info.ADMM_innerIterations_ = 2;
	optimization_info.ADMM_maxIterations_ = 20;
	optimization_info.ADMM_ConvergenceTolerance_ = 0.0005;
	interface->set_optimizationInfo(optimization_info);

	typeRNum Tsim = 20;

	// parameters for model
	typeRNum Omega = 1;
	typeRNum I = 1;
	typeRNum P_max = 0.1;
	typeRNum kappa = 1e-3;

	// parameters for cost function
	typeRNum P = 0.1;
	typeRNum Q = 1;
	typeRNum R = 0.01;

	// initial and desired states and controls
	std::vector<typeRNum> uinit = { 0 };
	std::vector<typeRNum> xinit = { 0, 0 };
	std::vector<typeRNum> xdes = { 0, 0 };
	std::vector<typeRNum> udes = {0};

	// register agents
	unsigned int n_agents = 3;
	auto agent = interface->agentInfo();
	agent.model_name_ = "smartGrid_agentModel";

	unsigned int agent_id = 0;
	typeRNum P0 = 0;
	typeRNum p = 1;
	agent.id_ = agent_id;
	agent.model_parameters_ = { I, Omega, kappa, P0, p };
	agent.cost_parameters_ = { 0, P, 0, Q, R };
	interface->register_agent(agent, xinit, uinit);
	interface->set_desiredAgentState(agent_id, xdes, udes);

	agent_id = 1;
	P0 = -0.01;
	p = 0;
	agent.id_ = agent_id;
	agent.model_parameters_ = { I, Omega, kappa, P0, p };
	agent.cost_parameters_ = { 0, P, 0, Q, R };
	interface->register_agent(agent, xinit, uinit);
	interface->set_desiredAgentState(agent_id, xdes, udes);

	// register couplings
	auto coupling_info = interface->couplingInfo();

	coupling_info.model_name_ = "smartGrid_couplingModel";
	coupling_info.model_parameters_ = { I, Omega, P_max };

	coupling_info.agent_id_ = 0;
	coupling_info.neighbor_id_ = 1;
	interface->register_coupling(coupling_info);

	coupling_info.agent_id_ = 1;
	coupling_info.neighbor_id_ = 0;
	interface->register_coupling(coupling_info);

	// run distributed controller
	interface->run_DMPC(0, Tsim);

	// register agent
	agent_id = 2;
	agent.id_ = agent_id;
	interface->register_agent(agent, xinit, uinit);
	interface->set_desiredAgentState(agent_id, xdes, udes);

	// register couplings
	coupling_info.agent_id_ = 1;
	coupling_info.neighbor_id_ = 2;
	interface->register_coupling(coupling_info);

	coupling_info.agent_id_ = 2;
	coupling_info.neighbor_id_ = 1;
	interface->register_coupling(coupling_info);

	// continue simulation
	interface->run_DMPC(Tsim, Tsim);

	interface->print_solution_to_file("all");

	return 0;
}
