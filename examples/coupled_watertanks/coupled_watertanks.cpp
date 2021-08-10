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
	optimization_info.COMMON_Nhor_ = 21;
	optimization_info.COMMON_Thor_ = 5;
	optimization_info.COMMON_dt_ = 0.1;
	optimization_info.GRAMPC_MaxGradIter_ = 10;
	optimization_info.GRAMPC_MaxMultIter_ = 2;
	optimization_info.ADMM_maxIterations_ = 10;
	optimization_info.ADMM_ConvergenceTolerance_ = 0.02;

	bool approx = true;
	optimization_info.APPROX_ApproximateCost_ = approx;
	optimization_info.APPROX_ApproximateConstraints_ = approx;
	optimization_info.APPROX_ApproximateDynamics_ = approx;

	interface->set_optimizationInfo(optimization_info);

	const typeRNum Tsim = 25;

	// parameters for cost function
	typeRNum P = 1;
	typeRNum Q = 1;
	typeRNum R = 0.1;

	// parameters for model
	typeRNum A = 0.1;
	typeRNum a = 0.005;
	typeRNum d = 0.01;

	// inital and desired states and controls
	std::vector<typeRNum> xinit(1, 0.5);
	std::vector<typeRNum> uinit(1, 0.0);
	std::vector<typeRNum> xdes(1, 2.0);
	std::vector<typeRNum> udes(1, 0.0);

	// register agents
	auto agent = interface->agentInfo();
	agent.model_name_ = "water_tank_agentModel";

	int agent_id = 1;
	agent.id_ = agent_id;
	agent.model_parameters_ = { A, 1, 0 };
	agent.cost_parameters_ = { 0, 0, R };
	interface->register_agent(agent, xinit, uinit);
	interface->set_desiredAgentState(agent_id, xdes, udes);

	agent_id = 2;
	agent.id_ = agent_id;
	agent.model_parameters_ = { A, 0, 0 };
	agent.cost_parameters_ = { 0, 0, 0 };
	interface->register_agent(agent, xinit, uinit);
	interface->set_desiredAgentState(agent_id, xdes, udes);

	agent_id = 3;
	agent.id_ = agent_id;
	agent.model_parameters_ = { A, 0, 0 };
	agent.cost_parameters_ = { 0, 0, 0 };
	interface->register_agent(agent, xinit, uinit);
	interface->set_desiredAgentState(agent_id, xdes, udes);

	agent_id = 4;
	agent.id_ = agent_id;
	agent.model_parameters_ = { A, 0, d };
	agent.cost_parameters_ = { P, Q, 0 };
	interface->register_agent(agent, xinit, uinit);
	interface->set_desiredAgentState(agent_id, xdes, udes);

	// register couplings
	auto coupling_info = interface->couplingInfo();
	coupling_info.model_name_ = "water_tank_couplingModel";
	coupling_info.model_parameters_ = { A, a };

	coupling_info.agent_id_ = 1;
	coupling_info.neighbor_id_ = 2;
	interface->register_coupling(coupling_info);

	coupling_info.agent_id_ = 1;
	coupling_info.neighbor_id_ = 3;
	interface->register_coupling(coupling_info);

	coupling_info.agent_id_ = 2;
	coupling_info.neighbor_id_ = 1;
	interface->register_coupling(coupling_info);

	coupling_info.agent_id_ = 2;
	coupling_info.neighbor_id_ = 4;
	interface->register_coupling(coupling_info);

	coupling_info.agent_id_ = 3;
	coupling_info.neighbor_id_ = 1;
	interface->register_coupling(coupling_info);

	coupling_info.agent_id_ = 3;
	coupling_info.neighbor_id_ = 4;
	interface->register_coupling(coupling_info);

	coupling_info.agent_id_ = 4;
	coupling_info.neighbor_id_ = 2;
	interface->register_coupling(coupling_info);

	coupling_info.agent_id_ = 4;
	coupling_info.neighbor_id_ = 3;
	interface->register_coupling(coupling_info);

	// activate progress bar
	interface->set_print_progressbar(true);

	// run distributed controller
	interface->run_DMPC(0, Tsim);

	// print solution to file
	interface->print_solution_to_file("all");

	return 0;
}
