/* This file is part of GRAMPC-D - (https://github.com/grampc-d/grampc-d.git)
 *
 * GRAMPC-D -- A software framework for distributed model predictive control (DMPC)
 * 
 *
 * Copyright 2023 by Daniel Burk, Maximilian Pierer von Esch, Andreas Voelz, Knut Graichen
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
	auto optimization_info = interface->optimization_info();


	optimization_info.COMMON_Nhor_ = 16;
	optimization_info.COMMON_Thor_ = 4;
	optimization_info.COMMON_dt_ = 0.5;
	optimization_info.GRAMPC_MaxGradIter_ = 40;
	optimization_info.GRAMPC_MaxMultIter_ = 2;
	optimization_info.COMMON_Solver_ = "ADMM";
	optimization_info.ADMM_innerIterations_ = 5;
	optimization_info.ADMM_maxIterations_ = 100;
	optimization_info.ADMM_ConvergenceTolerance_ = 0.009;
	optimization_info.COMMON_DebugCost_ = true;

	optimization_info.ASYNC_Active_ = false;
	optimization_info.ASYNC_Delay_ = 0;

	bool approx = true;
	optimization_info.APPROX_ApproximateCost_ = approx;
	optimization_info.APPROX_ApproximateConstraints_ = approx;
	optimization_info.APPROX_ApproximateDynamics_ = approx;

	interface->set_optimizationInfo(optimization_info);

	const typeRNum Tsim = 0.05;

	// parameters for cost function
	typeRNum P = 1;
	typeRNum Q = 1;
	typeRNum R = 0.1;

	// parameters for model
	typeRNum A = 0.1;
	typeRNum a = 0.005;
	typeRNum d = 0.01;

	// initial and desired states and controls
	std::vector<typeRNum> xinit(1, 0.5);
	std::vector<typeRNum> uinit(1, 0.0);
	std::vector<typeRNum> xdes(1, 2.0);
	std::vector<typeRNum> udes(1, 0.0);

	// register agent
	auto agent = interface->agent_info();
	agent.model_name_ = "water_tank_agentModel";

	int agent_id = 0;
	agent.id_ = agent_id;
	agent.model_parameters_ = { A, 1, 0 };
	agent.cost_parameters_ = { 0, 0, R };

	interface->register_agent(agent, xinit, uinit, xdes, udes);

	const int n_agents = 5;

	for (unsigned int i = 1; i < n_agents - 1; ++i)
	{
		agent.id_ = i;
		agent.model_parameters_ = { A, 0, 0 };
		agent.cost_parameters_ = { 0, 0, 0 };
		interface->register_agent(agent, xinit, uinit, xdes, udes);
	}

	agent_id = n_agents - 1;
	agent.id_ = agent_id;
	agent.model_parameters_ = { A, 0, d };
	agent.cost_parameters_ = { P, Q, 0 };
	interface->register_agent(agent, xinit, uinit, xdes, udes);

	// register couplings
	auto coupling_info = interface->coupling_info();
	coupling_info.model_name_ = "water_tank_couplingModel";
	coupling_info.model_parameters_ = { A, a };

	for (unsigned int i = 0; i < n_agents; ++i)
	{
		coupling_info.agent_id_ = i;
		if (i > 0)
		{
			coupling_info.neighbor_id_ = i - 1;
			interface->register_coupling(coupling_info);
		}
		if (i < n_agents - 1)
		{
			coupling_info.neighbor_id_ = i + 1;
			interface->register_coupling(coupling_info);
		}
	}

	// run distributed controller
	interface->run_DMPC(0, Tsim);

	interface->print_solution_to_file("all", "withApprox_");

	// deregister agents
	for (unsigned int i = 0; i < n_agents; ++i)
	{
		agent.id_ = i;
		interface->deregister_agent(agent);
	}

	approx = false;
	optimization_info.APPROX_ApproximateCost_ = approx;
	optimization_info.APPROX_ApproximateConstraints_ = approx;
	optimization_info.APPROX_ApproximateDynamics_ = approx;

	interface->set_optimizationInfo(optimization_info);

	// register agents
	agent_id = 0;
	agent.id_ = agent_id;
	agent.model_parameters_ = { A, 1, 0 };
	agent.cost_parameters_ = { 0, 0, R };

	interface->register_agent(agent, xinit, uinit, xdes, udes);

	for (unsigned int i = 1; i < n_agents - 1; ++i)
	{
		agent.id_ = i;
		agent.model_parameters_ = { A, 0, 0 };
		agent.cost_parameters_ = { 0, 0, 0 };
		interface->register_agent(agent, xinit, uinit, xdes, udes);
	}

	agent_id = n_agents - 1;
	agent.id_ = agent_id;
	agent.model_parameters_ = { A, 0, d };
	agent.cost_parameters_ = { P, Q, 0 };
	interface->register_agent(agent, xinit, uinit, xdes, udes);

	// register couplings
	for (unsigned int i = 0; i < n_agents; ++i)
	{
		coupling_info.agent_id_ = i;
		if (i > 0)
		{
			coupling_info.neighbor_id_ = i - 1;
			interface->register_coupling(coupling_info);
		}
		if (i < n_agents - 1)
		{
			coupling_info.neighbor_id_ = i + 1;
			interface->register_coupling(coupling_info);
		}
	}

	// run distributed controller
	interface->run_DMPC(0, Tsim);

	interface->print_solution_to_file("all", "withoutApprox_");

	return 0;
}
