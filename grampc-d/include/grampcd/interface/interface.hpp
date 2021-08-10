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

#pragma once

#include "grampcd/util/class_forwarding.hpp"

namespace grampcd
{
	class Interface
	{
	public:

		virtual ~Interface();

		/*Initialize the central communication interface.*/
		virtual void initialize_central_communicationInterface(const int number_of_threads = 0) = 0;
		/*Initialize the local communication interface.*/
		virtual void initialize_local_communicationInterface_as_agent(const CommunicationInfo& adress_coordinator) = 0;
		/*Initialize the local communication interface.*/
		virtual void initialize_local_communicationInterface_as_coordinator(const unsigned short port) = 0;

		/*Register an agent.*/
		virtual void register_agent(const AgentInfo& info, const std::vector<typeRNum>& x_init, const std::vector<typeRNum>& u_init) = 0;
		/*Set the desired agent state of an agent.*/
		virtual void set_desiredAgentState(const int agent_id, const std::vector<typeRNum>& x_des, const std::vector<typeRNum>& u_des) = 0;
		/*De-register an agent.*/
		virtual void deregister_agent(const AgentInfo& info) = 0;
		/*Set intial state*/
		virtual void set_initialState(const unsigned int agent_id, const std::vector<typeRNum>& x_init) = 0;

		/*Register a coupling between agents.*/
		virtual void register_coupling(const CouplingInfo& info) = 0;
		/*De-register a coupling between agents.*/
		virtual void deregister_coupling(const CouplingInfo& info) = 0;

		/*Run a centralized controller.*/
		virtual void run_MPC(ctypeRNum Tsim, ctypeRNum t0) = 0;
		/*Run a centralized controller in endless model*/
		virtual void run_MPC() = 0;
		/* Run a distributed controller.*/
		virtual void run_DMPC(ctypeRNum Tsim, ctypeRNum t0) = 0;
		/*Run a distributed controller in endless mode.*/
		virtual void run_DMPC() = 0;

		/*Returns the solution of an agent.*/
		virtual SolutionPtr get_solution(const unsigned int agent_id) const = 0;
		/*Returns the solutions of a set of agents.*/
		virtual std::vector< SolutionPtr > get_solution(const std::string& agents) const = 0;
		/*Resets the solution of an agent.*/
		virtual void reset_solution(const unsigned int agent_id) = 0;
		/*Resets the solutions of a set of agents.*/
		virtual void reset_solution(const std::string& agents) = 0;
		/*Print solution of an agent to a file using a specific prefix for the file name.*/
		virtual void print_solution_to_file(const unsigned int agent_id, const std::string& prefix = "Solution_agent") const = 0;
		/*Print solutions of a set of agents to a file using a specific prefix for the file name.*/
		virtual void print_solution_to_file(const std::string& agents, const std::string& prefix = "Solution_agent") const = 0;

		/*Returns the current optimization info.*/
		virtual OptimizationInfoPtr get_optimizationInfo() const = 0;
		/*Sets an optimization info.*/
		virtual void set_optimizationInfo(const OptimizationInfo& optimization_info) = 0;

		/*Wait for some seconds.*/
		virtual void wait_blocking_s(const unsigned int s) = 0;
		/*Wait for active connections.*/
		virtual void wait_for_connections(const int agents, const int couplings) = 0;
		/*Set as passive.*/
		virtual void set_passive() = 0;
		/*Send flags to a set of agents.*/
		virtual void send_flag_to_agents(const std::string& agents) const = 0;
		/*Send flags to a set of agents.*/
		virtual void send_flag_to_agents(const std::vector<int>& agent_ids) const = 0;
		/*Send flags to an agent.*/
		virtual void send_flag_to_agents(const int agent_id) const = 0;
		/*Wait for flag from coordinator*/
		virtual void waitFor_flag_from_coordinator() = 0;

		/*Activate to simulation real time.*/
		virtual void simulate_realtime(const bool realtime) = 0;
		/*Cap the stored data.*/
		virtual void cap_stored_data(const unsigned int data_points) = 0;

		/*Print messages of type base*/
		virtual void set_print_base(const bool print) = 0;
		/*Print messages of type message*/
		virtual void set_print_message(const bool print) = 0;
		/*Print messages of type warning*/
		virtual void set_print_warning(const bool print) = 0;
		/*Print messages of type error*/
		virtual void set_print_error(const bool print) = 0;
		/*Print messages of type progress bar*/
		virtual void set_print_progressbar(const bool print) = 0;
	};
}