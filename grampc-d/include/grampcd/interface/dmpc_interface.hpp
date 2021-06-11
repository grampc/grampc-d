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

#include "grampcd/interface/interface.hpp"

#include "grampcd/info/communication_info.hpp"
#include "grampcd/info/optimization_info.hpp"
#include "grampcd/info/agent_info.hpp"
#include "grampcd/info/coupling_info.hpp"

#include "grampcd/util/class_forwarding.hpp"

namespace grampcd
{
	/*@brief Interface for the DMPC-framework*/
	class DmpcInterface : public Interface
	{
	public:
		DmpcInterface();

		/*Create object of an agent info*/
		AgentInfo agentInfo() const;
		/*Create object of a coupling info*/
		CouplingInfo couplingInfo() const;
		/*Create object of an optimization info*/
		OptimizationInfo optimizationInfo() const;
		/*Create object of a communication info*/
		CommunicationInfo communicationInfo() const;

		/*Initialize the central communication interface.*/
		void initialize_central_communicationInterface(int number_of_threads = 0) override;
		/*Initialize the local communication interface.*/
		void initialize_local_communicationInterface_as_agent(CommunicationInfo adress_coordinator) override;
		/*Initialize the local communication interface.*/
		void initialize_local_communicationInterface_as_coordinator(unsigned short port) override;

		/*Register an agent.*/
		void register_agent(AgentInfo info, std::vector<typeRNum> x_init, std::vector<typeRNum> u_init) override;
		/*Set the desired agent state of an agent.*/
		void set_desiredAgentState(int agent_id, std::vector<typeRNum> x_des, std::vector<typeRNum> u_des) override;
		/*De-register an agent.*/
		void deregister_agent(AgentInfo info) override;
		/*Set initial state*/
		void set_initialState(const unsigned int agent_id, const std::vector<typeRNum>& x_init) override;

		/*Register a coupling between agents.*/
		void register_coupling(CouplingInfo info) override;
		/*De-register a coupling between agents.*/
		void deregister_coupling(CouplingInfo info) override;

		/*Run a centralized controller.*/
		void run_MPC(typeRNum Tsim, typeRNum t0) override;
		/*Run a centralized controller in endless model*/
		void run_MPC() override;
		/* Run a distributed controller.*/
		void run_DMPC(typeRNum Tsim, typeRNum t0) override;
		/*Run a distributed controller in endless mode.*/
		void run_DMPC() override;

		/*Returns the solution of an agent.*/
		SolutionPtr get_solution(unsigned int agent_id) const override;
		/*Returns the solutions of a set of agents.*/
		std::vector< SolutionPtr > get_solution(std::string agents) const override;
		/*Resets the solution of an agent.*/
		virtual void reset_solution(unsigned int agent_id) override;
		/*Resets the solutions of a set of agents.*/
		virtual void reset_solution(std::string agents) override;
		/*Print solution of an agent to a file.*/
		virtual void print_solution_to_file(const unsigned int agent_id, const std::string prefix = "Solution_agent") const override;
		/*Print solutions of a set of agents to a file.*/
		virtual void print_solution_to_file(const std::string agents, const std::string prefix = "Solution_agent") const override;

		/*Returns the current optimization info.*/
		OptimizationInfoPtr get_optimizationInfo() const override;
		/*Sets an optimization info.*/
		void set_optimizationInfo(OptimizationInfo optimization_info) override;

		/*Wait for some seconds.*/
		void wait_blocking_s(unsigned int s) override;
		/*Wait for active connections.*/
		void wait_for_connections(int agents, int couplings) override;
		/*Set as passive.*/
		void set_passive() override;
		/*Send flags to a set of agents.*/
		void send_flag_to_agents(const std::string agents) const override;
		/*Send flags to a set of agents.*/
		void send_flag_to_agents(const std::vector<int> agent_ids) const override;
		/*Send flags to an agent.*/
		void send_flag_to_agents(const int agent_id) const override;
		/*Wait for flag from coordinator*/
		void waitFor_flag_from_coordinator() override;

		/*Activate to simulation real time.*/
		void simulate_realtime(bool realtime) override;
		/*Cap the stored data.*/
		void cap_stored_data(unsigned int data_points) override;

		/*Print messages of type base*/
		void set_print_base(bool print) override;
		/*Print messages of type message*/
		void set_print_message(bool print) override;
		/*Print messages of type warning*/
		void set_print_warning(bool print) override;
		/*Print messages of type error*/
		void set_print_error(bool print) override;

	private:
		/*Run a centralized controller.*/
		void run_MPC(const std::vector<AgentPtr>& agents, SimulatorPtr simulator, const OptimizationInfo& oi, typeRNum Tsim, typeRNum t_0);
		/*Run a centralized controller in endless model*/
		void run_MPC(const std::vector<AgentPtr>& agents, SimulatorPtr simulator, const OptimizationInfo& oi);
		/*Run a distributed controller.*/
		void run_DMPC(SimulatorPtr simulator, const OptimizationInfo& oi, typeRNum Tsim, typeRNum t_0);
		/*Run a centralized distributed in endless model*/
		void run_DMPC(SimulatorPtr simulator, const OptimizationInfo& oi);

		LoggingPtr log_;

		OptimizationInfo optimizationInfo_;
		CommunicationInfo communication_info_;
		CommunicationInterfacePtr communication_interface_;
		CoordinatorPtr coordinator_;
		SimulatorPtr simulator_;
		ModelFactoryPtr factory_;
		std::vector< SolutionPtr> solutions_;

		std::vector< AgentPtr> agents_;

		bool realtime_ = false;
	};
}