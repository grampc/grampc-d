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

#include "grampcd/comm/communication_interface.hpp"

#include "grampcd/info/communication_info.hpp"

#include "grampcd/util/class_forwarding.hpp"

#include "asio.hpp"

namespace grampcd
{

	/*!
	 * @brief Interface for communication between agents
	 *
	 * Implementation of a communication interface that has direct access to all agents.
	 * This means that all agents run within the same process.
	*/
	class CommunicationInterfaceCentral : public CommunicationInterface
	{
	public:

		CommunicationInterfaceCentral(const LoggingPtr& log, const int number_of_threads = 0);

		void set_coordinator(const CoordinatorPtr& coordinator);
		void set_simulator(const SimulatorPtr& simulator);

		/*Register an agent.*/
		const bool register_agent(const AgentPtr& agent) override;
		/*De-register an agent.*/
		const bool deregister_agent(const AgentInfo& agent) override;
		/*Register coupling between agents.*/
		const bool register_coupling(const CouplingInfo& coupling) override;
		/*De-register coupling between agents.*/
		const bool deregister_coupling(const CouplingInfo& coupling) override;

		/*This function is called if a message arrived to de-register a coupling..*/
		const bool fromCommunication_deregistered_coupling(const CouplingInfo& coupling) override;

		/*Send number of neighbors to an agent.*/
		const bool send_numberOfNeighbors(const int number, const int from, const int to) override;
		/*Send an agent state to an agent.*/
		const bool send_agentState(const AgentState& state, const int from, const int to) override;
		/*Send desired agent state to an agent.*/
		const bool send_desiredAgentState(const AgentState& desired_state, const int from, const int to) override;
		/*Send coupling state to an agent.*/
		const bool send_couplingState(const CouplingState& state, const int from, const int to) override;
		/*Send two coupling states to an agent.*/
		const bool send_couplingState(const CouplingState& state, const CouplingState& state2, const int from, const int to) override;
		/*Send multiplier and penalty states to an agent.*/
		const bool send_multiplierState(const MultiplierState& state, PenaltyState penalty, const int from, const int to) override;
		/*Send convergence flag to the coordinator.*/
		const bool send_convergenceFlag(const bool converged, const int from) override;

		/*Configure optimization.*/
		const bool configure_optimization(const OptimizationInfo& info) override;

		/*Trigger a step of the ADMM algorithm.*/
		const bool trigger_step(const ADMMStep& step) override;
		/*Trigger simulation.*/
		void trigger_simulation(const std::string& Integrator, const typeRNum dt) override;

		/*Return the agent state of an agent.*/
		const AgentStatePtr get_agentState_from_agent(const int agentId) const override;
		/*Return the desired agent state of an agent.*/
		const AgentStatePtr get_desiredAgentState_from_agent(const int agentId) const override;
		/*Returns a map with agent infos.*/
		const std::shared_ptr< std::map<unsigned int, AgentInfoPtr > > get_agentInfo_from_coordinator() const override;
		/*Return a map with sending neighbors.*/
		const std::shared_ptr< std::map< unsigned int, std::vector< CouplingInfoPtr > > > get_sendingNeighbors_from_coordinator() const override;
		/*Return a map with receiving neighbors.*/
		const std::shared_ptr< std::map< unsigned int, std::vector< CouplingInfoPtr > > > get_receivingNeighbors_from_coordinator() const override;
		/*Return the agent model of an agent.*/
		const AgentModelPtr get_agentModel(const int agentId) const override;
		/*Return all coupling models of an agent.*/
		const std::shared_ptr< std::map<int, CouplingModelPtr> > get_couplingModels_from_agent(const int agentId) const override;

		/*Set simulated state of an agent.*/
		void set_simulatedState_for_agent
		(
			const int agentId, 
			const std::vector<typeRNum>& new_state,
			const typeRNum dt, 
			const typeRNum t0,
			const typeRNum cost
		) override;
		/*Returns the current solution of an agent.*/
		const SolutionPtr get_solution(const unsigned int agent_id) const override;
		/*Return the current solution of a set of agents.*/
		const std::vector< SolutionPtr > get_solution(const std::string& agents) const override;
		/*Resets the solution of an agent.*/
		void reset_solution(const unsigned int agent_id) override;
		/*Resets the solution of a set of agents.*/
		void reset_solution(const std::string& agents) override;

		/*Wait for connections.*/
		void waitFor_connection(const int agents, const int couplings) override;
		/*Wait for a flag of the coordinator.*/
		void waitFor_flag_from_coordinator() override;
		/*Send flag to a set of agents.*/
		void send_flag_to_agents(const std::string& agents) const override;
		/*Send flag to a set of agents.*/
		void send_flag_to_agents(const std::vector<int>& agent_ids) const override;
		/*Send flag to an agent.*/
		void send_flag_to_agents(const int agent_id) const override;
		/*Set to passive mode.*/
		void set_passive() override;

		/*Cap the stored data to a number of data points.*/
		void cap_stored_data(const unsigned int data_points) override;
		/*Returns the number of agents.*/
		const unsigned int get_numberOfAgents() const override;

	private:
		CommunicationInfo comm_info_;
		CoordinatorPtr coordinator_;
		SimulatorPtr simulator_;
		std::vector<AgentPtr> agents_;

		LoggingPtr log_;
		std::ostringstream stream_;

		asio::io_service ioService_;
		asio::io_service::work work_;
		std::vector< std::thread > threads_;
		std::condition_variable condVar_triggerStep_;
		int counter_triggerStep_ = 0;
		std::mutex mutex_triggerStep_;
	};

}
