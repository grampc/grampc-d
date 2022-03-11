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

    enum class ADMMStep;

    /**
     * @brief The communication interface is used by the coordinator and each agent
     * to communicate over the network. There is no direct interaction between agents
     * or the coordinator except over the communication interface.
     */
    class CommunicationInterface
    {
    public:
        virtual ~CommunicationInterface();

        /*Register an agent.*/
        virtual const bool register_agent(const AgentPtr& agent) = 0;
        /*De-register an agent.*/
	    virtual const bool deregister_agent(const AgentInfo& agent) = 0;
        /*Register coupling between agents.*/
        virtual const bool register_coupling(const CouplingInfo& coupling) = 0;
        /*De-register coupling between agents.*/
        virtual const bool deregister_coupling(const CouplingInfo& coupling) = 0;

        /*This function is called if a message arrived to de-register a coupling..*/
        virtual const bool fromCommunication_deregistered_coupling(const CouplingInfo& coupling) = 0;

        /*Send number of neighbors to an agent.*/
        virtual const bool send_numberOfNeighbors( const int number, const int from, const int to ) = 0;
        /*Send an agent state to an agent.*/
        virtual const bool send_agentState(const AgentState& state, const int from, const int to) = 0;
        /*Send desired agent state to an agent.*/
        virtual const bool send_desiredAgentState( const AgentState& desired_state, const int from, const int to ) = 0;
        /*Send coupling state to an agent.*/
        virtual const bool send_couplingState(const CouplingState& state, const int from, const int to) = 0;
        /*Send two coupling states to an agent.*/
        virtual const bool send_couplingState(const CouplingState& state, const CouplingState& state2, const int from, const int to) = 0;
        /*Send multiplier and penalty states to an agent.*/
        virtual const bool send_multiplierState(const MultiplierState& state, PenaltyState penalty, const int from, const int to) = 0;
        /*Send convergence flag to the coordinator.*/
        virtual const bool send_convergenceFlag(const bool converged, const int from) = 0;
        /*send the stopped admm flag to an agent*/
        virtual const bool send_flagStoppedAdmm(const bool flag, const int from, const int to) = 0;
        /*ssend the stoppped admm flag to the coordinator*/
        virtual const bool send_flagStoppedAdmm(const bool flag, const int from) = 0;

        /*Configure optimization.*/
        virtual const bool configure_optimization(const OptimizationInfo& info) = 0;

        /*Trigger a step of the ADMM algorithm.*/
        virtual const bool trigger_step(const ADMMStep& step) = 0;
        /*Trigger simulation.*/
        virtual void trigger_simulation(const std::string& Integrator, const typeRNum dt) = 0;

        /*Return the agent state of an agent.*/
        virtual const AgentStatePtr get_agentState_from_agent(const int agentId ) const = 0;
        /*Return the desired agent state of an agent.*/
        virtual const AgentStatePtr get_desiredAgentState_from_agent(const int agentId ) const = 0;
        /*Returns a map with agent infos.*/
        virtual const std::shared_ptr< std::map<unsigned int, AgentInfoPtr > > get_agentInfo_from_coordinator() const = 0;
        /*Return a map with sending neighbors.*/
        virtual const std::shared_ptr< std::map< unsigned int, std::vector< CouplingInfoPtr > > > get_sendingNeighbors_from_coordinator() const = 0;
        /*Return a map with receiving neighbors.*/
        virtual const std::shared_ptr< std::map< unsigned int, std::vector< CouplingInfoPtr > > > get_receivingNeighbors_from_coordinator() const = 0;
        /*Return the agent model of an agent.*/
        virtual const AgentModelPtr get_agentModel(const int agentId ) const = 0;
        /*Return all coupling models of an agent.*/
        virtual const std::shared_ptr< std::map<int, CouplingModelPtr> > get_couplingModels_from_agent(const int agentId ) const = 0;

        /*Set simulated state of an agent.*/
        virtual void set_simulatedState_for_agent
        ( 
            const int agentId, 
            const std::vector<typeRNum>& new_state, 
            const typeRNum dt, 
			const typeRNum t0,
			const typeRNum cost
        ) = 0;
        /*Returns the current solution of an agent.*/
	    virtual const SolutionPtr get_solution(const unsigned int agent_id) const = 0;
        /*Return the current solution of a set of agents.*/
	    virtual const std::vector< SolutionPtr > get_solution(const std::string& agents) const = 0;
        /*Resets the solution of an agent.*/
        virtual void reset_solution(const unsigned int agent_id) = 0;
        /*Resets the solution of a set of agents.*/
        virtual void reset_solution(const std::string& agents) = 0;

        /*Wait for connections.*/
	    virtual void waitFor_connection(const int agents, const int couplings) = 0;
        /*Wait for a flag of the coordinator.*/
	    virtual void waitFor_flag_from_coordinator() = 0;
        /*Send flag to a set of agents.*/
        virtual void send_flag_to_agents(const std::string& agents) const = 0;
        /*Send flag to a set of agents.*/
        virtual void send_flag_to_agents(const std::vector<int>& agent_ids) const = 0;
        /*Send flag to an agent.*/
        virtual void send_flag_to_agents(const int agent_id) const = 0;
        /*Set to passive mode.*/
	    virtual void set_passive() = 0;

        /*Cap the stored data to a number of data points.*/
	    virtual void cap_stored_data(const unsigned int data_points) = 0;
        /*Returns the number of agents.*/
        virtual const unsigned int get_numberOfAgents() const = 0;

    };

}