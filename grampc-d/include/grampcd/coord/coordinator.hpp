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

#include "grampcd/info/optimization_info.hpp"

namespace grampcd
{

    /**
     * @brief The coordinator is the central node that knows all agents.
     * Its central purpose is to synchronize the actions of the agents.
     */
    class Coordinator
    {
    public:
        Coordinator(const CommunicationInterfacePtr& communication_interface, bool simulation, const LoggingPtr& log);

        /* Register agent in the network */
        const bool register_agent(const AgentInfo& agent_info);
        /* Register coupling in the network */
        const bool register_coupling(const CouplingInfo& coupling_info);

        /* Remove coupling from the network */
        const bool deregister_agent(const AgentInfo& agent_info);
        const bool deregister_coupling(CouplingInfoPtr coupling_info);

        /* Print network structure */
        void print_network() const;

        /* Return the number of agents in the network */
        const unsigned int get_numberOfAgents() const;
        /*Returns map with agent infos.*/
        const std::map<unsigned int, AgentInfoPtr >& get_agentInfos() const;
        /*Returns map with sending neighbors.*/
        const std::map< unsigned int, std::vector< CouplingInfoPtr > >& get_sendingNeighbors() const;
        /*Returns map with receiving neighbors.*/
        const std::map< unsigned int, std::vector< CouplingInfoPtr > >& get_receivingNeighbors() const;

        /*************************************************************************
         coordination of alternating direction method of multipliers
         *************************************************************************/

        /* Initialize agents for alternating direction method of multipliers (ADMM) */
        void initialize_ADMM(const OptimizationInfo& oi);

        /* Solve optimization problem using alternating direction method of multipliers (ADMM) */
        const bool solve_ADMM(int outer_iterations = 1, int inner_iterations = 1);

        /* Received convergence flag from agent */
        void fromCommunication_received_convergenceFlag(bool converged, int from);

        /* Advance all agents to next sampling step */
        void trigger_simulation(const std::string& Integrator, typeRNum dt) const;

        /*Returns optimization info.*/
        const OptimizationInfo& get_optimizationInfo() const;

    private:
	    std::map<unsigned int, AgentInfoPtr > agents_;
        CommunicationInterfacePtr communication_interface_;
        std::map< unsigned int, std::vector< CouplingInfoPtr > > sending_neighbors_;
        std::map< unsigned int, std::vector< CouplingInfoPtr > > receiving_neighbors_;

        LoggingPtr log_;

        bool ADMM_converged_ = false;
        bool simulation_ = false;
        OptimizationInfo optimizationInfo_;
    };

}