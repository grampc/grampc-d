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

#pragma once

#include "grampcd/util/class_forwarding.hpp"

#include "grampcd/info/optimization_info.hpp"

#include <mutex>
#include <condition_variable>

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

        /*************************************************************************
        coordination sensitivity-based algorithm
        *************************************************************************/
        /* Initialize agents for sensitivity-based algorithm */
        void initialize_sensi(const OptimizationInfo& oi);

        /* Solve optimization problem using sensitivity-based algorithm */
        const bool solve_sensi(int iterations = 1);

        /*************************************************************************
       Algorithm generic steps 
       *************************************************************************/

        /* Advance all agents to next sampling step */
        void trigger_simulation(const std::string& Integrator, typeRNum dt) const;

        /*Returns optimization info.*/
        const OptimizationInfo& get_optimizationInfo() const;

        /* Received convergence flag from agent */
        void fromCommunication_received_convergenceFlag(bool converged, int from);

        /*received flag of agent which executed all algorithm steps*/
        void fromCommunication_received_flagStoppedAlg(bool flag, int from);

    private:
	    std::map<unsigned int, AgentInfoPtr > agents_;
        CommunicationInterfacePtr communication_interface_;
        std::map< unsigned int, std::vector< CouplingInfoPtr > > sending_neighbors_;
        std::map< unsigned int, std::vector< CouplingInfoPtr > > receiving_neighbors_;
        // Algorithm variables 
        std::vector< int> agents_thatStoppedAlg_;
        std::vector< int> agents_thatConverged_;
        std::mutex mutex_stop_alg_;
        std::condition_variable cond_var_stop_alg_;
        bool alg_converged_ = false;

        LoggingPtr log_;

       
        bool simulation_ = false;
        OptimizationInfo optimizationInfo_;
    };

}