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

#include "grampcd/optim/solver_central.hpp"

#include "grampcd/util/class_forwarding.hpp"
#include "grampcd/util/logging.hpp"

namespace grampcd
{

    class Simulator
    {
    public:
        Simulator(const CommunicationInterfacePtr& communication_interface, const LoggingPtr& log);

        /*Simulate the overall system based on a distributed setup.*/
	    void distributed_simulation(const std::string& Integrator, typeRNum dt);
	    /*Simulate the overall system based on a centralized setup.*/
        void centralized_simulation(const SolverCentral *solver, const std::string& Integrator, typeRNum dt);

        /*Set t0*/
        void set_t0(typeRNum t0);

    private:
        void simulate();
        const typeRNum evaluate_cost
        (
            const unsigned int agent_id,
            const AgentModelPtr& agent_model, 
            const std::shared_ptr< std::map<int, CouplingModelPtr> >& coupling_models
        ) const;

        LoggingPtr log_;
        typeRNum dt_ = 0.0;
        typeRNum t0_ = 0.0;
        std::string Integrator_;
		std::map<unsigned int, AgentStatePtr > agentStates_;
		std::map<unsigned int, AgentStatePtr > desired_agentStates_;
        CommunicationInterfacePtr communication_interface_;
    };

}