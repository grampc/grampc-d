/* This file is part of GRAMPC-D - (https://github.com/DanielBurk/GRAMPC-D.git)
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

#ifndef SIMULATOR_HPP
#define SIMULATOR_HPP

#include "dmpc/comm/communication_interface.hpp"
#include "dmpc/optim/solver_central.hpp"

#include "dmpc/util/types.hpp"

namespace dmpc
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

    LoggingPtr log_;
    typeRNum dt_ = 0.0;
    typeRNum t0_ = 0.0;
    std::string Integrator_ = "";
    std::map<unsigned int, AgentStatePtr > agentStates_ = std::map<unsigned int, AgentStatePtr >();
    CommunicationInterfacePtr communication_interface_;
};

typedef std::shared_ptr<Simulator> SimulatorPtr;

}

#endif // SIMULATOR_HPP
