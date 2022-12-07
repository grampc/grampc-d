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

#include "grampcd/info/optimization_info.hpp"
#include "grampcd/comm/communication_interface.hpp"
#include "grampcd/optim/problem_description_local_sensi.hpp"
#include "grampcd/optim/solver_local.hpp"

namespace grampcd
{

    class SolverLocalSensi : public SolverLocal
    {
    public:
        SolverLocalSensi(Agent* agent, const OptimizationInfo& info, const LoggingPtr& log, const CommunicationInterfacePtr& communication_interface);

        /***********************************************
       Generic functions
       ************************************************/
        /*updates the agent States*/
        void update_agentStates() override;
        /*sends the agent states to receiving neighbors*/
        void send_agentStates() override;
       /*checks if the algorithm is converged*/
        const bool is_converged()  override;
        /*prints the debug costs to the solution file*/
        void print_debugCost() override;
        /*Sends the convergence flag*/
        void send_convergenceFlag() override;
        /*sends a flag if the agents finished its algorithm iterations*/
        void send_stoppedAlgFlag() override;

        /***********************************************
       Sensi functions
       ************************************************/
        /*evaluates the sensitivities for all receiving neighbors*/
        void update_sensiStates() override;
        /*Initializes the local sensi solver*/
        void initialize_Sensi() override;

    protected: 
        ProblemDescriptionLocalSensi sensi_problem_description_;
        SolverPtr solver_;
        LoggingPtr log_;
        Agent* agent_; 

        typeRNum cost_ = 0; 
        typeRNum previous_cost_ = 0;

        OptimizationInfo info_;
        CommunicationInterfacePtr communication_interface_;    
    };


}