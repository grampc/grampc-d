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

namespace grampcd
{

    /**
     * @brief abstract solver class for the distributed optimization 
     */
    class SolverLocal
    {
    public:
        virtual ~SolverLocal();

        /***********************************************
        Generic functions
        ************************************************/
        virtual void update_agentStates() = 0;
        virtual void send_agentStates() = 0;

        virtual const bool is_converged() = 0;
        virtual void print_debugCost() = 0;
        virtual void send_convergenceFlag() = 0;
        virtual void send_stoppedAlgFlag() = 0;

        /***********************************************
        ADMM functions 
        ************************************************/
       virtual void initialize_ADMM();
       virtual void update_couplingStates();
       virtual void send_couplingStates();
       virtual void update_multiplierStates();
       virtual void send_multiplierStates();
       virtual void send_numberofNeighbors();

        /***********************************************
        Sensi functions
        ************************************************/
       virtual void initialize_Sensi();
       virtual void update_sensiStates();
        
    };

}