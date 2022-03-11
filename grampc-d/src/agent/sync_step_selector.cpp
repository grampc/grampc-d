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

#include"grampcd/agent/sync_step_selector.hpp"

#include "grampcd/optim/solver_local.hpp"

namespace grampcd
{
    SyncStepSelector::SyncStepSelector(SolverLocalPtr& local_solver)
        :
        local_solver_(local_solver)

    {

    }


    void SyncStepSelector::execute_admmStep(const ADMMStep& step)
    {
        switch (step)
        {
        case ADMMStep::INITIALIZE:
        {
            local_solver_->initialize_ADMM();
            break;
        }

        case ADMMStep::UPDATE_AGENT_STATE:
        {
            local_solver_->update_agentStates();
            break;
        }

        case ADMMStep::UPDATE_COUPLING_STATE:
        {
            local_solver_->update_couplingStates();
            break;
        }

        case ADMMStep::UPDATE_MULTIPLIER_STATE:
        {
            local_solver_->update_multiplierStates();
            break;
        }

        case ADMMStep::SEND_AGENT_STATE:
        {
            local_solver_->send_agentStates();
            break;
        }

        case ADMMStep::SEND_COUPLING_STATE:
        {
            local_solver_->send_couplingStates();

            // increase ADMM iterations
            ++admmIter_;
            break;
        }

        case ADMMStep::SEND_MULTIPLIER_STATE:
        {
            local_solver_->send_multiplierStates();
            break;
        }

        case ADMMStep::SEND_CONVERGENCE_FLAG:
        {
            local_solver_->send_convergenceFlag();
            break;
        }

        case ADMMStep::PRINT:
        {
            local_solver_->print_debugCost();
            break;
        }

        default:
            break;
        }
    }

    
    const unsigned int SyncStepSelector::get_admmIter()
    {
        return admmIter_;
    }


}


