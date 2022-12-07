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

#include"grampcd/agent/sync_step_selector.hpp"

#include "grampcd/optim/solver_local_admm.hpp"
#include "grampcd/optim/solver_local_sensi.hpp"

#include "grampcd/agent/agent.hpp"

namespace grampcd
{
    SyncStepSelector::SyncStepSelector(SolverLocalPtr& local_solver, Agent* agent)
        :
        agent_(agent),
        local_solver_(local_solver)
    {

    }


    void SyncStepSelector::execute_algStep(const AlgStep& step)
    {
        switch (step)
        {
        /*****************************
        ADMM STEPS
        ******************************/
        case AlgStep::ADMM_INITIALIZE:
        {
            local_solver_->initialize_ADMM();
            break;
        }

        case AlgStep::ADMM_UPDATE_AGENT_STATE:
        {
            local_solver_->update_agentStates();
            break;
        }

        case AlgStep::ADMM_UPDATE_COUPLING_STATE:
        {
            local_solver_->update_couplingStates();
            break;
        }

        case AlgStep::ADMM_UPDATE_MULTIPLIER_STATE:
        {
            local_solver_->update_multiplierStates();
            break;
        }

        case AlgStep::ADMM_SEND_AGENT_STATE:
        {
            local_solver_->send_agentStates();
            break;
        }

        case AlgStep::ADMM_SEND_COUPLING_STATE:
        {
            local_solver_->send_couplingStates();

            // increase ADMM iterations
            ++admmIter_;
            break;
        }

        case AlgStep::ADMM_SEND_MULTIPLIER_STATE:
        {
            local_solver_->send_multiplierStates();
            break;
        }
        
        /*****************************
       Sensi STEPS
       ******************************/

        case AlgStep::SENSI_INITIALIZE:
        {
            local_solver_->initialize_Sensi();
            break;
        }

        case AlgStep::SENSI_UPDATE_AGENT_STATE:
        {
            local_solver_->update_agentStates();
            break;
        }

        case AlgStep::SENSI_UPDATE_SENSI_STATE:
        {
            local_solver_->update_sensiStates();
            break;
        }

        case AlgStep::SENSI_SEND_AGENT_STATE:
        {
            local_solver_->send_agentStates();

            // increase sensi Iterations 
            ++sensiIter_;
            break;
        }
       /*****************************
      GENERIC STEPS
      ******************************/
        case AlgStep::GEN_SEND_CONVERGENCE_FLAG:
        {
            local_solver_->send_convergenceFlag();
            break;
        }

        case AlgStep::GEN_PRINT:
        {
            local_solver_->print_debugCost();
            break;
        }

        default:
            break;
        }
    }

    
    const unsigned int SyncStepSelector::get_algIter()
    {
        if (agent_->get_optimizationInfo().COMMON_Solver_ == "ADMM")
            return admmIter_;
        else if (agent_->get_optimizationInfo().COMMON_Solver_ == "Sensi")
            return sensiIter_;
        else
            return 0;
    }


}


