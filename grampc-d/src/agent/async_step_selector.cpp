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

#include"grampcd/agent/async_step_selector.hpp"

#include "grampcd/optim/solver_local_admm.hpp"
#include "grampcd/optim/solver_local_sensi.hpp"

#include "grampcd/agent/agent.hpp"

namespace grampcd
{
	AsyncStepSelector::AsyncStepSelector(SolverLocalPtr& local_solver, Agent* agent)
		:
		local_solver_(local_solver),
        agent_(agent),
        flagStopAlg_(true)
	{

	}

	void AsyncStepSelector::execute_algStep(const AlgStep& step)
	{
        switch (step)
        {
         /*****************************
         ADMM STEPS
         ******************************/
        case AlgStep::ADMM_INITIALIZE:
        {
            admmIter_ = 0;
            local_solver_->initialize_ADMM();
            set_flagStopAlg(false);

            // set steps
            previous_step_ = AlgStep::ADMM_UPDATE_MULTIPLIER_STATE;
            current_step_ = AlgStep::ADMM_INITIALIZE;
            next_step_ = AlgStep::ADMM_UPDATE_AGENT_STATE;
            break;
        }

        case AlgStep::ADMM_START_ASYNC_ADMM:
        {
            current_step_ = AlgStep::ADMM_START_ASYNC_ADMM;

            break;
        }

        case AlgStep::ADMM_UPDATE_AGENT_STATE:
        {   
            // set current and previous Step
            current_step_ = step;
            previous_step_ = AlgStep::ADMM_UPDATE_MULTIPLIER_STATE;

            // stopping criteria
            if (get_flagStopAlg())
                break;
            if (next_step_ != AlgStep::ADMM_UPDATE_AGENT_STATE)
                break;
            if (!check_delays(current_step_))
                break;

            // set next step
            next_step_ = AlgStep::ADMM_UPDATE_COUPLING_STATE;

            // minimize local OCP w.r.t. local variables
            local_solver_->update_agentStates();

            // increase the relative delay of the neighbors 
            agent_->increase_all_delays(previous_step_);

            // update debug cost 
            if (agent_->get_optimizationInfo().COMMON_DebugCost_)
                local_solver_->print_debugCost();

            // send agent state
            local_solver_->send_agentStates();
            break;
        }

        case AlgStep::ADMM_UPDATE_COUPLING_STATE:
        {
            // set previous and current step
            current_step_ = step;
            previous_step_ = AlgStep::ADMM_UPDATE_AGENT_STATE;
            
            // stopping criteria
            if (get_flagStopAlg())
                break;
            if (next_step_ != AlgStep::ADMM_UPDATE_COUPLING_STATE)
                break;
            if (!check_delays(current_step_))
                break;

            // set next step 
            next_step_ = AlgStep::ADMM_UPDATE_MULTIPLIER_STATE;

            local_solver_->update_couplingStates();

            // increase the relative delay of the neighbors 
            agent_->increase_all_delays(previous_step_);

            // send coupling state
            local_solver_->send_couplingStates();

            break;
        }

        case AlgStep::ADMM_UPDATE_MULTIPLIER_STATE:
        {   
            // set previous and current step
            current_step_ = step;
            previous_step_ = AlgStep::ADMM_UPDATE_COUPLING_STATE;
           

            // stopping criteria
            if (get_flagStopAlg())
                break;
            if (next_step_ != AlgStep::ADMM_UPDATE_MULTIPLIER_STATE)
                break;
            if (!check_delays(current_step_))
                break;
                

            // set next step
            next_step_ = AlgStep::ADMM_UPDATE_AGENT_STATE;

            local_solver_->update_multiplierStates();

            // increase the relative delay of the neighbors 
            agent_->increase_all_delays(previous_step_);

            // increase admm iterations 
            ++admmIter_;

            // check if maximum iterations have been reached and set flag if so 
            check_algIterations();

            // check if convergence has been achieved
            check_convergence();

            // send multiplier state
            local_solver_->send_multiplierStates();

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
            sensiIter_ = 0;
            local_solver_->initialize_Sensi();
            set_flagStopAlg(false);

            // set steps
            previous_step_ = AlgStep::SENSI_UPDATE_AGENT_STATE;
            current_step_ = AlgStep::SENSI_INITIALIZE;
            next_step_ = AlgStep::SENSI_UPDATE_AGENT_STATE;
            break;
        }

        case AlgStep::SENSI_START_ASYNC_SENSI:
        {
            current_step_ = AlgStep::SENSI_START_ASYNC_SENSI;
            break;
        }

        case AlgStep::SENSI_UPDATE_AGENT_STATE:
        {

            // set previous and current step
            current_step_ = step;
            previous_step_ = AlgStep::SENSI_UPDATE_AGENT_STATE;

            // stopping criteria
            if (get_flagStopAlg())
                break;
            if (next_step_ != AlgStep::SENSI_UPDATE_AGENT_STATE)
                break;
            if (!check_delays(current_step_))
                break;

            // set next step
            next_step_ = AlgStep::SENSI_UPDATE_AGENT_STATE;

           // sensitivities do not have delays
           // calculate sensitivities for neighbors 
            local_solver_->update_sensiStates();

            // minimize local OCP w.r.t. local variables
            local_solver_->update_agentStates();

            // increase the relative delay of the neighbors 
            agent_->increase_all_delays(previous_step_);

            // update debug cost 
            if (agent_->get_optimizationInfo().COMMON_DebugCost_)
                local_solver_->print_debugCost();

            // increase admm iterations 
            ++sensiIter_;

            // check if maximum iterations have been reached and set flag if so 
            check_algIterations();

            // check if convergence has been achieved
            check_convergence();

            // directly send agent state
            local_solver_->send_agentStates();
            break;
        }


        case AlgStep::SENSI_SEND_AGENT_STATE:
        {
            local_solver_->send_agentStates();
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

        // check how to continue
        check_continuation();
    }

   void AsyncStepSelector::set_flagStopAlg(const bool flag)
    {
        flagStopAlg_ = flag;
    }

    const bool AsyncStepSelector::get_flagStopAlg()
    {
        return flagStopAlg_;
    }

    const unsigned int AsyncStepSelector::get_algIter()
    {
        if (agent_->get_optimizationInfo().COMMON_Solver_ == "ADMM")
            return admmIter_;
        else if (agent_->get_optimizationInfo().COMMON_Solver_ == "Sensi")
            return sensiIter_;
        else
            return 0;
      
    }

    const bool AsyncStepSelector::check_delays(const AlgStep& step)
    {
        switch (step)
        {

        case AlgStep::ADMM_UPDATE_AGENT_STATE:
            if (agent_->get_delay_sending_neighbors(AlgStep::ADMM_UPDATE_MULTIPLIER_STATE) > check_criterionRegardingDelay())
                return false;
            break;

        case AlgStep::ADMM_SEND_AGENT_STATE:
            return false;
            break;

        case AlgStep::ADMM_UPDATE_COUPLING_STATE:
            if (agent_->get_delay_receiving_neighbors(AlgStep::ADMM_UPDATE_AGENT_STATE) > check_criterionRegardingDelay())
                return false;
            break;

        case AlgStep::ADMM_SEND_COUPLING_STATE:
            return false;
            break;

        case AlgStep::ADMM_UPDATE_MULTIPLIER_STATE:
            if (agent_->get_delay_sending_neighbors(AlgStep::ADMM_UPDATE_COUPLING_STATE) > check_criterionRegardingDelay())
                return false;
            break;

        case AlgStep::ADMM_SEND_MULTIPLIER_STATE:
            return false;
            break;

        case AlgStep::SENSI_UPDATE_AGENT_STATE:
                if (agent_->get_delay_sending_neighbors(AlgStep::SENSI_UPDATE_AGENT_STATE) > check_criterionRegardingDelay())
                    return false;
                break;

        case AlgStep::SENSI_SEND_AGENT_STATE:
            return false; 
            break; 

        default:
            return false;
            break;

        }
        return true;

    }

    const int AsyncStepSelector::check_criterionRegardingDelay()
    {
        return agent_->get_optimizationInfo().ASYNC_Delay_;
    }


    void  AsyncStepSelector::check_continuation()
    {
        // only continue if there is no flag to stop or initialize step
        if (!get_flagStopAlg() && current_step_!=AlgStep::ADMM_INITIALIZE && current_step_ != AlgStep::SENSI_INITIALIZE)
        {
          if (check_delays(next_step_))
                execute_algStep(next_step_);
        }

        // set step after initializing
        if (current_step_ == AlgStep::ADMM_INITIALIZE)
            current_step_ = AlgStep::ADMM_START_ASYNC_ADMM;
        if (current_step_ == AlgStep::SENSI_INITIALIZE)
            current_step_ = AlgStep::SENSI_START_ASYNC_SENSI;

    }

    void AsyncStepSelector::check_algIterations()
    {
        // check if maximum iterations have been reached for ADMM and Sensi case 
        if ((agent_->get_optimizationInfo().COMMON_Solver_ == "ADMM" && get_algIter() == agent_->get_optimizationInfo().ADMM_maxIterations_ ) 
            || (agent_->get_optimizationInfo().COMMON_Solver_ == "Sensi" && get_algIter() == agent_->get_optimizationInfo().SENSI_maxIterations_))
        {
            set_flagStopAlg(true);

            // send flag to neighbors and coordinator that the maximum ADMM iterations have been reached
            local_solver_->send_stoppedAlgFlag();
        }
    }

    void AsyncStepSelector::check_convergence()
    {
        if (local_solver_->is_converged())
        {
            local_solver_->send_convergenceFlag();
        }
    }
}