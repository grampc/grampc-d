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

#include"grampcd/agent/async_step_selector.hpp"

#include "grampcd/optim/solver_local.hpp"
#include "grampcd/agent/agent.hpp"

namespace grampcd
{
	AsyncStepSelector::AsyncStepSelector(SolverLocalPtr& local_solver, Agent* agent)
		:
		local_solver_(local_solver),
        agent_(agent),
        previous_step_(ADMMStep::UPDATE_MULTIPLIER_STATE),
        current_step_(ADMMStep::UPDATE_AGENT_STATE),
        next_step_(ADMMStep::UPDATE_COUPLING_STATE),
        flagStopAdmm_(true)
	{

	}

	void AsyncStepSelector::execute_admmStep(const ADMMStep& step)
	{
        switch (step)
        {
        case ADMMStep::INITIALIZE:
        {
            admmIter_ = 0;
            local_solver_->initialize_ADMM();
            set_flagStopAdmm(false);

            // set steps
            previous_step_ = ADMMStep::UPDATE_MULTIPLIER_STATE;
            current_step_ = ADMMStep::INITIALIZE;
            next_step_ = ADMMStep::UPDATE_AGENT_STATE;
            break;
        }

        case ADMMStep::START_ASYNC_ADMM:
        {
            current_step_ = ADMMStep::START_ASYNC_ADMM;

            break;
        }

        case ADMMStep::UPDATE_AGENT_STATE:
        {   
            // set current and previous Step
            current_step_ = step;
            previous_step_ = ADMMStep::UPDATE_MULTIPLIER_STATE;

            // stopping criteria
            if (get_flagStopAdmm())
                break;
            if (next_step_ != ADMMStep::UPDATE_AGENT_STATE)
                break;
            if (!check_delays(current_step_))
                break;

            // set next step
            next_step_ = ADMMStep::UPDATE_COUPLING_STATE;

            // minimize local OCP w.r.t. local variables
            local_solver_->update_agentStates();

            // increase the relative delay of the neighbors 
            agent_->increase_all_delays(previous_step_);

            // update debug cost 
            if (agent_->get_optimizationInfo().ADMM_DebugCost_)
                local_solver_->print_debugCost();

            // send agent state
            local_solver_->send_agentStates();
            break;
        }

        case ADMMStep::UPDATE_COUPLING_STATE:
        {
            // set previous and current step
            current_step_ = step;
            previous_step_ = ADMMStep::UPDATE_AGENT_STATE;
            
            // stopping criteria
            if (get_flagStopAdmm())
                break;
            if (next_step_ != ADMMStep::UPDATE_COUPLING_STATE)
                break;
            if (!check_delays(current_step_))
                break;

            // set next step 
            next_step_ = ADMMStep::UPDATE_MULTIPLIER_STATE;

            local_solver_->update_couplingStates();

            // increase the relative delay of the neighbors 
            agent_->increase_all_delays(previous_step_);

            // send coupling state
            local_solver_->send_couplingStates();

            break;
        }

        case ADMMStep::UPDATE_MULTIPLIER_STATE:
        {   
            // set previous and current step
            current_step_ = step;
            previous_step_ = ADMMStep::UPDATE_COUPLING_STATE;
           

            // stopping criteria
            if (get_flagStopAdmm())
                break;
            if (next_step_ != ADMMStep::UPDATE_MULTIPLIER_STATE)
                break;
            if (!check_delays(current_step_))
                break;
                

            // set next step
            next_step_ = ADMMStep::UPDATE_AGENT_STATE;

            local_solver_->update_multiplierStates();

            // increase the relative delay of the neighbors 
            agent_->increase_all_delays(previous_step_);

            // increase admm iterations 
            ++admmIter_;

            // check if maximum iterations have been reached and set flag if so 
            check_admmIterations();

            // check if convergence has been achieved
            check_convergence();

            // send multiplier state
            local_solver_->send_multiplierStates();

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

        // check how to continue
        check_continuation();
    }

   void AsyncStepSelector::set_flagStopAdmm(const bool flag)
    {
        flagStopAdmm_ = flag;
    }

    const bool AsyncStepSelector::get_flagStopAdmm()
    {
        return flagStopAdmm_;
    }

    const unsigned int AsyncStepSelector::get_admmIter()
    {
        return admmIter_;

    }

    const bool AsyncStepSelector::check_delays(const ADMMStep& step)
    {
        switch (step)
        {

        case grampcd::ADMMStep::UPDATE_AGENT_STATE:
            if (agent_->get_delay_sending_neighbors(ADMMStep::UPDATE_MULTIPLIER_STATE) > check_criterionRegardingDelay())
                return false;
            break;

        case grampcd::ADMMStep::SEND_AGENT_STATE:
            return false;
            break;

        case grampcd::ADMMStep::UPDATE_COUPLING_STATE:
            if (agent_->get_delay_recieving_neighbors(ADMMStep::UPDATE_AGENT_STATE) > check_criterionRegardingDelay())
                return false;
            break;

        case grampcd::ADMMStep::SEND_COUPLING_STATE:
            return false;
            break;

        case grampcd::ADMMStep::UPDATE_MULTIPLIER_STATE:
            if (agent_->get_delay_sending_neighbors(ADMMStep::UPDATE_COUPLING_STATE) > check_criterionRegardingDelay())
                return false;
            break;

        case grampcd::ADMMStep::SEND_MULTIPLIER_STATE:
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
        if (!get_flagStopAdmm() && current_step_!=ADMMStep::INITIALIZE )
        {
          if (check_delays(next_step_))
                execute_admmStep(next_step_);
        }

        // set step after initializing
        if (current_step_ == ADMMStep::INITIALIZE)
            current_step_ = ADMMStep::START_ASYNC_ADMM;

    }

    void AsyncStepSelector::check_admmIterations()
    {
        // check if maximum iterations have been reached
        if (get_admmIter() == agent_->get_optimizationInfo().ADMM_maxIterations_)
        {
            set_flagStopAdmm(true);

            // send flag to neighbors and coordinator that the maximum ADMM iterations have been reached
            local_solver_->send_flagStoppedAdmm();

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