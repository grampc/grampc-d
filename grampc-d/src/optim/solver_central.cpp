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

#include "grampcd/optim/solver_central.hpp"

#include "grampcd/optim/optim_util.hpp"

#include "grampcd/agent/agent.hpp"
#include "grampcd/agent/neighbor.hpp"

#include "grampcd/model/agent_model.hpp"

namespace grampcd
{

    SolverCentral::SolverCentral(const std::vector<AgentPtr>& agents,
                                         const OptimizationInfo& info)
        : agents_(agents),
          problem_description_(agents_),
          solver_(new grampc::Grampc(&problem_description_)),
          optimizationInfo_(info)
    {
        configureSolver(solver_, info);
    }

    void SolverCentral::solve()
    {
        const auto Nx = solver_->getParameters()->Nx;
        const auto Nu = solver_->getParameters()->Nu;
        const unsigned int Nhor = solver_->getOptions()->Nhor;
        typeRNum t0 = 0.0;
        if(agents_.size() > 0)
            t0 = agents_[0]->get_agentState().t0_;

        // set initial state x0 and initial control solution
        for(const AgentPtr& agent : agents_)
        {
            const auto Nxi = agent->get_agentModel()->get_Nxi();
            const auto Nui = agent->get_agentModel()->get_Nui();
            const auto x_index = problem_description_.get_x_index(agent->get_id());
            const auto u_index = problem_description_.get_u_index(agent->get_id());

            const AgentState& state = agent->get_agentState();
            // set initial time
            solver_->setparam_real("t0", state.t0_);
            for(unsigned int j = 0; j < Nxi; ++j)
                solver_->getParameters()->x0[x_index + j] = state.x_[j];

            for(unsigned int i = 0; i < Nhor; ++i)
            {
                for(unsigned int j = 0; j < Nui; ++j)
                    solver_->getWorkspace()->u[i * Nu + u_index + j] = state.u_[i * Nui + j];
            }
        }

        // construct control limits
        std::vector<typeRNum> umin(solver_->getParameters()->Nu, -INF);
        std::vector<typeRNum> umax(solver_->getParameters()->Nu, INF);
        for(const AgentPtr& agent : agents_)
        {
            const int control_idx = problem_description_.get_u_index(agent->get_id());
            std::vector<typeRNum> uimin( agent->get_agentModel()->get_Nui(), -INF );
            std::vector<typeRNum> uimax( agent->get_agentModel()->get_Nui(), INF );
            agent->get_agentModel()->get_controlLimits(uimin, uimax);
            std::copy(uimin.begin(), uimin.end(), &umin[control_idx]);
            std::copy(uimax.begin(), uimax.end(), &umax[control_idx]);
        }
        solver_->setparam_real_vector("umin", &umin[0]);
        solver_->setparam_real_vector("umax", &umax[0]);
        solver_->setparam_real("t0", t0);

        // solve
        solver_->run();

        // update predicted state and control trajectories of agents
        for(const AgentPtr& agent : agents_)
        {
            const auto Nxi = agent->get_agentModel()->get_Nxi();
            const auto Nui = agent->get_agentModel()->get_Nui();
            const auto x_index = problem_description_.get_x_index(agent->get_id());
            const auto u_index = problem_description_.get_u_index(agent->get_id());

            AgentState state = agent->get_agentState();
            for(unsigned int i = 0; i < Nhor; ++i)
            {
                for(unsigned int j = 0; j < Nxi; ++j)
                    state.x_[i * Nxi + j] = solver_->getWorkspace()->x[i * Nx + x_index + j];

                for(unsigned int j = 0; j < Nui; ++j)
                    state.u_[i * Nui + j] = solver_->getWorkspace()->u[i * Nu + u_index + j];
            }
            agent->set_agentState(state);
        }
    }

    const std::vector<AgentPtr>& SolverCentral::get_agents() const
    {
        return agents_;
    }

    const bool SolverCentral::is_converged() const
    {
        return optimizationInfo_.GRAMPC_MaxGradIter_ != solver_->getSolution()->iter[optimizationInfo_.GRAMPC_MaxMultIter_ - 1];
    }
}
