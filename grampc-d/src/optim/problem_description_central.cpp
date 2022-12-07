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

#include "grampcd/optim/optim_util.hpp"
#include "grampcd/optim/problem_description_central.hpp"

#include "grampcd/agent/agent.hpp"
#include "grampcd/agent/neighbor.hpp"

#include "grampcd/model/agent_model.hpp"
#include "grampcd/model/coupling_model.hpp"

namespace grampcd
{

        ProblemDescriptionCentral::ProblemDescriptionCentral(const std::vector<AgentPtr>& agents)
      : agents_(agents)
    {
        // create mapping from agent id to indices of xi and ui within x and u
        int x_index = 0;
        int u_index = 0;

        for(const auto& agent : agents_)
        {
            const auto i = agent->get_id();

            // adjust size of index vectors
            x_index_.resize(i + 1);
            u_index_.resize(i + 1);

            x_index_[i] = x_index;
            u_index_[i] = u_index;

            x_index += agent->get_agentModel()->get_Nxi();
            u_index += agent->get_agentModel()->get_Nui();
        }

        // final indices correspond to total number of states and controls
        Nx_ = x_index;
        Nu_ = u_index;

        // determine number of constraints
        Ng_ = 0;
        Nh_ = 0;
        for(const auto& agent : agents_)
        {
            Ng_ += agent->get_agentModel()->get_Ngi();
            Nh_ += agent->get_agentModel()->get_Nhi();
            for( const auto& neighbor : agent->get_sendingNeighbors() )
            {
                Ng_ += neighbor->get_couplingModel()->get_Ngij();
                Nh_ += neighbor->get_couplingModel()->get_Nhij();
            }
        }
    }

        const int ProblemDescriptionCentral::get_x_index(int agent_id) const
    {
        return x_index_[agent_id];
    }

        const int ProblemDescriptionCentral::get_u_index(int agent_id) const
    {
        return u_index_[agent_id];
    }

    void ProblemDescriptionCentral::ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT)
    {
        *Nx = Nx_;
        *Nu = Nu_;
        *Np = 0;
        *Ng = Ng_;
        *Nh = Nh_;
        *NgT = 0;
        *NhT = 0;
    }

    void ProblemDescriptionCentral::ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
    {
        MatSetScalar(out, 0, 1, Nx_);
        for(const AgentPtr& agent : agents_)
        {
            // agent dynamics
            const unsigned int i = agent->get_id();
            agent->get_agentModel()->ffct(out + x_index_[i], t, x + x_index_[i], u + u_index_[i]);

            for(const NeighborPtr& neighbor : agent->get_sendingNeighbors())
            {
                const unsigned int j = neighbor->get_id();

                // coupling dynamics f_{ij}(x_i, u_i, x_j, u_j)
                neighbor->get_couplingModel()->ffct(out + x_index_[i], t, x + x_index_[i], u + u_index_[i],
                                                                         x + x_index_[j], u + u_index_[j]);
            }
        }
    }

    void ProblemDescriptionCentral::dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p)
    {
	    MatSetScalar(out, 0, 1, Nx_);
        for(const AgentPtr& agent : agents_)
        {
            // agent dynamics
            const unsigned int i = agent->get_id();
            agent->get_agentModel()->dfdx_vec(out + x_index_[i], t, x + x_index_[i], u + u_index_[i], vec + x_index_[i]);

            for(const NeighborPtr& neighbor : agent->get_sendingNeighbors())
            {
                const unsigned int j = neighbor->get_id();

                // coupling dynamics \partial f_{ij}(x_i, u_i, x_j, u_j) / \partial x_i
                neighbor->get_couplingModel()->dfdxi_vec(out + x_index_[i], t, x + x_index_[i], u + u_index_[i],
                                                                              x + x_index_[j], u + u_index_[j], vec + x_index_[i]);
                // coupling dynamics \partial f_{ij}(x_i, u_i, x_j, u_j) / \partial x_j
                neighbor->get_couplingModel()->dfdxj_vec(out + x_index_[j], t, x + x_index_[i], u + u_index_[i],
                                                                              x + x_index_[j], u + u_index_[j], vec + x_index_[i]);
            }
        }
    }

    void ProblemDescriptionCentral::dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p)
    {
	    MatSetScalar(out, 0, 1, Nu_);
        for(const AgentPtr& agent : agents_)
        {
            // agent dynamics
            const unsigned int i = agent->get_id();
            agent->get_agentModel()->dfdu_vec(out + u_index_[i], t, x + x_index_[i], u + u_index_[i], vec + x_index_[i]);

            for(const NeighborPtr& neighbor : agent->get_sendingNeighbors())
            {
                const unsigned int j = neighbor->get_id();

                // coupling dynamics \partial f_{ij}(x_i, u_i, x_j, u_j) / \partial u_i
                neighbor->get_couplingModel()->dfdui_vec(out + u_index_[i], t, x + x_index_[i], u + u_index_[i],
                                                                              x + x_index_[j], u + u_index_[j], vec + x_index_[i]);
                // coupling dynamics \partial f_{ij}(x_i, u_i, x_j, u_j) / \partial u_j
                neighbor->get_couplingModel()->dfduj_vec(out + u_index_[j], t, x + x_index_[i], u + u_index_[i],
                                                                              x + x_index_[j], u + u_index_[j], vec + x_index_[i]);
            }
        }
    }

    void ProblemDescriptionCentral::lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
    {
	    MatSetScalar(out, 0, 1, 1);
        for(const AgentPtr& agent : agents_)
        {
            const unsigned int i = agent->get_id();

            interpolateState(agent->get_desiredAgentState(), t, desired_);
            agent->get_agentModel()->lfct(out, t, x + x_index_[i], u + u_index_[i], &desired_.x_[0]);

            for (const auto& neighbor : agent->get_sendingNeighbors())
            {
                const unsigned int j = neighbor->get_id();

                neighbor->get_couplingModel()->lfct(out, t, x + x_index_[i], u + u_index_[i], x + x_index_[j], u + u_index_[j]);
            }
        }
    }

    void ProblemDescriptionCentral::dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
    {
	    MatSetScalar(out, 0, 1, Nx_);
        for(const AgentPtr& agent : agents_)
        {
            const unsigned int i = agent->get_id();

            interpolateState(agent->get_desiredAgentState(), t, desired_);
			agent->get_agentModel()->dldx(out + x_index_[i], t, x + x_index_[i], u + u_index_[i], &desired_.x_[0]);

			for (const auto& neighbor : agent->get_sendingNeighbors())
			{
				const unsigned int j = neighbor->get_id();

				neighbor->get_couplingModel()->dldxi(out + x_index_[i], t, x + x_index_[i], u + u_index_[i], x + x_index_[j], u + u_index_[j]);

				neighbor->get_couplingModel()->dldxj(out + x_index_[j], t, x + x_index_[i], u + u_index_[i], x + x_index_[j], u + u_index_[j]);
			}
        }
    }

    void ProblemDescriptionCentral::dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
    {
	    MatSetScalar(out, 0, 1, Nu_);
        for(const AgentPtr& agent : agents_)
        {
            const unsigned int i = agent->get_id();

            interpolateState(agent->get_desiredAgentState(), t, desired_);
			agent->get_agentModel()->dldu(out + u_index_[i], t, x + x_index_[i], u + u_index_[i], &desired_.x_[0]);

			for (const auto& neighbor : agent->get_sendingNeighbors())
			{
				const unsigned int j = neighbor->get_id();

				neighbor->get_couplingModel()->dldui(out + u_index_[i], t, x + x_index_[i], u + u_index_[i], x + x_index_[j], u + u_index_[j]);

				neighbor->get_couplingModel()->dlduj(out + u_index_[j], t, x + x_index_[i], u + u_index_[i], x + x_index_[j], u + u_index_[j]);
			}
        }
    }

    void ProblemDescriptionCentral::Vfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes)
    {
	    MatSetScalar(out, 0, 1, 1);
        for(const AgentPtr& agent : agents_)
        {
            const unsigned int i = agent->get_id();

            interpolateState(agent->get_desiredAgentState(), t, desired_);
			agent->get_agentModel()->Vfct(out, t, x + x_index_[i], &desired_.x_[0]);

			for (const auto& neighbor : agent->get_sendingNeighbors())
			{
				const unsigned int j = neighbor->get_id();

				neighbor->get_couplingModel()->Vfct(out, t, x + x_index_[i], x + x_index_[j]);
			}
        }
    }

    void ProblemDescriptionCentral::dVdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes)
    {
	    MatSetScalar(out, 0, 1, Nx_);
        for(const AgentPtr& agent : agents_)
        {
            const unsigned int i = agent->get_id();

            interpolateState(agent->get_desiredAgentState(), t, desired_);
			agent->get_agentModel()->dVdx(out + x_index_[i], t, x + x_index_[i], &desired_.x_[0]);

			for (const auto& neighbor : agent->get_sendingNeighbors())
			{
				const unsigned int j = neighbor->get_id();

				neighbor->get_couplingModel()->dVdxi(out + x_index_[i], t, x + x_index_[i], x + x_index_[j]);

				neighbor->get_couplingModel()->dVdxj(out + x_index_[j], t, x + x_index_[i], x + x_index_[j]);
			}
        }
    }

    void ProblemDescriptionCentral::gfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
    {
	    MatSetScalar(out, 0, 1, Ng_);
        int idx = 0;

        for(const AgentPtr& agent : agents_)
        {
            const unsigned int i = agent->get_id();

            // equality constraints g_i(x_i, u_i) 
            agent->get_agentModel()->gfct(out + idx, t, x + x_index_[i], u + u_index_[i]);
            idx += agent->get_agentModel()->get_Ngi();

            for(const NeighborPtr& neighbor : agent->get_sendingNeighbors())
            {
                const unsigned int j = neighbor->get_id();

                // equality constraints g_{ij}(x_i, u_i, x_j, u_j) 
                neighbor->get_couplingModel()->gfct(out + idx, t, x + x_index_[i], u + u_index_[i],
                                                                 x + x_index_[j], u + u_index_[j]);
                idx += neighbor->get_couplingModel()->get_Ngij();
            }
        }
    }

    void ProblemDescriptionCentral::dgdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
    {
	    MatSetScalar(out, 0, 1, Nx_);
        int idx = 0;

        for(const AgentPtr& agent : agents_)
        {
            const unsigned int i = agent->get_id();

            // equality constraints \partial g_i(x_i, u_i) / \partial x_i
            agent->get_agentModel()->dgdx_vec(out + x_index_[i], t, x + x_index_[i], u + u_index_[i], vec + idx);
            idx += agent->get_agentModel()->get_Ngi();

            for(const NeighborPtr& neighbor : agent->get_sendingNeighbors())
            {
                const unsigned int j = neighbor->get_id();

                // equality constraints \partial g_{ij}(x_i, u_i, x_j, u_j) / \partial x_i
                neighbor->get_couplingModel()->dgdxi_vec(out + x_index_[i], t, x + x_index_[i], u + u_index_[i],
                                                                              x + x_index_[j], u + u_index_[j], vec + idx);
                // equality constraints \partial g_{ij}(x_i, u_i, x_j, u_j) / \partial x_j
                neighbor->get_couplingModel()->dgdxj_vec(out + x_index_[j], t, x + x_index_[i], u + u_index_[i],
                                                                              x + x_index_[j], u + u_index_[j], vec + idx);
                idx += neighbor->get_couplingModel()->get_Ngij();
            }
        }
    }

    void ProblemDescriptionCentral::dgdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
    {
	    MatSetScalar(out, 0, 1, Nu_);
        int idx = 0;
    
        for(const AgentPtr& agent : agents_)
        {
            const unsigned int i = agent->get_id();

            // equality constraints \partial g_i(x_i, u_i) / \partial u_i
            agent->get_agentModel()->dgdu_vec(out + u_index_[i], t, x + x_index_[i], u + u_index_[i], vec + idx);
            idx += agent->get_agentModel()->get_Ngi();

            for(const NeighborPtr& neighbor : agent->get_sendingNeighbors())
            {
                const unsigned int j = neighbor->get_id();

                // equality constraints \partial g_{ij}(x_i, u_i, x_j, u_j) / \partial u_i
                neighbor->get_couplingModel()->dgdui_vec(out + u_index_[i], t, x + x_index_[i], u + u_index_[i],
                                                                              x + x_index_[j], u + u_index_[j], vec + idx);
                // equality constraints \partial g_{ij}(x_i, u_i, x_j, u_j) / \partial u_j
                neighbor->get_couplingModel()->dgduj_vec(out + u_index_[j], t, x + x_index_[i], u + u_index_[i],
                                                                              x + x_index_[j], u + u_index_[j], vec + idx);
                idx += neighbor->get_couplingModel()->get_Ngij();
            }
        }
    }

    void ProblemDescriptionCentral::hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
    {
	    MatSetScalar(out, 0, 1, Nh_);
        int idx = 0;

        for(const AgentPtr& agent : agents_)
        {
            const unsigned int i = agent->get_id();

            // inequality constraints h_i(x_i, u_i) <= 0
            agent->get_agentModel()->hfct(out + idx, t, x + x_index_[i], u + u_index_[i]);
            idx += agent->get_agentModel()->get_Nhi();

            for(const NeighborPtr& neighbor : agent->get_sendingNeighbors())
            {
                const unsigned int j = neighbor->get_id();

                // inequality constraints h_{ij}(x_i, u_i, x_j, u_j) <= 0
                neighbor->get_couplingModel()->hfct(out + idx, t, x + x_index_[i], u + u_index_[i],
                                                                 x + x_index_[j], u + u_index_[j]);
                idx += neighbor->get_couplingModel()->get_Nhij();
            }
        }
    }

    void ProblemDescriptionCentral::dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
    {
	    MatSetScalar(out, 0, 1, Nx_);
        int idx = 0;

        for(const AgentPtr& agent : agents_)
        {
            const unsigned int i = agent->get_id();

            // inequality constraints \partial h_i(x_i, u_i) / \partial x_i
            agent->get_agentModel()->dhdx_vec(out + x_index_[i], t, x + x_index_[i], u + u_index_[i], vec + idx);
            idx += agent->get_agentModel()->get_Nhi();

            for(const NeighborPtr& neighbor : agent->get_sendingNeighbors())
            {
                const unsigned int j = neighbor->get_id();

                // inequality constraints \partial h_{ij}(x_i, u_i, x_j, u_j) / \partial x_i
                neighbor->get_couplingModel()->dhdxi_vec(out + x_index_[i], t, x + x_index_[i], u + u_index_[i],
                                                                              x + x_index_[j], u + u_index_[j], vec + idx);
                // inequality constraints \partial h_{ij}(x_i, u_i, x_j, u_j) / \partial x_j
                neighbor->get_couplingModel()->dhdxj_vec(out + x_index_[j], t, x + x_index_[i], u + u_index_[i],
                                                                              x + x_index_[j], u + u_index_[j], vec + idx);
                idx += neighbor->get_couplingModel()->get_Nhij();
            }
        }
    }

    void ProblemDescriptionCentral::dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
    {
	    MatSetScalar(out, 0, 1, Nu_);
        int idx = 0;

        for(const AgentPtr& agent : agents_)
        {
            const unsigned int i = agent->get_id();

            // inequality constraints \partial h_i(x_i, u_i) / \partial u_i
            agent->get_agentModel()->dhdu_vec(out + u_index_[i], t, x + x_index_[i], u + u_index_[i], vec + idx);
            idx += agent->get_agentModel()->get_Nhi();

            for(const NeighborPtr& neighbor : agent->get_sendingNeighbors())
            {
                const unsigned int j = neighbor->get_id();

                // inequality constraints \partial h_{ij}(x_i, u_i, x_j, u_j) / \partial u_i
                neighbor->get_couplingModel()->dhdui_vec(out + u_index_[i], t, x + x_index_[i], u + u_index_[i],
                                                                              x + x_index_[j], u + u_index_[j], vec + idx);
                // inequality constraints \partial h_{ij}(x_i, u_i, x_j, u_j) / \partial u_j
                neighbor->get_couplingModel()->dhduj_vec(out + u_index_[j], t, x + x_index_[i], u + u_index_[i],
                                                                              x + x_index_[j], u + u_index_[j], vec + idx);
                idx += neighbor->get_couplingModel()->get_Nhij();
            }
        }
    }

}
