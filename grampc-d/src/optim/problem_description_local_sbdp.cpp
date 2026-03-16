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

#include "grampcd/agent/agent.hpp"
#include "grampcd/agent/neighbor.hpp"

#include "grampcd/optim/problem_description_local_sbdp.hpp"
#include "grampcd/optim/optim_util.hpp"

#include "grampcd/model/agent_model.hpp"
#include "grampcd/model/coupling_model.hpp"

#include <cmath>

namespace grampcd
{

    ProblemDescriptionLocalSBDP::ProblemDescriptionLocalSBDP(Agent* agent)
      : agent_(agent),
        Nx_( agent->get_agentModel()->get_Nxi() ),
        Nu_( agent->get_agentModel()->get_Nui() ),
        Ng_( agent->get_agentModel()->get_Ngi() ),
        Nh_( agent->get_agentModel()->get_Nhi() )
    {

        // determine number of constraints
        Ng_ = agent_->get_agentModel()->get_Ngi();
        Nh_ = agent_->get_agentModel()->get_Nhi();
        for (const NeighborPtr& neighbor : agent_->get_sendingNeighbors())
        {         
                Ng_ += neighbor->get_couplingModel()->get_Ngij();
                Nh_ += neighbor->get_couplingModel()->get_Nhij();
        }
    }


    void ProblemDescriptionLocalSBDP::ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT)
    {
        *Nx = Nx_;
        *Nu = Nu_;
        *Np = 0;
        *Ng = Ng_;
        *Nh = Nh_;
        *NgT = 0;
        *NhT = 0;
    }

    void ProblemDescriptionLocalSBDP::ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
    {
        MatSetScalar(out, 0, 1, Nx_);

        // agent dynamics f_i(x_i, u_i)
        agent_->get_agentModel()->ffct(out, t, x, u);
        for (const NeighborPtr& neighbor : agent_->get_sendingNeighbors())
        {
            interpolateState(neighbor->get_neighbors_agentState(), t - agent_->get_agentState().t0_, neighbors_agentState_);

            // coupling dynamics f_{ij}(x_i, u_i, x_j, u_j)
            neighbor->get_couplingModel()->ffct(out, t, x, u, neighbors_agentState_.x_.data(), neighbors_agentState_.u_.data());
        } 
    }

    void ProblemDescriptionLocalSBDP::dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p)
    {
        MatSetScalar(out, 0, 1, Nx_);

        // agent dynamics \partial f_i(x_i, u_i) / \partial x_i
        agent_->get_agentModel()->dfdx_vec(out, t, x, u, vec);

        for (const NeighborPtr& neighbor : agent_->get_sendingNeighbors())
        {
            interpolateState(neighbor->get_neighbors_agentState(), t - agent_->get_agentState().t0_, neighbors_agentState_);

            // coupling dynamics \partial f_{ij}(x_i, u_i, x_j, u_j) / \partial x_i
            neighbor->get_couplingModel()->dfdxi_vec(out, t, x, u, neighbors_agentState_.x_.data(), neighbors_agentState_.u_.data(), vec);
        }
       
    }

    void ProblemDescriptionLocalSBDP::dfdu_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* vec, ctypeRNum* u, ctypeRNum* p)
    {
        MatSetScalar(out, 0, 1, Nu_);

        // agent dynamics \partial f_i(x_i, u_i) / \partial u_i
        agent_->get_agentModel()->dfdu_vec(out, t, x, u, vec);
        // add coupling part of neighbors
        for (const NeighborPtr& neighbor : agent_->get_sendingNeighbors())
        {
            interpolateState(neighbor->get_neighbors_agentState(), t - agent_->get_agentState().t0_, neighbors_agentState_);

            // coupling dynamics \partial f_{ij}(x_i, u_i, x_j, u_j) / \partial u_i
            neighbor->get_couplingModel()->dfdui_vec(out, t, x, u, neighbors_agentState_.x_.data(), neighbors_agentState_.u_.data(), vec);
        }

    }

    void ProblemDescriptionLocalSBDP::lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
    {
        MatSetScalar(out, 0, 1, 1);

        // consider cost function l
        {         
            typeRNum l_i = 0.0;
            typeRNum l_ij = 0.0;
            interpolateState(agent_->get_desiredAgentState(), t, desired_state_);

            // l_i( x_i, u_i )
            agent_->get_agentModel()->lfct(&l_i, t, x, u, &desired_state_.x_[0]);

            // l_ij (x_i,u_i,x_j,u_j)
            for (const auto& neighbor : agent_->get_sendingNeighbors())
            {
                interpolateState(neighbor->get_neighbors_agentState(), t, neighbors_agentState_);
            
                neighbor->get_couplingModel()->lfct(&l_ij, t, x, u, neighbors_agentState_.x_.data(), neighbors_agentState_.u_.data());
            }

            out[0] += l_i + l_ij;
        }

        // consider sensitivities and convexivity term 
        {
            // get old states and inputs
            interpolateState(agent_->get_previous_agentState(), t, previous_agentState_);

            // consider first-order sensitivities
            for (const NeighborPtr& neighbor : agent_->get_receivingNeighbors())
            {
                interpolateState(neighbor->get_SBDPState(), t, sbdpState_);
                // first-order sensitivity w.r.t. x 
                for (unsigned int k = 0; k < Nx_; ++k)
                {
                    out[0] += sbdpState_.psi_x_[k] * (x[k] - previous_agentState_.x_[k]);
                }

                // first order sensitivity w.r.t u 
                for (unsigned int k = 0; k < Nu_; ++k)
                {
                    out[0] += sbdpState_.psi_u_[k] * (u[k] - previous_agentState_.u_[k]);
                }

                // consider higher order terms of Tayler approximation
                if (agent_->get_optimizationInfo().SBDP_higherOrder_)
                {
                    // second-order sensitivity w.r.t x
                    for (unsigned int k = 0; k < Nx_; ++k)
                    {
                        for (unsigned int j = 0; j < Nx_; ++j)
                        {
                            out[0] += 0.5 * (x[k] - previous_agentState_.x_[k]) * sbdpState_.psi_xx_[k + Nx_ * j] * (x[j] - previous_agentState_.x_[j]);
                        }
                    }
                    // second-order sensitivity w.r.t u
                    for (unsigned int k = 0; k < Nu_; ++k)
                    {
                        for (unsigned int j = 0; j < Nu_; ++j)
                        {
                            out[0] += 0.5 * (u[k] - previous_agentState_.u_[k]) * sbdpState_.psi_uu_[k + Nu_ * j] * (u[j] - previous_agentState_.u_[j]);
                        }
                    }
                    // second-order sensitivity w.r.t x and u 
                    for (unsigned int k = 0; k < Nx_; ++k)
                    {
                        for (unsigned int j = 0; j < Nu_; ++j)
                        {
                            out[0] += (x[k] - previous_agentState_.x_[k]) * sbdpState_.psi_xu_[k + Nu_ * j] * (u[j] - previous_agentState_.u_[j]);
                        }
                    }
                }
                // consider convexification term 
                if (agent_->get_optimizationInfo().SBDP_ConvexivityTerm_)
                {
                    // convexivity term w.r.t. x
                    for (unsigned int k = 0; k < Nx_; ++k)
                    {
                        out[0] += (x[k] - previous_agentState_.x_[k]) * agent_->get_optimizationInfo().SBDP_ConvexFactor_x_[k] * (x[k] - previous_agentState_.x_[k]);
                    }

                    // convexivity term w.r.t. u
                    for (unsigned int k = 0; k < Nu_; ++k)
                    {
                        out[0] += (u[k] - previous_agentState_.u_[k]) * agent_->get_optimizationInfo().SBDP_ConvexFactor_u_[k] * (u[k] - previous_agentState_.u_[k]);
                    }
                }
            }
        }
    }

    void ProblemDescriptionLocalSBDP::dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
    {
        MatSetScalar(out, 0, 1, Nx_);

        // consider cost function l
        {
            std::vector<typeRNum> l_i(Nx_, 0.0);
            std::vector<typeRNum> l_ij(Nx_, 0.0);
            interpolateState(agent_->get_desiredAgentState(), t, desired_state_);

            // \partial l_i( x_i, u_i ) / \partial x_i
            agent_->get_agentModel()->dldx(&l_i[0], t, x, u, &desired_state_.x_[0]);

            // partial l_ij (x_i,u_i,x_j,u_j) / \partial x_i
            for (const auto& neighbor : agent_->get_sendingNeighbors())
            {
                interpolateState(neighbor->get_neighbors_agentState(), t, neighbors_agentState_);

                neighbor->get_couplingModel()->dldxi(&l_ij[0], t, x, u, neighbors_agentState_.x_.data(), neighbors_agentState_.u_.data());
            }

            for (unsigned int k = 0; k < Nx_; ++k)
                out[k] += l_i[k] + l_ij[k];
        }
        // consider sensitivities and convexivity term 
        { 
            interpolateState(agent_->get_previous_agentState(), t, previous_agentState_);
            // consider sensitivities
            for (const NeighborPtr& neighbor : agent_->get_receivingNeighbors())
            {
                interpolateState(neighbor->get_SBDPState(), t, sbdpState_);

                // first-order sensitivity w.r.t. x 
                for (unsigned int k = 0; k < Nx_; ++k)
                {
                    out[k] += sbdpState_.psi_x_[k];
                }
                // consider higher order terms of Tayler approximation
                if (agent_->get_optimizationInfo().SBDP_higherOrder_)
                {
                    // second-order sensitivity w.r.t.x
                    for (unsigned int k = 0; k < Nx_; ++k)
                    {
                        for (unsigned int j = 0; j < Nx_; ++j)
                        {
                            out[k] += (x[j] - previous_agentState_.x_[j]) * sbdpState_.psi_xx_[k * Nx_ + j];
                        }
                    }
                    // second-order sensitivity w.r.t x and u 
                    for (unsigned int k = 0; k < Nx_; ++k)
                    {
                        for (unsigned int j = 0; j < Nu_; ++j)
                        {
                            out[k] += sbdpState_.psi_xu_[k + Nu_ * j] * (u[j] - previous_agentState_.u_[j]);
                        }
                    }
                }
            }
            // consider convexification term 
            if (agent_->get_optimizationInfo().SBDP_ConvexivityTerm_)
            {   
                // convexivity term w.r.t. x
                for (unsigned int k = 0; k < Nx_; ++k)
                {
                    out[k] += 2 * agent_->get_optimizationInfo().SBDP_ConvexFactor_x_[k] * (x[k] - previous_agentState_.x_[k]);
                }
            }
        }
    }

    void ProblemDescriptionLocalSBDP::dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
    {
        MatSetScalar(out, 0, 1, Nu_);

        // consider cost function l
        {
            std::vector<typeRNum> l_i(Nu_, 0.0);
            std::vector<typeRNum> l_ij(Nu_, 0.0);

            interpolateState(agent_->get_desiredAgentState(), t, desired_state_);

            // \partial l_i( x_i, u_i ) / \partial u_i
            agent_->get_agentModel()->dldu(&l_i[0], t, x, u, &desired_state_.x_[0]);

            // partial l_ij (x_i,u_i,x_j,u_j) / \partial x_i
            for (const auto& neighbor : agent_->get_sendingNeighbors())
            {
                interpolateState(neighbor->get_neighbors_agentState(), t, neighbors_agentState_);

                neighbor->get_couplingModel()->dldui(&l_ij[0], t, x, u, neighbors_agentState_.x_.data(), neighbors_agentState_.u_.data());
            }

            for (unsigned int k = 0; k < Nu_; ++k)
                out[k] += l_i[k] + l_ij[k];;
        }

        // consider sensitivities and convexivity term 
        {
            interpolateState(agent_->get_previous_agentState(), t, previous_agentState_);
            // consider sensitivities
            for (const NeighborPtr& neighbor : agent_->get_receivingNeighbors())
            {
                interpolateState(neighbor->get_SBDPState(), t, sbdpState_);
                // sensitivity w.r.t. u 
                for (unsigned int k = 0; k < Nu_; ++k)
                {
                    out[k] += sbdpState_.psi_u_[k];
                }
                // consider higher order terms of Tayler approximation
                if (agent_->get_optimizationInfo().SBDP_higherOrder_)
                {
                    // second-order sensitivity w.r.t.u
                    for (unsigned int k = 0; k < Nu_; ++k)
                    {
                        for (unsigned int j = 0; j < Nu_; ++j)
                        {
                            out[k] += (u[j] - previous_agentState_.u_[j]) * sbdpState_.psi_uu_[k * Nu_ + j] ;
                        }
                    }
                    // second-order sensitivity w.r.t x and u 
                    for (unsigned int k = 0; k < Nx_; ++k)
                    {
                        for (unsigned int j = 0; j < Nu_; ++j)
                        {
                            out[j] += (x[k] - previous_agentState_.x_[k]) * sbdpState_.psi_xu_[k * Nu_ + j];
                        }
                    }
                }
            }
            // consider convexification term 
            if (agent_->get_optimizationInfo().SBDP_ConvexivityTerm_)
            {
                // convexivity term w.r.t. x
                for (unsigned int k = 0; k < Nu_; ++k)
                {
                    out[k] += 2 * agent_->get_optimizationInfo().SBDP_ConvexFactor_u_[k] * (u[k] - previous_agentState_.u_[k]);
                }
            }
        }
    }

    void ProblemDescriptionLocalSBDP::Vfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes)
    {
        MatSetScalar(out, 0, 1, 1);
        interpolateState(agent_->get_desiredAgentState(), t, desired_state_);

        agent_->get_agentModel()->Vfct(out, t, x, &desired_state_.x_[0]);

        // Directly consider coupled costs 
        for (const auto& neighbor : agent_->get_sendingNeighbors())
        {
            interpolateState(neighbor->get_neighbors_agentState(), t, neighbors_agentState_);
            neighbor->get_couplingModel()->Vfct(out, t, x, neighbors_agentState_.x_.data());
        }

        // consider sensitivity of coupled terminal cost 
        // get old states and inputs
        interpolateState(agent_->get_previous_agentState(), t, previous_agentState_);

        // consider sensitivities of coupled terminal costs
        for (const NeighborPtr& neighbor : agent_->get_receivingNeighbors())
        {
            // Terminal cost only invloves sensitivity w.r.t. x 
            for (unsigned int k = 0; k < agent_->get_Nxi(); ++k)
            {
                out[0] += neighbor->get_SBDPState().psi_V_[k] * (x[k] - previous_agentState_.x_[k]);

                // consider higher order terms of Taylor approximation
                if (agent_->get_optimizationInfo().SBDP_higherOrder_)
                {
                    for (unsigned int j = 0; j < agent_->get_Nxi(); ++j)
                    {
                        out[0] += 0.5 * (x[k] - previous_agentState_.x_[k]) * neighbor->get_SBDPState().psi_VV_[k + Nx_ * j] * (x[j] - previous_agentState_.x_[j]);
                    }
                }
            }
        }    
    }

    void ProblemDescriptionLocalSBDP::dVdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes)
    {
        MatSetScalar(out, 0, 1, Nx_);
        interpolateState(agent_->get_desiredAgentState(), t, desired_state_);

        agent_->get_agentModel()->dVdx(out, t, x, &desired_state_.x_[0]);

        // Directly consider coupled costs 
        for (const auto& neighbor : agent_->get_sendingNeighbors())
        {
            interpolateState(neighbor->get_neighbors_agentState(), t, neighbors_agentState_);
            neighbor->get_couplingModel()->dVdxi(out, t, x, neighbors_agentState_.x_.data());
        }

        // consider sensitivity of coupled terminal cost 
        // get old states and inputs
        interpolateState(agent_->get_previous_agentState(), t, previous_agentState_);

        // consider sensitivities of coupled terminal costs
        for (const NeighborPtr& neighbor : agent_->get_receivingNeighbors())
        {
            // Terminal cost only invloves sensitivity w.r.t. x 
            for (unsigned int k = 0; k < agent_->get_Nxi(); ++k)
            {
                out[k] += neighbor->get_SBDPState().psi_V_[k];

                // consider higher order terms of Tayler approximation
                if (agent_->get_optimizationInfo().SBDP_higherOrder_)
                {
                    for (unsigned int j = 0; j < agent_->get_Nxi(); ++j)
                    {
                        out[k] += (x[j] - previous_agentState_.x_[j]) * neighbor->get_SBDPState().psi_VV_[k * Nx_ + j];
                    }
                }
            }
        }      
    }

    void ProblemDescriptionLocalSBDP::gfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
    {
        MatSetScalar(out, 0, 1, Ng_);

        // equality constraints g_i(x_i, u_i) = 0
        agent_->get_agentModel()->gfct(out, t, x, u);
        int idx = agent_->get_agentModel()->get_Ngi();
        for (const NeighborPtr& neighbor : agent_->get_neighbors())
        {
            const auto j = neighbor->get_id();
            interpolateState(neighbor->get_neighbors_agentState(), t, neighbors_agentState_);

            // equality constraints g_{ij}(x_i, u_i, x_j, u_j) = 0
            if (neighbor->is_sendingNeighbor())
            {
                neighbor->get_couplingModel()->gfct(out + idx, t, x, u, neighbors_agentState_.x_.data(), neighbors_agentState_.u_.data());
                idx += neighbor->get_couplingModel()->get_Ngij();
            }
        }  
    }

    void ProblemDescriptionLocalSBDP::dgdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
    {
        MatSetScalar(out, 0, 1, Nx_);

        // equality constraints \partial g_i(x_i, u_i) / \partial x_i
        agent_->get_agentModel()->dgdx_vec(out, t, x, u, vec);
        int idx = agent_->get_agentModel()->get_Ngi();

        for (const NeighborPtr& neighbor : agent_->get_neighbors())
        {
            const auto j = neighbor->get_id();
            interpolateState(neighbor->get_neighbors_agentState(), t, neighbors_agentState_);

            // equality constraints \partial g_{ij}(x_i, u_i, x_j, u_j) / \partial x_i
            if (neighbor->is_sendingNeighbor())
            {
                neighbor->get_couplingModel()->dgdxi_vec(out, t, x, u, neighbors_agentState_.x_.data(), neighbors_agentState_.u_.data(), vec + idx);
                idx += neighbor->get_couplingModel()->get_Ngij();
            }
        }   
    }

    void ProblemDescriptionLocalSBDP::dgdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
    {
        MatSetScalar(out, 0, 1, Nu_);
        // equality constraints \partial g_i(x_i, u_i) / \partial u_i
        agent_->get_agentModel()->dgdu_vec(out, t, x, u, vec);
        int idx = agent_->get_agentModel()->get_Ngi();

        for (const NeighborPtr& neighbor : agent_->get_neighbors())
        {
            const auto j = neighbor->get_id();
            interpolateState(neighbor->get_neighbors_agentState(), t, neighbors_agentState_);

            if (neighbor->is_sendingNeighbor())
            {
                // equality constraints \partial g_{ij}(x_i, u_i, x_j, u_j) / \partial u_i
                neighbor->get_couplingModel()->dgdui_vec(out, t, x, u, neighbors_agentState_.x_.data(), neighbors_agentState_.u_.data(), vec + idx);

                // increase index
                idx += neighbor->get_couplingModel()->get_Ngij();
            }
        }      
    }

    void ProblemDescriptionLocalSBDP::hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
    {
        MatSetScalar(out, 0, 1, Nh_);
        // inequality constraints h_i(x_i, u_i) <= 0
        agent_->get_agentModel()->hfct(out, t, x, u);
        int idx = agent_->get_agentModel()->get_Nhi();

        for (const NeighborPtr& neighbor : agent_->get_neighbors())
        {
            const auto j = neighbor->get_id();
            interpolateState(neighbor->get_neighbors_agentState(), t, neighbors_agentState_);

            // inequality constraints h_{ij}(x_i, u_i, x_j, u_j) <= 0
            if (neighbor->is_sendingNeighbor())
            {
                neighbor->get_couplingModel()->hfct(out + idx, t, x, u, neighbors_agentState_.x_.data(), neighbors_agentState_.u_.data());
                idx += neighbor->get_couplingModel()->get_Nhij();
            }
        }  
    }

    void ProblemDescriptionLocalSBDP::dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
    {
        MatSetScalar(out, 0, 1, Nx_);
        // inequality constraints \partial h_i(x_i, u_i) / \partial x_i
        agent_->get_agentModel()->dhdx_vec(out, t, x, u, vec);
        int idx = agent_->get_agentModel()->get_Nhi();

        for (const NeighborPtr& neighbor : agent_->get_neighbors())
        {
            const auto j = neighbor->get_id();
            interpolateState(neighbor->get_neighbors_agentState(), t, neighbors_agentState_);

            // inequality constraints \partial h_{ij}(x_i, u_i, x_j, u_j) / \partial x_i
            if (neighbor->is_sendingNeighbor())
            {
                neighbor->get_couplingModel()->dhdxi_vec(out, t, x, u, neighbors_agentState_.x_.data(), neighbors_agentState_.u_.data(), vec + idx);
                idx += neighbor->get_couplingModel()->get_Nhij();
            }

        }
    
    }

    void ProblemDescriptionLocalSBDP::dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
    {
        MatSetScalar(out, 0, 1, Nu_);
        // inequality constraints \partial h_i(x_i, u_i) / \partial u_i
        agent_->get_agentModel()->dhdu_vec(out, t, x, u, vec);
        int idx = agent_->get_agentModel()->get_Nhi();

        for (const NeighborPtr& neighbor : agent_->get_neighbors())
        {
            const auto j = neighbor->get_id();
            interpolateState(neighbor->get_neighbors_agentState(), t, neighbors_agentState_);

            if (neighbor->is_sendingNeighbor())
            {
                // inequality constraints \partial h_{ij}(x_i, u_i, x_j, u_j) / \partial u_i
                neighbor->get_couplingModel()->dhdui_vec(out, t, x, u, neighbors_agentState_.x_.data(), neighbors_agentState_.u_.data(), vec + idx);

                // increase index
                idx += neighbor->get_couplingModel()->get_Nhij();
            }
        }
    }
}
