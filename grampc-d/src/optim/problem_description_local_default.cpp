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

#include "grampcd/optim/problem_description_local_default.hpp"
#include "grampcd/optim/optim_util.hpp"

#include "grampcd/model/agent_model.hpp"
#include "grampcd/model/coupling_model.hpp"

#include <cmath>

namespace grampcd
{

        ProblemDescriptionLocalDefault::ProblemDescriptionLocalDefault(Agent* agent)
      : agent_(agent),
        Nx_( agent->get_agentModel()->get_Nxi() ),
        Nu_( agent->get_agentModel()->get_Nui() ),
        Ng_( agent->get_agentModel()->get_Ngi() ),
        Nh_( agent->get_agentModel()->get_Nhi() )
    {
        // create mapping from neighbor number j to indices of xj and uj within u
        int u_index = Nu_;
        int x_index = Nx_;

        // consider local copies of neighbors as controls
        for(const NeighborPtr& neighbor : agent_->get_neighbors())
        {
            const int j = neighbor->get_id();
            if( neighbor->is_sendingNeighbor() || neighbor->is_approximating() )
            {
                Nu_ += neighbor->get_Nxj();
                Nu_ += neighbor->get_Nuj();

                // adjust size of index vectors
                if(j >= u_index_xji_.size())
                {
                    u_index_uji_.resize(j + 1);
                    u_index_xji_.resize(j + 1);
                }

                u_index_xji_.at(j) = u_index;
                u_index += neighbor->get_Nxj();

                u_index_uji_.at(j) = u_index;
                u_index += neighbor->get_Nuj();
            }
        }

        // determine number of constraints
        Ng_ = agent_->get_agentModel()->get_Ngi();
        Nh_ = agent_->get_agentModel()->get_Nhi();
        for(const NeighborPtr& neighbor : agent_->get_neighbors())
        {
            if( neighbor->is_sendingNeighbor())
            {
                Ng_ += neighbor->get_couplingModel()->get_Ngij();
                Nh_ += neighbor->get_couplingModel()->get_Nhij();
            }
            if( agent->is_approximatingConstraints() )
            {
                // approximate neighbors agent constraints
                Ng_ += neighbor->get_agentModel()->get_Ngi();
                Nh_ += neighbor->get_agentModel()->get_Nhi();
                if(neighbor->is_receivingNeighbor())
                {
                    // approximate neighbors coupling constraints
                    Ng_ += neighbor->get_copied_couplingModel()->get_Ngij();
                    Nh_ += neighbor->get_copied_couplingModel()->get_Nhij();
                }
            }
        }
    }

    const std::vector<int>& ProblemDescriptionLocalDefault::get_u_index_uji() const
    {
        return u_index_uji_;
    }

    const std::vector<int>& ProblemDescriptionLocalDefault::get_u_index_xji() const
    {
        return u_index_xji_;
    }

    const int ProblemDescriptionLocalDefault::get_u_index_uji(int agent_id) const
    {
        return u_index_uji_.at(agent_id);
    }

    const int ProblemDescriptionLocalDefault::get_u_index_xji(int agent_id) const
    {
        return u_index_xji_.at(agent_id);
    }

    const int ProblemDescriptionLocalDefault::get_u_index_vji(int agent_id) const
    {
        return u_index_vji_.at(agent_id);
    }

    const int ProblemDescriptionLocalDefault::get_x_index_xji(int agent_id) const
    {
        return x_index_xji_.at(agent_id);
    }

    void ProblemDescriptionLocalDefault::ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT)
    {
        *Nx = Nx_;
        *Nu = Nu_;
        *Np = 0;
        *Ng = Ng_;
        *Nh = Nh_;
        *NgT = 0;
        *NhT = 0;
    }

    void ProblemDescriptionLocalDefault::ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
    {
	    MatSetScalar(out, 0, 1, Nx_);

        // agent dynamics f_i(x_i, u_i)
        agent_->get_agentModel()->ffct(out, t, x, u);

        for(const NeighborPtr& neighbor : agent_->get_sendingNeighbors())
        {
            const auto j = neighbor->get_id();

            // coupling dynamics f_{ij}(x_i, u_i, x_j, u_j)
            neighbor->get_couplingModel()->ffct(out, t, x, u, u + u_index_xji_[j], u + u_index_uji_[j]);
        }
    }

    void ProblemDescriptionLocalDefault::dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p)
    {
        MatSetScalar(out, 0, 1, Nx_);

        // agent dynamics \partial f_i(x_i, u_i) / \partial x_i
        agent_->get_agentModel()->dfdx_vec(out, t, x, u, vec);

        for(const NeighborPtr& neighbor : agent_->get_sendingNeighbors())
        {
            const auto j = neighbor->get_id();

            // coupling dynamics \partial f_{ij}(x_i, u_i, x_j, u_j) / \partial x_i
            neighbor->get_couplingModel()->dfdxi_vec(out, t, x, u, u + u_index_xji_[j], u + u_index_uji_[j], vec);
        }
    }

    void ProblemDescriptionLocalDefault::dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p)
    {
	    MatSetScalar(out, 0, 1, Nu_);

        // agent dynamics \partial f_i(x_i, u_i) / \partial u_i
        agent_->get_agentModel()->dfdu_vec(out, t, x, u, vec);

        for(const NeighborPtr& neighbor : agent_->get_sendingNeighbors())
        {
            const auto j = neighbor->get_id();

            // coupling dynamics \partial f_{ij}(x_i, u_i, x_j, u_j) / \partial u_i
            neighbor->get_couplingModel()->dfdui_vec(out, t, x, u, u + u_index_xji_[j], u + u_index_uji_[j], vec);

            // coupling dynamics \partial f_{ij}(x_i, u_i, x_j, u_j) / \partial x_j
            neighbor->get_couplingModel()->dfdxj_vec(out + u_index_xji_[j], t, x, u, u + u_index_xji_[j], u + u_index_uji_[j], vec);

            // coupling dynamics \partial f_{ij}(x_i, u_i, x_j, u_j) / \partial u_j
            neighbor->get_couplingModel()->dfduj_vec(out + u_index_uji_[j], t, x, u, u + u_index_xji_[j], u + u_index_uji_[j], vec);
        }
    }

    void ProblemDescriptionLocalDefault::lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
    {
        MatSetScalar(out, 0, 1, 1);

        // consider cost function l
        {
			typeRNum l_i = 0.0;
			typeRNum l_ij = 0.0;

            interpolateState(agent_->get_desiredAgentState(), t, desired_state_);

            // l_i( x_i, u_i )
            agent_->get_agentModel()->lfct(&l_i, t, x, u, &desired_state_.x_[0]);

            for (const auto& neighbor : agent_->get_sendingNeighbors())
            {
                const auto j = neighbor->get_id();

                neighbor->get_couplingModel()->lfct(&l_ij, t, x, u, u + u_index_xji_[j], u + u_index_uji_[j]);
            }

            // consider cost approximation
            if (!agent_->is_approximatingCost())
            {
                out[0] += l_i + l_ij;
            }
            else
			{
				// rescale agent cost
				out[0] += l_i / (1.0 + agent_->get_neighbors().size()) + l_ij / 2.0;

                // add the neighbors approximated cost
                for( const NeighborPtr& neighbor : agent_->get_neighbors() )
                {
                    const auto j = neighbor->get_id();
					typeRNum l_j = 0.0;
					typeRNum l_ji = 0.0;

                    interpolateState(neighbor->get_neighbors_desiredAgentState(), t, desired_state_);

                    neighbor->get_agentModel()->lfct( &l_j, t, u + u_index_xji_[j], u + u_index_uji_[j], &desired_state_.x_[0] );

                    if(neighbor->is_receivingNeighbor())
                        neighbor->get_copied_couplingModel()->lfct(&l_ji, t, u + u_index_xji_[j], u + u_index_uji_[j], x, u);

                    out[0] += l_j / ( 1.0 + neighbor->get_numberOfNeighbors() ) + l_ji / 2.0;
				}
            }
        }
        {
            interpolateState(agent_->get_couplingState(), t, couplingState_);
            interpolateState(agent_->get_multiplierState(), t, multiplierState_);
            interpolateState(agent_->get_penaltyState(), t, penaltyState_);

            const PenaltyState& penalty = agent_->get_penaltyState();

            // consistency constraints ( zx_i - x_i )
            for(unsigned int k = 0; k < agent_->get_Nxi(); ++k)
            {
                const typeRNum zx_min_x = couplingState_.z_x_[k] - x[k];
                out[0] += multiplierState_.mu_x_[k] * zx_min_x + 0.5 * penalty.rho_x_[k] * std::pow(zx_min_x, 2);
            }

            // consistency constraints ( zu_i - u_i )
            for(unsigned int k = 0; k < agent_->get_Nui(); ++k)
            {
                const typeRNum zu_min_u = couplingState_.z_u_[k] - u[k];
                out[0] += multiplierState_.mu_u_[k] * zu_min_u + 0.5 * penalty.rho_u_[k] * std::pow(zu_min_u, 2);
            }
        }

        // consistency constraints for neighbors
        for(const NeighborPtr& neighbor : agent_->get_sendingNeighbors())
        {
            const auto j = neighbor->get_id();

            interpolateState(neighbor->get_neighbors_couplingState(), t, couplingState_);
            interpolateState(neighbor->get_coupled_multiplierState(), t, multiplierState_);
            interpolateState(neighbor->get_coupled_penaltyState(), t, penaltyState_);

            // consistency constraints (zx_j - x_{ji})
            for(unsigned int k = 0; k < neighbor->get_Nxj(); ++k)
            {
                const typeRNum zx_min_x = couplingState_.z_x_[k] - (u + u_index_xji_[j])[k];
                out[0] += multiplierState_.mu_x_[k] * zx_min_x + 0.5 * penaltyState_.rho_x_[k] * std::pow(zx_min_x, 2);
            }

            // consistency constraints (zu_j - u_{ji})
            for(unsigned int k = 0; k < neighbor->get_Nuj(); ++k)
            {
                const typeRNum zu_min_u = couplingState_.z_u_[k] - (u + u_index_uji_[j])[k];
                out[0] += multiplierState_.mu_u_[k] * zu_min_u + 0.5 * penaltyState_.rho_u_[k] * std::pow(zu_min_u, 2);
            }
        }
    }

    void ProblemDescriptionLocalDefault::dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
    {
        MatSetScalar(out, 0, 1, Nx_);

        // consider cost function l
        {
			std::vector<typeRNum> l_i(Nx_, 0.0);
			std::vector<typeRNum> l_ij(Nx_, 0.0);

			interpolateState(agent_->get_desiredAgentState(), t, desired_state_);

            // \partial l_i( x_i, u_i ) / \partial x_i
            agent_->get_agentModel()->dldx(&l_i[0], t, x, u, &desired_state_.x_[0]);

            for (const auto& neighbor : agent_->get_sendingNeighbors())
            {
                const auto j = neighbor->get_id();

                neighbor->get_couplingModel()->dldxi(&l_ij[0], t, x, u, u + u_index_xji_[j], u + u_index_uji_[j]);
            }

            // consider neighbor approximation
            if (!agent_->is_approximatingCost())
            {
                for (unsigned int k = 0; k < Nx_; ++k)
                    out[k] += l_i[k] + l_ij[k];
            }
            else
            {
				std::vector<typeRNum> l_ji(Nx_, 0.0);
                for (const auto& neighbor : agent_->get_receivingNeighbors())
                {
					const auto j = neighbor->get_id();

                    neighbor->get_copied_couplingModel()->dldxj(&l_ji[0], t, u + u_index_xji_[j], u + u_index_uji_[j], x, u);
				}

				for (unsigned int k = 0; k < Nx_; ++k)
					out[k] += l_i[k] / (1.0 + agent_->get_neighbors().size()) + l_ij[k] / 2.0 + l_ji[k] / 2.0;
            }
        }

        // consider constraints
        {
            interpolateState(agent_->get_couplingState(), t, couplingState_);
            interpolateState(agent_->get_multiplierState(), t, multiplierState_);
            interpolateState(agent_->get_penaltyState(), t, penaltyState_);

            // consistency constraint ( zx_i - x_i )
            // derivative w.r.t. x_i
            for(unsigned int k = 0; k < agent_->get_Nxi(); ++k)
            {
                const typeRNum zx_min_x = couplingState_.z_x_[k] - x[k];
                out[k] += (multiplierState_.mu_x_[k] + penaltyState_.rho_x_[k] * zx_min_x) * (-1);
            }
        }
    }

    void ProblemDescriptionLocalDefault::dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
    {

        MatSetScalar(out, 0, 1, Nu_);

        // consider cost function l
        {
			std::vector<typeRNum> l_i(Nu_, 0.0);
			std::vector<typeRNum> l_ij(Nu_, 0.0);

            interpolateState(agent_->get_desiredAgentState(), t, desired_state_);

            // \partial l_i( x_i, u_i ) / \partial u_i
			agent_->get_agentModel()->dldu(&l_i[0], t, x, u, &desired_state_.x_[0]);

			for (const auto& neighbor : agent_->get_sendingNeighbors())
			{
				const auto j = neighbor->get_id();
                const auto& coupling_model = neighbor->get_couplingModel();

                coupling_model->dldui(&l_ij[0], t, x, u, u + u_index_xji_[j], u + u_index_uji_[j]);

				coupling_model->dldxj(&l_ij[0] + u_index_xji_[j], t, x, u, u + u_index_xji_[j], u + u_index_uji_[j]);

				coupling_model->dlduj(&l_ij[0] + u_index_uji_[j], t, x, u, u + u_index_xji_[j], u + u_index_uji_[j]);
			}

            // consider cost approximation
            if (!agent_->is_approximatingCost())
            {
                for (unsigned int k = 0; k < Nu_; ++k)
                    out[k] += l_i[k] + l_ij[k];
            }
            else
            {
                // rescale cost
				for (unsigned int k = 0; k < Nu_; ++k)
					out[k] += l_i[k] / (1.0 + agent_->get_neighbors().size()) + l_ij[k] / 2.0;

                // add approximated cost
                for( const auto& neighbor : agent_->get_neighbors() )
                {
                    interpolateState(neighbor->get_neighbors_desiredAgentState(), t, desired_state_);

                    const auto j = neighbor->get_id();
                    const auto Nxj = neighbor->get_Nxj();
					const auto Nuj = neighbor->get_Nuj();

					std::vector<typeRNum> l_j(Nu_, 0.0);
					std::vector<typeRNum> l_ji(Nu_, 0.0);

                    // consider local copies x_{ji} as control
                    neighbor->get_agentModel()->dldx( &l_j[0] + u_index_xji_[j], t, u + u_index_xji_[j], u + u_index_uji_[j], &desired_state_.x_[0] );

                    // consider local copies u_{ji} as control
                    neighbor->get_agentModel()->dldu( &l_j[0] + u_index_uji_[j], t, u + u_index_xji_[j], u + u_index_uji_[j], &desired_state_.x_[0] );

                    if (neighbor->is_receivingNeighbor())
                    {
						neighbor->get_copied_couplingModel()->dldxi(&l_ji[0] + u_index_xji_[j], t, u + u_index_xji_[j], u + u_index_uji_[j], x, u);

						neighbor->get_copied_couplingModel()->dldui(&l_ji[0] + u_index_uji_[j], t, u + u_index_xji_[j], u + u_index_uji_[j], x, u);

						neighbor->get_copied_couplingModel()->dlduj(&l_ji[0], t, u + u_index_xji_[j], u + u_index_uji_[j], x, u);
                    }

                    // rescale
                    for(unsigned int k = 0; k < Nu_; ++k)
                        out[k] += l_j[k] / (1.0 + neighbor->get_numberOfNeighbors()) + l_ji[k] / 2.0;
                }
            }
        }

        // consider own constraints
        {
            interpolateState(agent_->get_couplingState(), t, couplingState_);
            interpolateState(agent_->get_multiplierState(), t, multiplierState_);
            interpolateState(agent_->get_penaltyState(), t, penaltyState_);

            // consistency constraints ( zu_i - u_i )
            // derivative w.r.t. u_i
            for(unsigned int k = 0; k < agent_->get_Nui(); ++k)
            {
                const typeRNum zu_min_u = couplingState_.z_u_[k] - u[k];
                out[k] += (multiplierState_.mu_u_[k] + penaltyState_.rho_u_[k] * zu_min_u) * (-1);
            }
        }

        // consistency constraints for neighbors
        for(const NeighborPtr& neighbor : agent_->get_sendingNeighbors())
        {
            interpolateState(neighbor->get_neighbors_couplingState(), t, couplingState_);
            interpolateState(neighbor->get_coupled_multiplierState(), t, multiplierState_);
            interpolateState(neighbor->get_coupled_penaltyState(), t, penaltyState_);

            const auto j = neighbor->get_id();

            // consistency constraint ( zx_j - x_{ji} )
            // derivative w.r.t. x_{ji}
            for(unsigned int k = 0; k < couplingState_.z_x_.size(); ++k)
            {
                const typeRNum zx_min_x = couplingState_.z_x_[k] - (u + u_index_xji_[j])[k];
                (out + u_index_xji_[j])[k] += (multiplierState_.mu_x_[k] + penaltyState_.rho_x_[k] * zx_min_x) * (-1);
            }
            // consistency constraint ( zu_j - u_{ji} )
            // derivative w.r.t. u_{ji}
            for(unsigned int k = 0; k < couplingState_.z_u_.size(); ++k)
            {
                const typeRNum zu_min_u = couplingState_.z_u_[k] - (u + u_index_uji_[j])[k];
                (out + u_index_uji_[j])[k] += (multiplierState_.mu_u_[k] + penaltyState_.rho_u_[k] * zu_min_u) * (-1);
            }
        }
    }

    void ProblemDescriptionLocalDefault::Vfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes)
    {
        MatSetScalar(out, 0, 1, 1);
        interpolateState(agent_->get_desiredAgentState(), t, desired_state_);

        agent_->get_agentModel()->Vfct(out, t, x, &desired_state_.x_[0]);
    }

    void ProblemDescriptionLocalDefault::dVdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes)
    {
        MatSetScalar(out, 0, 1, Nx_);
        interpolateState(agent_->get_desiredAgentState(), t, desired_state_);

        agent_->get_agentModel()->dVdx(out, t, x, &desired_state_.x_[0]);
    }

    void ProblemDescriptionLocalDefault::gfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
    {
        MatSetScalar(out, 0, 1, Ng_);

        // equality constraints g_i(x_i, u_i) = 0
        agent_->get_agentModel()->gfct(out, t, x, u);
        int idx = agent_->get_agentModel()->get_Ngi();

        for(const NeighborPtr& neighbor : agent_->get_neighbors())
        {
            const auto j = neighbor->get_id();

            // equality constraints g_{ij}(x_i, u_i, x_j, u_j) = 0
           if( neighbor->is_sendingNeighbor() )
           {
               neighbor->get_couplingModel()->gfct(out + idx, t, x, u, u + u_index_xji_[j], u + u_index_uji_[j]);
               idx += neighbor->get_couplingModel()->get_Ngij();
           }

           if( agent_->is_approximatingConstraints() && neighbor->is_receivingNeighbor())
           {
               // equality constraints g_j(x_{ji}, u_{ji})
               neighbor->get_agentModel()->gfct(out + idx, t, u + u_index_xji_[j], u + u_index_uji_[j]);
               idx += neighbor->get_agentModel()->get_Ngi();

               // equality constraints g_{ji}(x_{ji}, u_{ji}, x_i, u_i)
               neighbor->get_copied_couplingModel()->gfct(out + idx, t, u + u_index_xji_[j], u + u_index_uji_[j], x, u);
               idx += neighbor->get_copied_couplingModel()->get_Ngij();
           }
        }
    }

    void ProblemDescriptionLocalDefault::dgdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
    {
        MatSetScalar(out, 0, 1, Nx_);

        // equality constraints \partial g_i(x_i, u_i) / \partial x_i
        agent_->get_agentModel()->dgdx_vec(out, t, x, u, vec);
        int idx = agent_->get_agentModel()->get_Ngi();

        for(const NeighborPtr& neighbor : agent_->get_neighbors())
        {
            const auto j = neighbor->get_id();

            // equality constraints \partial g_{ij}(x_i, u_i, x_j, u_j) / \partial x_i
            if(neighbor->is_sendingNeighbor())
            {
                neighbor->get_couplingModel()->dgdxi_vec(out, t, x, u, u + u_index_xji_[j], u + u_index_uji_[j], vec + idx);
                idx += neighbor->get_couplingModel()->get_Ngij();
            }

            if(agent_->is_approximatingConstraints() && neighbor->is_receivingNeighbor())
            {
                // equality constraints \partial g_{ji}(x_{ji}, u_{ji}, x_i, u_i) / \partial x_i
                neighbor->get_copied_couplingModel()->dgdxj_vec(out, t, u + u_index_xji_[j], u + u_index_uji_[j], x, u, vec + idx);
                idx += neighbor->get_copied_couplingModel()->get_Ngij();
            }
        }
    }

    void ProblemDescriptionLocalDefault::dgdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
    {
        MatSetScalar(out, 0, 1, Nu_);
        // equality constraints \partial g_i(x_i, u_i) / \partial u_i
        agent_->get_agentModel()->dgdu_vec(out, t, x, u, vec);
        int idx = agent_->get_agentModel()->get_Ngi();

        for(const NeighborPtr& neighbor : agent_->get_neighbors())
        {
            const auto j = neighbor->get_id();
            if( neighbor->is_sendingNeighbor() )
            {
                // equality constraints \partial g_{ij}(x_i, u_i, x_j, u_j) / \partial u_i
                neighbor->get_couplingModel()->dgdui_vec(out, t, x, u, u + u_index_xji_[j], u + u_index_uji_[j], vec + idx);

                // equality constraints \partial g_{ij}(x_i, u_i, x_j, u_j) / \partial x_j
                neighbor->get_couplingModel()->dgdxj_vec(out + u_index_xji_[j], t, x, u, u + u_index_xji_[j], u + u_index_uji_[j], vec + idx);

                // equality constraints \partial g_{ij}(x_i, u_i, x_j, u_j) / \partial u_j
                neighbor->get_couplingModel()->dgduj_vec(out + u_index_uji_[j], t, x, u, u + u_index_xji_[j], u + u_index_uji_[j], vec + idx);

                // increase index
                idx += neighbor->get_couplingModel()->get_Ngij();
            }
            if( agent_->is_approximatingConstraints() && neighbor->is_receivingNeighbor())
            {
                // equality constraints \partial g_j(x_{ji}, u_{ji}) / \partial x_{ji}
                neighbor->get_agentModel()->dgdx_vec(out + u_index_xji_[j], t, u + u_index_xji_[j], u + u_index_uji_[j], vec + idx);

                // equality constraints \partial g_j(x_{ji}, u_{ji}) / \partial u_{ji
                neighbor->get_agentModel()->dgdu_vec(out + u_index_uji_[j], t, u + u_index_xji_[j], u + u_index_uji_[j], vec + idx);
                idx += neighbor->get_agentModel()->get_Ngi();

                // equality constraints \partial g_{ji}(x_{ji}, u_{ji}, x_i, u_i) / \partial u_i
                neighbor->get_copied_couplingModel()->dgduj_vec(out, t, u + u_index_xji_[j], u + u_index_uji_[j], x, u, vec + idx);

                // equality constraints \partial g_{ji}(x_{ji}, u_{ji}, x_i, u_i) / \partial x_{ji}
                neighbor->get_copied_couplingModel()->dgdxi_vec(out + u_index_xji_[j], t, u + u_index_xji_[j], u + u_index_uji_[j], x, u, vec + idx);

                // equality constraints \partial g_{ji}(x_{ji}, u_{ji}, x_i, u_i) / \partial u_{ji}
                neighbor->get_copied_couplingModel()->dgdui_vec(out + u_index_uji_[j], t, u + u_index_xji_[j], u + u_index_uji_[j], x, u, vec + idx);

                // increase index
                idx += neighbor->get_copied_couplingModel()->get_Ngij();
            }
        }
    }

    void ProblemDescriptionLocalDefault::hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
    {
        MatSetScalar(out, 0, 1, Nh_);
        // inequality constraints h_i(x_i, u_i) <= 0
        agent_->get_agentModel()->hfct(out, t, x, u);
        int idx = agent_->get_agentModel()->get_Nhi();

        for(const NeighborPtr& neighbor : agent_->get_neighbors())
        {
            const auto j = neighbor->get_id();

            // inequality constraints h_{ij}(x_i, u_i, x_j, u_j) <= 0
            if(neighbor->is_sendingNeighbor())
            {
                neighbor->get_couplingModel()->hfct(out + idx, t, x, u, u + u_index_xji_[j], u + u_index_uji_[j]);
                idx += neighbor->get_couplingModel()->get_Nhij();
            }

            if(agent_->is_approximatingConstraints() && neighbor->is_receivingNeighbor())
            {
                // inequality constraint h_j(x_{ji}, u_{ji})
                neighbor->get_agentModel()->hfct(out + idx, t, u + u_index_xji_[j], u + u_index_uji_[j]);
                idx += neighbor->get_agentModel()->get_Nhi();

                // inequality constraint h_{ji}(x_{ji}, u_{ji}, x_i, u_i)
                neighbor->get_copied_couplingModel()->hfct(out + idx, t, u + u_index_xji_[j], u + u_index_uji_[j], x, u);
                idx += neighbor->get_copied_couplingModel()->get_Nhij();
            }
        }
    }

    void ProblemDescriptionLocalDefault::dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
    {
        MatSetScalar(out, 0, 1, Nx_);
        // inequality constraints \partial h_i(x_i, u_i) / \partial x_i
        agent_->get_agentModel()->dhdx_vec(out, t, x, u, vec);
        int idx = agent_->get_agentModel()->get_Nhi();

        for(const NeighborPtr& neighbor : agent_->get_neighbors())
        {
            const auto j = neighbor->get_id();

            // inequality constraints \partial h_{ij}(x_i, u_i, x_j, u_j) / \partial x_i
            if(neighbor->is_sendingNeighbor())
            {
                neighbor->get_couplingModel()->dhdxi_vec(out, t, x, u, u + u_index_xji_[j], u + u_index_uji_[j], vec + idx);
                idx += neighbor->get_couplingModel()->get_Nhij();
            }

            if(agent_->is_approximatingConstraints() && neighbor->is_receivingNeighbor())
            {
                // inequality constraint \partial h_{ji}(x_{ji}, u_{ji}, x_i, u_i) / \partial x_i
                neighbor->get_copied_couplingModel()->dhdxj_vec(out, t, u + u_index_xji_[j], u + u_index_uji_[j], x, u, vec + idx);
                idx += neighbor->get_copied_couplingModel()->get_Nhij();
            }
        }
    }

    void ProblemDescriptionLocalDefault::dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
    {
        MatSetScalar(out, 0, 1, Nu_);
        // inequality constraints \partial h_i(x_i, u_i) / \partial u_i
        agent_->get_agentModel()->dhdu_vec(out, t, x, u, vec);
        int idx = agent_->get_agentModel()->get_Nhi();

        for(const NeighborPtr& neighbor : agent_->get_neighbors())
        {
            const auto j = neighbor->get_id();

            if(neighbor->is_sendingNeighbor())
            {
                // inequality constraints \partial h_{ij}(x_i, u_i, x_j, u_j) / \partial u_i
                neighbor->get_couplingModel()->dhdui_vec(out, t, x, u, u + u_index_xji_[j], u + u_index_uji_[j], vec + idx);

                // inequality constraints \partial h_{ij}(x_i, u_i, x_j, u_j) / \partial x_j
                neighbor->get_couplingModel()->dhdxj_vec(out + u_index_xji_[j], t, x, u, u + u_index_xji_[j], u + u_index_uji_[j], vec + idx);

                // inequality constraints \partial h_{ij}(x_i, u_i, x_j, u_j) / \partial u_j
                neighbor->get_couplingModel()->dhduj_vec(out + u_index_uji_[j], t, x, u, u + u_index_xji_[j], u + u_index_uji_[j], vec + idx);

                // increase index
                idx += neighbor->get_couplingModel()->get_Nhij();
            }

            if(agent_->is_approximatingConstraints() && neighbor->is_receivingNeighbor())
            {
                // inequality constraint \partial h_j(x_{ji}, u_{ji}) / \partial x_{ji}
                neighbor->get_agentModel()->dhdx_vec(out + u_index_xji_[j], t, u + u_index_xji_[j], u + u_index_uji_[j], vec + idx);

                // inequality constraint \partial h_j(x_{ji}, u_{ji}) / \partial u_{ji}
                neighbor->get_agentModel()->dhdu_vec(out + u_index_uji_[j], t, u + u_index_xji_[j], u + u_index_uji_[j], vec + idx);

                // increase index
                idx += neighbor->get_agentModel()->get_Nhi();

                // inequality constraint \partial h_{ji}(x_{ji}, u_{ji}, x_i, u_i) / \partial u_i
                neighbor->get_copied_couplingModel()->dhduj_vec(out, t, u + u_index_xji_[j], u + u_index_uji_[j], x, u, vec + idx);

                // inequality constraint \partial h_{ji}(x_{ji}, u_{ji}, x_i, u_i) / \partial x_{ji}
                neighbor->get_copied_couplingModel()->dhdxi_vec(out + u_index_xji_[j], t, u + u_index_xji_[j], u + u_index_uji_[j], x, u, vec + idx);

                // inequality constraint \partial h_{ji}(x_{ji}, u_{ji}, x_i, u_i) / \partial u_{ji}
                neighbor->get_copied_couplingModel()->dhdui_vec(out + u_index_uji_[j], t, u + u_index_xji_[j], u + u_index_uji_[j], x, u, vec + idx);

                // increase index
                idx += neighbor->get_copied_couplingModel()->get_Nhij();
            }
        }
    }

}
