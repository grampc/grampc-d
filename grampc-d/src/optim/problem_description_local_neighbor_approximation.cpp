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

#include "grampcd/optim/optim_util.hpp"
#include "grampcd/optim/problem_description_local_neighbor_approximation.hpp"
#include "grampcd/optim/approximate_neighbor.hpp"

#include "grampcd/model/agent_model.hpp"
#include "grampcd/model/coupling_model.hpp"

namespace grampcd
{

    ProblemDescriptionLocalNeighborApproximation::ProblemDescriptionLocalNeighborApproximation(Agent* agent, const OptimizationInfo& optimization_info)
      : agent_(agent),
        Nx_(agent->get_agentModel()->get_Nxi()),
        Nu_(agent->get_agentModel()->get_Nui()),
        Ng_(agent->get_agentModel()->get_Ngi()),
        Nh_(agent->get_agentModel()->get_Nhi()),
        optimizationInfo_(optimization_info)
    {
        // create mapping from neighbor number j to indices of xj and uj within u
        int u_index = Nu_;
        int x_index = Nx_;

        // consider local copies of neighbors as controls
        for(const auto& neighbor : agent_->get_neighbors())
        {
            const int j = neighbor->get_id();

            Nu_ += neighbor->get_Nxj();
            Nu_ += neighbor->get_Nuj();
            Nx_ += neighbor->get_Nxj();

            // adjust size of index vectors
            if(j >= x_index_xji_.size())
            {
                u_index_uji_.resize(j + 1);
                u_index_vji_.resize(j + 1);
                x_index_xji_.resize(j + 1);
            }

            u_index_uji_.at(j) = u_index;
            u_index += neighbor->get_Nuj();

            u_index_vji_.at(j) = u_index;
            u_index += neighbor->get_Nxj();

            x_index_xji_.at(j) = x_index;
            x_index += neighbor->get_Nxj();
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

    const std::vector<int>& ProblemDescriptionLocalNeighborApproximation::get_x_index_xji() const
    {
        return x_index_xji_;
    }

    const std::vector<int>& ProblemDescriptionLocalNeighborApproximation::get_u_index_uji() const
    {
        return u_index_uji_;
    }

    const std::vector<int>& ProblemDescriptionLocalNeighborApproximation::get_u_index_vji() const
    {
        return u_index_vji_;
    }

    const int ProblemDescriptionLocalNeighborApproximation::get_x_index_xji(int agent_id) const
    {
        return x_index_xji_[agent_id];
    }

    const int ProblemDescriptionLocalNeighborApproximation::get_u_index_vji(int agent_id) const
    {
        return u_index_vji_[agent_id];
    }

    const int ProblemDescriptionLocalNeighborApproximation::get_u_index_uji(int agent_id) const
    {
        return u_index_uji_[agent_id];
    }

    void ProblemDescriptionLocalNeighborApproximation::ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT)
    {
        *Nx = Nx_;
        *Nu = Nu_;
        *Np = 0;
        *Ng = Ng_;
        *Nh = Nh_;
        *NgT = 0;
        *NhT = 0;
    }

    void ProblemDescriptionLocalNeighborApproximation::ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
    {
        // reset out
        MatSetScalar(out, 0, 1, Nx_);

        // consider agent dynamics f_i(x_i, u_i)
        agent_->get_agentModel()->ffct(out, t, x, u);

        for(const auto& neighbor : agent_->get_neighbors())
        {
            const auto j = neighbor->get_id();

		    // consider coupling dynamics f_{ij}(x_i, u_i, x_j, u_j)
            if( neighbor->is_sendingNeighbor() )
                neighbor->get_couplingModel()->ffct(out, t, x, u, x + x_index_xji_[j], u + u_index_uji_[j]);

            // approximate neighbors dynamics
            neighbor->get_neighborApproximation()->ffct(out, t, x, u, x_index_xji_, u_index_uji_, u_index_vji_);
        }
    }

    void ProblemDescriptionLocalNeighborApproximation::dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p)
    {
        // reset out
        MatSetScalar(out, 0, 1, Nx_);

        // consider agent dynamics \partial f_i(x_i, u_i) / \partial x_i
        agent_->get_agentModel()->dfdx_vec(out, t, x, u, vec);

        for(const auto& neighbor : agent_->get_neighbors())
        {
            const auto j = neighbor->get_id();

            if(neighbor->is_sendingNeighbor())
            {
                // consider coupling dynamics \partial f_{ij}(x_i, u_i, x_j, u_j) / \partial x_i
                neighbor->get_couplingModel()->dfdxi_vec(out, t, x, u, x + x_index_xji_[j], u + u_index_uji_[j], vec);

                // consider coupling dynamics \partial f_{ij}(x_i, u_i, x_j, u_j) / \partial x_{ji}
                neighbor->get_couplingModel()->dfdxj_vec(out + x_index_xji_[j], t, x, u, x + x_index_xji_[j], u + u_index_uji_[j], vec);
            }

            // approximate neighbors dynamics
            neighbor->get_neighborApproximation()->dfdx_vec(out, t, x, u, vec, x_index_xji_, u_index_uji_, u_index_vji_);
        }
    }

    void ProblemDescriptionLocalNeighborApproximation::dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p)
    {
        // reset out
        MatSetScalar(out, 0, 1, Nu_);

        // consider agent dynamics \partial f_i(x_i, u_i) / \partial u_i
        agent_->get_agentModel()->dfdu_vec(out, t, x, u, vec);

        for(const auto& neighbor : agent_->get_neighbors())
        {
            if(neighbor->is_sendingNeighbor())
            {
                const auto j = neighbor->get_id();

                // consider coupling dynamics \partial f_{ij}(x_i, u_i, x_j, u_j) / \partial u_i
                neighbor->get_couplingModel()->dfdui_vec(out, t, x, u, x + x_index_xji_[j], u + u_index_uji_[j], vec);

                // consider coupling dynamics \partial f_{ij}(x_i, u_i, x_j, u_j) / \partial u_j
                neighbor->get_couplingModel()->dfduj_vec(out + u_index_uji_[j], t, x, u, x + x_index_xji_[j], u + u_index_uji_[j], vec);
            }

            // approximate neighbors dynamics
            neighbor->get_neighborApproximation()->dfdu_vec(out, t, x, u, vec, x_index_xji_, u_index_uji_, u_index_vji_);
        }
    }

    void ProblemDescriptionLocalNeighborApproximation::lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
    {
        // reset out
        MatSetScalar(out, 0, 1, 1);
        const int Nxi = agent_->get_Nxi();

        typeRNum l_i = 0.0;
        typeRNum l_ij = 0.0;

        // *************************
        // own cost
        //*************************

        interpolateState(agent_->get_desiredAgentState(), t, desired_state_);

        // l_i( x_i, u_i )
        agent_->get_agentModel()->lfct(&l_i, t, x, u, &desired_state_.x_[0]);

        for (const auto& neighbor : agent_->get_sendingNeighbors())
        {
            const auto j = neighbor->get_id();

            neighbor->get_couplingModel()->lfct(&l_ij, t, x, u, x + x_index_xji_[j], u + u_index_uji_[j]);
        }

        // *************************
        // Approx cost function
        //*************************

        if (!agent_->is_approximatingCost())
        {
            out[0] += l_i + l_ij;
        }
        else
        {
			// rescale agent cost
            out[0] += l_i / (1.0 + agent_->get_neighbors().size()) + l_ij / 2.0;

            // add the neighbors approximated cost
            for (const auto& neighbor : agent_->get_neighbors())
            {
				const auto j = neighbor->get_id(); 

                typeRNum l_j = 0.0;
				typeRNum l_ji = 0.0;

                interpolateState(neighbor->get_neighbors_desiredAgentState(), t, desired_state_);

                // l_j(x_i, u_j)
                neighbor->get_agentModel()->lfct(&l_j, t, x + x_index_xji_[j], u + u_index_uji_[j], &desired_state_.x_[0]);

                if (neighbor->is_receivingNeighbor())
                    neighbor->get_copied_couplingModel()->lfct(&l_ji, t, x + x_index_xji_[j], u + u_index_uji_[j], x, u);

                // add neighbors cost with respect to normalizing factor
                out[0] += l_j / (1.0 + neighbor->get_numberOfNeighbors()) + l_ji / 2.0;
            }
        }

        // *************************
        // Consider own constraints
        //*************************

        interpolateState(agent_->get_couplingState(), t, couplingState_);
        interpolateState(agent_->get_multiplierState(), t, multiplierState_);
        interpolateState(agent_->get_penaltyState(), t, penaltyState_);

        // consistency constraints ( zu_i - u_i )
        for (unsigned int k = 0; k < agent_->get_Nui(); ++k)
        {
            const typeRNum zu_min_u_ = couplingState_.z_u_[k] - u[k];
            out[0] += multiplierState_.mu_u_[k] * zu_min_u_ + 0.5 * penaltyState_.rho_u_[k] * zu_min_u_ * zu_min_u_;
        }

        // *************************
        // Consider neighbor constraints
        //*************************

        for(const auto& neighbor : agent_->get_neighbors())
        {
            const auto j = neighbor->get_id();

            // consistency constraint (zv_{ij} - v_{ij})
            interpolateState(neighbor->get_externalInfluence_couplingState(), t, couplingState_);
            interpolateState(neighbor->get_externalInfluence_multiplierState(), t, multiplierState_);
            interpolateState(neighbor->get_externalInfluence_penaltyState(), t, penaltyState_);

            // consider constraints (zv_i - v_{ij})
            std::vector<typeRNum> v(Nxi, 0.0);
            neighbor->get_neighborApproximation()->vfct(&v[0], t, x, u, x_index_xji_, u_index_uji_, u_index_vji_);
            for(unsigned int k = 0; k < neighbor->get_Nxj(); ++k)
            {
                const typeRNum zv_min_v_ = couplingState_.z_v_[k] - v[k];
                out[0] += multiplierState_.mu_v_[k] * zv_min_v_ + 0.5 * penaltyState_.rho_v_[k] * zv_min_v_ * zv_min_v_;
            }

            // consistency constraint (zv_{ji} - v_{ji})
            interpolateState(neighbor->get_neighbors_externalInfluence_couplingState(), t, couplingState_);
            interpolateState(neighbor->get_coupled_multiplierState(), t, multiplierState_);
            interpolateState(neighbor->get_coupled_penaltyState(), t, penaltyState_);

            for(unsigned int k = 0; k < neighbor->get_Nxj(); ++k)
            {
                const typeRNum zv_min_v_ = couplingState_.z_v_[k] - (u + u_index_vji_[j])[k];
                out[0] += multiplierState_.mu_v_[k] * zv_min_v_ + 0.5 * penaltyState_.rho_v_[k] * zv_min_v_ * zv_min_v_;
            }

            // consistency constraint (zu_j - u_{ji})
            interpolateState(neighbor->get_neighbors_couplingState(), t, couplingState_);
            for(unsigned int k = 0; k < neighbor->get_Nuj(); ++k)
            {
                const typeRNum zu_min_u_ = couplingState_.z_u_[k] - (u + u_index_uji_[j])[k];
                out[0] += multiplierState_.mu_u_[k] * zu_min_u_ + 0.5 * penaltyState_.rho_u_[k] * zu_min_u_ * zu_min_u_;
            }
        }
    }

    void ProblemDescriptionLocalNeighborApproximation::dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
    {
        // reset out
        MatSetScalar(out, 0, 1, Nx_);
        const unsigned int Nxi = agent_->get_Nxi();

		std::vector<typeRNum> l_i(Nx_, 0.0);
		std::vector<typeRNum> l_ij(Nx_, 0.0);

        // *************************
        // own cost
        //*************************

        interpolateState(agent_->get_desiredAgentState(), t, desired_state_);

        // \partial l_i(x_i, u_i) / \partial x_i
        agent_->get_agentModel()->dldx(&l_i[0], t, x, u, &desired_state_.x_[0]);

        for (const auto& neighbor : agent_->get_sendingNeighbors())
        {
            const auto j = neighbor->get_id();

            neighbor->get_couplingModel()->dldxi(&l_ij[0], t, x, u, x + x_index_xji_[j], u + u_index_uji_[j]);

			neighbor->get_couplingModel()->dldxj(&l_ij[0] + x_index_xji_[j], t, x, u, x + x_index_xji_[j], u + u_index_uji_[j]);
        }

        // *************************
        // approximate neighbors cost
        //*************************

        if (!agent_->is_approximatingCost())
        {
            for (unsigned int k = 0; k < Nx_; ++k)
                out[k] += l_i[k] + l_ij[k];
        }
        else
        {
            // rescale cost
            for (unsigned int k = 0; k < Nx_; ++k)
                out[k] += l_i[k] / (1.0 + agent_->get_neighbors().size()) + l_ij[k] / 2.0;

            // consider approximated cost
            for (const auto& neighbor : agent_->get_neighbors())
            {
                const auto j = neighbor->get_id();
				std::vector<typeRNum> l_j(Nx_, 0.0);
				std::vector<typeRNum> l_ji(Nx_, 0.0);

                interpolateState(neighbor->get_neighbors_desiredAgentState(), t, desired_state_);

                // \partial l_j(x_j, u_i) / \partial d_x_j
                neighbor->get_agentModel()->dldx(&l_j[0] + x_index_xji_[j], t, x + x_index_xji_[j], u + u_index_uji_[j], &desired_state_.x_[0]);

                if (neighbor->is_receivingNeighbor())
                {
                    neighbor->get_copied_couplingModel()->dldxj(&l_ji[0], t, x + x_index_xji_[j], u + u_index_uji_[j], x, u);

					neighbor->get_copied_couplingModel()->dldxi(&l_ji[0] + x_index_xji_[j], t, x + x_index_xji_[j], u + u_index_uji_[j], x, u);
                }

                for (unsigned int k = 0; k < Nx_; ++k)
                    out[k] += l_j[k] / (1.0 + neighbor->get_numberOfNeighbors()) + l_ji[k] / 2.0;
            }
        }

        // *************************
        // constraints
        //*************************
        std::vector<typeRNum> vec(0, 0.0);
        for(const auto& neighbor : agent_->get_neighbors())
        {
            interpolateState(neighbor->get_externalInfluence_couplingState(), t, couplingState_);
            interpolateState(neighbor->get_externalInfluence_multiplierState(), t, multiplierState_);
            interpolateState(neighbor->get_externalInfluence_penaltyState(), t, penaltyState_);

            // evaluate v_{ij}( x_i, u_i, x_{ji}, u_{ji} )
            std::vector<typeRNum> vij( Nxi, 0.0 );
            neighbor->get_neighborApproximation()->vfct( &vij[0], t, x, u, x_index_xji_, u_index_uji_, u_index_vji_ );

            // d \mu_{ij}*( zv_{ij} - v_{ij}( x_i, u_i, x_{ji}, u_{ji} ) ) / dx
            vec = multiplierState_.mu_v_;
            for(unsigned int k = 0; k < Nxi; ++k )
                vec[k] = vec[k] * (-1);

            neighbor->get_neighborApproximation()->dvdx_vec( out, t, x, u, &vec[0], x_index_xji_, u_index_uji_, u_index_vji_ );

            // d 0.5*rho*( zv_{ij} - v_{ij}( x_i, u_i, x_{ji}, u_{ji} ) )^2 / dx
           for(unsigned int k = 0; k < Nxi; ++k )
               vec[k] = penaltyState_.rho_v_[k] * (couplingState_.z_v_[k] - vij[k]) * (-1);

           neighbor->get_neighborApproximation()->dvdx_vec( out, t, x, u, &vec[0], x_index_xji_, u_index_uji_, u_index_vji_ );
        }
    }

    void ProblemDescriptionLocalNeighborApproximation::dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
    {
        // reset out
        MatSetScalar(out, 0, 1, Nu_);
        const unsigned int Nui = agent_->get_Nui();
		const unsigned int Nxi = agent_->get_Nxi();

		std::vector<typeRNum> l_i(Nu_, 0.0);
		std::vector<typeRNum> l_ij(Nu_, 0.0);

        // *************************
        // own cost
        //*************************

        // \partial l_i( x_i, u_i ) / \partial u_i
        interpolateState(agent_->get_desiredAgentState(), t, desired_state_);
        agent_->get_agentModel()->dldu(&l_i[0], t, x, u, &desired_state_.x_[0]);

        for (const auto& neighbor : agent_->get_sendingNeighbors())
        {
            const auto j = neighbor->get_id();

            neighbor->get_couplingModel()->dldui(&l_ij[0], t, x, u, x + x_index_xji_[j], u + u_index_uji_[j]);

			neighbor->get_couplingModel()->dlduj(&l_ij[0] + u_index_uji_[j], t, x, u, x + x_index_xji_[j], u + u_index_uji_[j]);
        }

        // *************************
        // approximate neighbors cost
        //*************************

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
            for (const auto& neighbor : agent_->get_neighbors())
            {
                interpolateState(neighbor->get_neighbors_desiredAgentState(), t, desired_state_);

                const auto j = neighbor->get_id();
				std::vector<typeRNum> l_j(Nu_, 0.0);
				std::vector<typeRNum> l_ji(Nu_, 0.0);

                // consider local copies u_{ji} as control
                neighbor->get_agentModel()->dldu(&l_j[0] + u_index_uji_[j], t, x + x_index_xji_[j], u + u_index_uji_[j], &desired_state_.x_[0]);

                if (neighbor->is_receivingNeighbor())
                {
                    neighbor->get_copied_couplingModel()->dlduj(&l_ji[0], t, x + x_index_xji_[j], u + u_index_uji_[j], x, u);

					neighbor->get_copied_couplingModel()->dldui(&l_ji[0] + u_index_uji_[j], t, x + x_index_xji_[j], u + u_index_uji_[j], x, u);
                }

                for (unsigned int k = 0; k < Nu_; ++k)
                    out[k] += l_j[k] / (1.0 + neighbor->get_numberOfNeighbors()) + l_ji[k] / 2.0;
            }
        }

        // *************************
        // own constraints
        //*************************


        interpolateState(agent_->get_couplingState(), t, couplingState_);
        interpolateState(agent_->get_multiplierState(), t, multiplierState_);
        interpolateState(agent_->get_penaltyState(), t, penaltyState_);

        // consistency constraints ( zu_i - u_i )
        // derivative w.r.t. u_i
        for (unsigned int k = 0; k < agent_->get_Nui(); ++k)
        {
            const typeRNum zu_min_u_ = couplingState_.z_u_[k] - u[k];
            out[k] += (multiplierState_.mu_u_[k] + penaltyState_.rho_u_[k] * zu_min_u_) * (-1);
        }

        // *************************
        // neighbors constraints
        //*************************

        // consistency constraints for neighbors
        for(const auto& neighbor : agent_->get_neighbors())
        {
            const auto j = neighbor->get_id();

            // consistency constraint ( zv_{ji} - v_{ji} )
            // derivative w.r.t. v_{ji}
            interpolateState(neighbor->get_neighbors_externalInfluence_couplingState(), t, couplingState_);
            interpolateState(neighbor->get_coupled_multiplierState(), t, multiplierState_);
            interpolateState(neighbor->get_coupled_penaltyState(), t, penaltyState_);

            for(unsigned int k = 0; k < couplingState_.z_v_.size(); ++k)
            {
                const typeRNum zv_min_v_ = couplingState_.z_v_[k] - (u + u_index_vji_[j])[k];
                (out + u_index_vji_[j])[k] += (multiplierState_.mu_v_[k] + penaltyState_.rho_v_[k] * zv_min_v_) * (-1);
            }

            // consistency constraint ( zu_j - u_{ji} )
            // derivative w.r.t. u_{ji}
            interpolateState(neighbor->get_neighbors_couplingState(), t, couplingState_);

            for(unsigned int k = 0; k < couplingState_.z_u_.size(); ++k)
            {
                const typeRNum zu_min_u_ = couplingState_.z_u_[k] - (u + u_index_uji_[j])[k];
                (out + u_index_uji_[j])[k] += (multiplierState_.mu_u_[k] + penaltyState_.rho_u_[k] * zu_min_u_) * (-1);
            }

            // consistency constraint ( zv_{ij} - v_{ij}(x_i, u_i, x_{ji}, u_{ji}) )
            interpolateState(neighbor->get_externalInfluence_couplingState(), t, couplingState_);
            interpolateState(neighbor->get_externalInfluence_multiplierState(), t, multiplierState_);
            interpolateState(neighbor->get_externalInfluence_penaltyState(), t, penaltyState_);

            // evaluate v_{ij}
            std::vector<typeRNum> v( Nxi, 0.0 );
            neighbor->get_neighborApproximation()->vfct( &v[0], t, x, u, x_index_xji_, u_index_uji_, u_index_vji_ );

            // d \mu_{ij}*( zv_{ij} - v_{ij}(x_i, u_i, x_{ji}, u_{ji}) ) / du
            std::vector<typeRNum> vec = multiplierState_.mu_v_;
            for (unsigned int k = 0; k < Nxi; ++k)
                vec[k] = vec[k] * (-1);

            neighbor->get_neighborApproximation()->dvdu_vec( out, t, x, u, &vec[0], x_index_xji_, u_index_uji_, u_index_vji_ );

            // d 0.5*rho*( zv_{ij} - v_{ij}(x_i, u_i, x_{ji}, u_{ji}) ) / du
            for (unsigned int k = 0; k < Nxi; ++k)
                vec[k] = penaltyState_.rho_v_[k] * (couplingState_.z_v_[k] - v[k]) * (-1);

            neighbor->get_neighborApproximation()->dvdu_vec( out, t, x, u, &vec[0], x_index_xji_, u_index_uji_, u_index_vji_ );
        }
    }

    void ProblemDescriptionLocalNeighborApproximation::Vfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes)
    {
        MatSetScalar(out, 0, 1, 1);
        interpolateState(agent_->get_desiredAgentState(), t, desired_state_);

        agent_->get_agentModel()->Vfct(out, t, x, &desired_state_.x_[0]);

        // consider approximated cost
        if(agent_->is_approximatingCost())
        {
            // rescale cost
            out[0] /= 1.0 + agent_->get_neighbors().size();

            for( const auto& neighbor : agent_->get_neighbors() )
            {
                const auto j = neighbor->get_id();
                typeRNum lj = 0.0;
                interpolateState(neighbor->get_desiredAgentState(), t, desired_state_);

                neighbor->get_agentModel()->Vfct(&lj, t, x + x_index_xji_[j], &desired_state_.x_[0]);
                out[0] += lj / ( 1 + neighbor->get_numberOfNeighbors() );
            }
        }
    }

    void ProblemDescriptionLocalNeighborApproximation::dVdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes)
    {
        MatSetScalar(out, 0, 1, Nx_);
        interpolateState(agent_->get_desiredAgentState(), t, desired_state_);

        agent_->get_agentModel()->dVdx(out, t, x, &desired_state_.x_[0]);

        if(agent_->is_approximatingCost())
        {
            // rescale cost
            for(unsigned int k = 0; k < agent_->get_Nxi(); ++k)
                out[k] /= 1 + agent_->get_neighbors().size();

            for(const auto& neighbor : agent_->get_neighbors())
            {
                const auto j = neighbor->get_id();
                const auto Nxj_ = neighbor->get_Nxi();

                std::vector<typeRNum> l(Nxj_, 0.0);
                interpolateState(neighbor->get_desiredAgentState(), t, desired_state_);

                neighbor->get_agentModel()->dVdx(&l[0], t, x + x_index_xji_[j], &desired_state_.x_[0]);
                for(unsigned int k = 0; k < Nxj_; ++k)
                    (out + x_index_xji_[j])[k] += l[k] / (1 + neighbor->get_numberOfNeighbors());
            }
        }
    }

    void ProblemDescriptionLocalNeighborApproximation::gfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
    {
        MatSetScalar(out, 0, 1, Ng_);

        // equality constraints g_i(x_i, u_i) = 0
        agent_->get_agentModel()->gfct(out, t, x, u);
        auto idx_ = agent_->get_agentModel()->get_Ngi();

        for(const auto& neighbor : agent_->get_neighbors())
        {
            const auto j = neighbor->get_id();

            // equality constraints g_{ij}(x_i, u_i, x_j, u_j) = 0
           if( neighbor->is_sendingNeighbor() )
           {
               neighbor->get_couplingModel()->gfct(out + idx_, t, x, u, x + x_index_xji_[j], u + u_index_uji_[j]);
               idx_ += neighbor->get_couplingModel()->get_Ngij();
           }

           if( agent_->is_approximatingConstraints() )
           {
               // equality constraints g_j(x_{ji}, u_{ji})
               neighbor->get_agentModel()->gfct(out + idx_, t, x + x_index_xji_[j], u + u_index_uji_[j]);
               idx_ += neighbor->get_agentModel()->get_Ngi();

               if( neighbor->is_receivingNeighbor() )
               {
                   // equality constraints g_{ji}(x_{ji}, u_{ji}, x_i, u_i)
                   neighbor->get_copied_couplingModel()->gfct(out + idx_, t, x + x_index_xji_[j], u + u_index_uji_[j], x, u);
                   idx_ += neighbor->get_copied_couplingModel()->get_Ngij();
               }
           }
        }
    }

    void ProblemDescriptionLocalNeighborApproximation::dgdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
    {
        MatSetScalar(out, 0, 1, Nx_);

        // equality constraints \partial g_i(x_i, u_i) / \partial x_i
        agent_->get_agentModel()->dgdx_vec(out, t, x, u, vec);
        auto idx_ = agent_->get_agentModel()->get_Ngi();

        for(const auto& neighbor : agent_->get_neighbors())
        {
            const auto j = neighbor->get_id();

            if(neighbor->is_sendingNeighbor())
            {
                // equality constraints \partial g_{ij}(x_i, u_i, x_j, u_j) / \partial x_i
                neighbor->get_couplingModel()->dgdxi_vec(out, t, x, u, x + x_index_xji_[j], u + u_index_uji_[j], vec + idx_);

                // equality constraints \partial g_{ij}(x_i, u_i, x_j, u_j) / \partial x_{ji}
                neighbor->get_couplingModel()->dgdxj_vec(out + x_index_xji_[j], t, x, u, x + x_index_xji_[j], u + u_index_uji_[j], vec + idx_);

                // index is increased only once as it's the same constraint twice in this loop
                idx_ += neighbor->get_couplingModel()->get_Ngij();
            }

            if(agent_->is_approximatingConstraints())
            {
                // equality constraints \partial g_j(x_{ji}, u_{ji}) / \partial x_{ji}
                neighbor->get_agentModel()->dgdx_vec(out + x_index_xji_[j], t, x + x_index_xji_[j], u + u_index_uji_[j], vec + idx_);
                idx_ += neighbor->get_agentModel()->get_Ngi();

                if(neighbor->is_receivingNeighbor())
                {
                    // equality constraints \partial g_{ji}(x_{ji}, u_{ji}, x_i, u_i) / \partial x_i
                    neighbor->get_copied_couplingModel()->dgdxj_vec(out, t, x + x_index_xji_[j], u + u_index_uji_[j], x, u, vec + idx_);

                    // equality constraints \partial g_{ji}(x_{ji}, u_{ji}, x_i, u_i) / \partial x_{ji}
                    neighbor->get_copied_couplingModel()->dgdxi_vec(out + x_index_xji_[j], t, x + x_index_xji_[j], u + u_index_uji_[j], x, u, vec + idx_);

                    // index is increased only once as it's the same constraint twice in this loop
                    idx_ += neighbor->get_copied_couplingModel()->get_Ngij();
                }
            }
        }
    }

    void ProblemDescriptionLocalNeighborApproximation::dgdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
    {
        MatSetScalar(out, 0, 1, Nu_);

        // equality constraints \partial g_i(x_i, u_i) / \partial u_i
        agent_->get_agentModel()->dgdu_vec(out, t, x, u, vec);
        auto idx_ = agent_->get_agentModel()->get_Ngi();

        for(const auto& neighbor : agent_->get_neighbors())
        {
            const auto j = neighbor->get_id();

            if( neighbor->is_sendingNeighbor() )
            {
                // equality constraints \partial g_{ij}(x_i, u_i, x_j, u_j) / \partial u_i
                neighbor->get_couplingModel()->dgdui_vec(out, t, x, u, x + x_index_xji_[j], u + u_index_uji_[j], vec + idx_);

                // equality constraints \partial g_{ij}(x_i, u_i, x_j, u_j) / \partial u_j
                neighbor->get_couplingModel()->dgduj_vec(out + u_index_uji_[j], t, x, u, x + x_index_xji_[j], u + u_index_uji_[j], vec + idx_);

                // index is increased only once as it's the same constraint twice in this loop
                idx_ += neighbor->get_couplingModel()->get_Ngij();
            }
            if( agent_->is_approximatingConstraints() )
            {
                // equality constraints \partial g_j(x_{ji}, u_{ji}) / \partial u_{ji}
                neighbor->get_agentModel()->dgdu_vec(out + u_index_uji_[j], t, x + x_index_xji_[j], u + u_index_uji_[j], vec + idx_);

                // increase index
                idx_ += neighbor->get_agentModel()->get_Ngi();

                if(neighbor->is_receivingNeighbor())
                {
                    // equality constraints \partial g_{ji}(x_{ji}, u_{ji}, x_i, u_i) / \partial u_i
                    neighbor->get_copied_couplingModel()->dgduj_vec(out, t, x + x_index_xji_[j], u + u_index_uji_[j], x, u, vec + idx_);

                    // equality constraints \partial g_{ji}(x_{ji}, u_{ji}, x_i, u_i) / \partial u_{ji}
                    neighbor->get_copied_couplingModel()->dgdui_vec(out + u_index_uji_[j], t, x + x_index_xji_[j], u + u_index_uji_[j], x, u, vec + idx_);

                    // index is increased only once as it's the same constraint twice in this loop
                    idx_ += neighbor->get_copied_couplingModel()->get_Ngij();
                }
            }
        }
    }

    void ProblemDescriptionLocalNeighborApproximation::hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
    {
        MatSetScalar(out, 0, 1, Nh_);

        // inequality constraints h_i(x_i, u_i) <= 0
        agent_->get_agentModel()->hfct(out, t, x, u);
        auto idx_ = agent_->get_agentModel()->get_Nhi();

        for(const auto& neighbor : agent_->get_neighbors())
        {
            const auto j = neighbor->get_id();

            if(neighbor->is_sendingNeighbor())
		    {
			    // inequality constraints h_{ij}(x_i, u_i, x_j, u_j) <= 0
                neighbor->get_couplingModel()->hfct(out + idx_, t, x, u, x + x_index_xji_[j], u + u_index_uji_[j]);
                idx_ += neighbor->get_couplingModel()->get_Nhij();
            }

            if(agent_->is_approximatingConstraints())
            {
                // inequality constraint h_j(x_{ji}, u_{ji})
                neighbor->get_agentModel()->hfct(out + idx_, t, x + x_index_xji_[j], u + u_index_uji_[j]);
                idx_ += neighbor->get_agentModel()->get_Nhi();

                if(neighbor->is_receivingNeighbor())
                {
                    // inequality constraint h_{ji}(x_{ji}, u_{ji}, x_i, u_i)
                    neighbor->get_copied_couplingModel()->hfct(out + idx_, t, x + x_index_xji_[j], u + u_index_uji_[j], x, u);
                    idx_ += neighbor->get_copied_couplingModel()->get_Nhij();
                }
            }
        }
    }

    void ProblemDescriptionLocalNeighborApproximation::dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
    {
        MatSetScalar(out, 0, 1, Nx_);

        // inequality constraints \partial h_i(x_i, u_i) / \partial x_i
        agent_->get_agentModel()->dhdx_vec(out, t, x, u, vec);
        auto idx_ = agent_->get_agentModel()->get_Nhi();

        for(const auto& neighbor : agent_->get_neighbors())
        {
            const auto j = neighbor->get_id();

            if(neighbor->is_sendingNeighbor())
		    {
			    // inequality constraints \partial h_{ij}(x_i, u_i, x_j, u_j) / \partial x_i
                neighbor->get_couplingModel()->dhdxi_vec(out, t, x, u, x + x_index_xji_[j], u + u_index_uji_[j], vec + idx_);

                // inequality constraints \partial h_{ij}(x_i, u_i, x_j, u_j) / \partial x_{ji}
                neighbor->get_couplingModel()->dhdxj_vec(out + x_index_xji_[j], t, x, u, x + x_index_xji_[j], u + u_index_uji_[j], vec + idx_);

                // index is increased only once as it's the same constraint twice in this loop
                idx_ += neighbor->get_couplingModel()->get_Nhij();
            }

            if(agent_->is_approximatingConstraints())
            {
                // inequality constraint \partial h_j(x_{ji}, u_{ji}) / \partial x_{ji}
                neighbor->get_agentModel()->dhdx_vec(out + x_index_xji_[j], t, x + x_index_xji_[j], u + u_index_uji_[j], vec + idx_);
                idx_ += neighbor->get_agentModel()->get_Nhi();

                if(neighbor->is_receivingNeighbor())
                {
                    // inequality constraint \partial h_{ji}(x_{ji}, u_{ji}, x_i, u_i) / \partial x_i
                    neighbor->get_copied_couplingModel()->dhdxj_vec(out, t, x + x_index_xji_[j], u + u_index_uji_[j], x, u, vec + idx_);

                    // inequality constraint \partial h_{ji}(x_{ji}, u_{ji}, x_i, u_i) / \partial x_{ji}
                    neighbor->get_copied_couplingModel()->dhdxi_vec(out + x_index_xji_[j], t, x + x_index_xji_[j], u + u_index_uji_[j], x, u, vec + idx_);

                    // index is increased only once as it's the same constraint twice in this loop
                    idx_ += neighbor->get_copied_couplingModel()->get_Nhij();
                }
            }
        }
    }

    void ProblemDescriptionLocalNeighborApproximation::dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
    {
	    MatSetScalar(out, 0, 1, Nu_);

        // inequality constraints \partial h_i(x_i, u_i) / \partial u_i
        agent_->get_agentModel()->dhdu_vec(out, t, x, u, vec);
        auto idx_ = agent_->get_agentModel()->get_Nhi();

        for(const auto& neighbor : agent_->get_neighbors())
        {
            const auto j = neighbor->get_id();

            if(neighbor->is_sendingNeighbor())
            {
                // inequality constraints \partial h_{ij}(x_i, u_i, x_j, u_j) / \partial u_i
                neighbor->get_couplingModel()->dhdui_vec(out, t, x, u, x + x_index_xji_[j], u + u_index_uji_[j], vec + idx_);

                // inequality constraints \partial h_{ij}(x_i, u_i, x_j, u_j) / \partial u_j
                neighbor->get_couplingModel()->dhduj_vec(out + u_index_uji_[j], t, x, u, x + x_index_xji_[j], u + u_index_uji_[j], vec + idx_);

                // index is increased only once as it's the same constraint twice in this loop
                idx_ += neighbor->get_couplingModel()->get_Nhij();
            }

            if(agent_->is_approximatingConstraints())
            {
                // inequality constraint \partial h_j(x_{ji}, u_{ji}) / \partial u_{ji}
                neighbor->get_agentModel()->dhdu_vec(out + u_index_uji_[j], t, x + x_index_xji_[j], u + u_index_uji_[j], vec + idx_);
                idx_ += neighbor->get_agentModel()->get_Nhi();

                if(neighbor->is_receivingNeighbor())
                {
                    // inequality constraint \partial h_{ji}(x_{ji}, u_{ji}, x_i, u_i) / \partial u_i
                    neighbor->get_copied_couplingModel()->dhduj_vec(out, t, x + x_index_xji_[j], u + u_index_uji_[j], x, u, vec + idx_);

                    // inequality constraint \partial h_{ji}(x_{ji}, u_{ji}, x_i, u_i) / \partial u_{ji}
                    neighbor->get_copied_couplingModel()->dhdui_vec(out + u_index_uji_[j], t, x + x_index_xji_[j], u + u_index_uji_[j], x, u, vec + idx_);

                    // index is increased only once as it's the same constraint twice in this loop
                    idx_ += neighbor->get_copied_couplingModel()->get_Nhij();
                }
            }
        }
    }

}
