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

#include "grampcd/optim/approximate_neighbor.hpp"

#include "grampcd/agent/agent.hpp"
#include "grampcd/agent/neighbor.hpp"

#include "grampcd/model/agent_model.hpp"
#include "grampcd/model/coupling_model.hpp"

namespace grampcd
{

    ApproximateNeighbor::ApproximateNeighbor(Agent *agent, Neighbor *neighbor)
        : agent_( agent),
          neighbor_(neighbor)
    {
    }

    void ApproximateNeighbor::ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u,
                                            const std::vector<int>& x_index_xji_ptr, const std::vector<int>& u_index_uji_ptr, const std::vector<int>& u_index_vji_ptr)
    {
        // approximate dynamics \dot x_{ji} = v_{ji} + f_j(x_{ji}, u_{ji}) + f_{ji}(x_{ji}, u_{ji}, x_i, u_i)
        const auto x_index_xji_ = x_index_xji_ptr[neighbor_->get_id()];
        const auto u_index_vji_ = u_index_vji_ptr[neighbor_->get_id()];
        const auto u_index_uji_ = u_index_uji_ptr[neighbor_->get_id()];

        // part v_{ji}
        MatCopy(out + x_index_xji_, u + u_index_vji_, 1, neighbor_->get_Nxj());

        // part f_j( x_{ji}, u_{ji} )
        neighbor_->get_agentModel()->ffct(out + x_index_xji_, t, x + x_index_xji_, u + u_index_uji_);

        // part f_{ji}( x_{ji}, u_{ji}, x_i, u_i )
        if( neighbor_->is_receivingNeighbor() )
            neighbor_->get_copied_couplingModel()->ffct(out + x_index_xji_, t, x + x_index_xji_, u + u_index_uji_, x, u);
    }

    void ApproximateNeighbor::dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *vec,
                                                const std::vector<int>& x_index_xji_ptr, const std::vector<int>& u_index_uji_ptr, const std::vector<int>& u_index_vji_ptr)
    {
        const auto x_index_xji_ = x_index_xji_ptr[neighbor_->get_id()];
        const auto u_index_uji_ = u_index_uji_ptr[neighbor_->get_id()];

        // \partial f_j(x_{ji}, u_{ji}) / \partial x_{ji}
        neighbor_->get_agentModel()->dfdx_vec(out + x_index_xji_, t, x + x_index_xji_, u + u_index_uji_, vec + x_index_xji_);

        if( neighbor_->is_receivingNeighbor() )
        {
            // \partial f_{ji}(x_{ji}, u_{ji}, x_i, u_i) / \partial x_i
            neighbor_->get_copied_couplingModel()->dfdxj_vec(out, t, x + x_index_xji_, u + u_index_uji_, x, u, vec + x_index_xji_);

            // \partial f_{ji}(x_{ji}, u_{ji}, x_i, u_i) / \partial x_{ji}
            neighbor_->get_copied_couplingModel()->dfdxi_vec(out + x_index_xji_, t, x + x_index_xji_, u + u_index_uji_, x, u, vec + x_index_xji_);
        }
    }

    void ApproximateNeighbor::dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *vec,
                                                const std::vector<int>& x_index_xji_ptr, const std::vector<int>& u_index_uji_ptr, const std::vector<int>& u_index_vji_ptr)
    {
        // approximate dynamics \dot x_{ji} = v_{ji} + f_j(x_{ji}, u_{ji}) + f_{ji}(x_{ji}, u_{ji}, x_i, u_i)
        const auto x_index_xji_ = x_index_xji_ptr[neighbor_->get_id()];
        const auto u_index_vji_ = u_index_vji_ptr[neighbor_->get_id()];
        const auto u_index_uji_ = u_index_uji_ptr[neighbor_->get_id()];

        // \partial v_{ji} / \partial v_{ji}
        MatCopy(out + u_index_vji_, vec + x_index_xji_, 1, neighbor_->get_Nxj());

        // \partial f_j(x_{ji}, u_{ji}) / \partial u_{ji}
        neighbor_->get_agentModel()->dfdu_vec(out + u_index_uji_, t, x + x_index_xji_, u + u_index_uji_, vec + x_index_xji_);

        if( neighbor_->is_receivingNeighbor() )
        {
            // \partial f_{ji}(x_{ji}, u_{ji}, x_i, u_i) / \partial u_i
            neighbor_->get_copied_couplingModel()->dfduj_vec(out, t, x + x_index_xji_, u + u_index_uji_, x, u, vec + x_index_xji_);

            // \partial f_{ji}(x_{ji}, u_{ji}, x_i, u_i) / \partial u_{ji}
            neighbor_->get_copied_couplingModel()->dfdui_vec(out + u_index_uji_, t, x + x_index_xji_, u + u_index_uji_, x, u, vec + x_index_xji_);
        }
    }

    void ApproximateNeighbor::vfct(typeRNum *out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u,
                                            const std::vector<int>& x_index_xji_ptr, const std::vector<int>& u_index_uji_ptr, const std::vector<int>& u_index_vji_ptr)
    {
        // auxiliary variable v_{ij} = \sum_{s \in N_i \ j} f_{is}(x_i, u_i, x_{si}, u_{si})
        const auto j_ = neighbor_->get_id();

        // f_{is}(x_i, u_i, x_{si}, u_{si})
        for(const auto& neighbor : agent_->get_sendingNeighbors())
        {
            const auto s_ = neighbor->get_id();
            if(s_ != j_)
			    neighbor->get_couplingModel()->ffct(out, t, x, u, x + x_index_xji_ptr[s_], u + u_index_uji_ptr[s_]);
        }
    }

    void ApproximateNeighbor::dvdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *vec,
                                                const std::vector<int>& x_index_xji_ptr, const std::vector<int>& u_index_uji_ptr, const std::vector<int>& u_index_vji_ptr)
    {
        // auxiliary variable v_{ij} = \sum_{s \in N \ j} f_{is}(x_i, u_i, x_{si}, u_{si})
        const auto j_ = neighbor_->get_id();

        for(const auto& neighbor : agent_->get_sendingNeighbors())
        {
            const auto s_ = neighbor->get_id();
            if(s_ != j_)
            {
                // dynamic coupling \partial f_{is}(x_i, u_i, x_{si}, u_{si}) / \partial x_i
                neighbor->get_couplingModel()->dfdxi_vec(out, t, x, u, x + x_index_xji_ptr[s_], u + u_index_uji_ptr[s_], vec);

                // dynamic coupling \partial f_{is}(x_i, u_i, x_{si}, u_{si}) / \partial x_{is}
                neighbor->get_couplingModel()->dfdxj_vec(out + x_index_xji_ptr[s_], t, x, u, x + x_index_xji_ptr[s_], u + u_index_uji_ptr[s_], vec);
            }
        }
    }

    void ApproximateNeighbor::dvdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *vec,
                                                const std::vector<int>& x_index_xji_ptr, const std::vector<int>& u_index_uji_ptr, const std::vector<int>& u_index_vji_ptr)
    {
        // auxiliary variable v_{ij} = \sum_{s in N \ j} f_{is}(x_i, u_i, x_{si}, u_{si})
        const auto j_ = neighbor_->get_id();

        for(const auto& neighbor : agent_->get_sendingNeighbors())
        {
            const auto s_ = neighbor->get_id();
            if(s_ != j_)
            {
                // dynamic coupling \partial f_{is}(x_i, u_i, x_{si}, u_{si}) / \partial u_i
                neighbor->get_couplingModel()->dfdui_vec(out, t, x, u, x + x_index_xji_ptr[s_], u + u_index_uji_ptr[s_], vec);

                // dynamic coupling \partial f_{is}(x_i, u_i, x_{si}, u_{si}) / \partial u_{is}
                neighbor->get_couplingModel()->dfduj_vec(out + u_index_uji_ptr[s_], t, x, u, x + x_index_xji_ptr[s_], u + u_index_uji_ptr[s_], vec);
            }
        }
    }

}
