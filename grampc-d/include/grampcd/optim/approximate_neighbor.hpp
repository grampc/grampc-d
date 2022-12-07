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

#pragma once

#include "grampcd/util/class_forwarding.hpp"

#include "grampcd/state/agent_state.hpp"

namespace grampcd 
{

    /* @brief Approximation of neighbors dynamics. */
    class ApproximateNeighbor
    {
    public:
        ApproximateNeighbor(Agent *agent, Neighbor *neighbor);

        /*Dynamics*/
        void ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u,
                          const std::vector<int>& x_index_xji_ptr, const std::vector<int>& u_index_uji_ptr, const std::vector<int>& u_index_vji_ptr );

        /*Partial derivate of the dynamics with respect to state multiplied with Lagrangian multipliers*/
        void dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *vec,
                              const std::vector<int>& x_index_xji_ptr, const std::vector<int>& u_index_uji_ptr, const std::vector<int>& u_index_vji_ptr);

        /*Partial derivate of the dynamics with respect to controls multiplied with Lagrangian multipliers*/
        void dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *vec,
                              const std::vector<int>& x_index_xji_ptr, const std::vector<int>& u_index_uji_ptr, const std::vector<int>& u_index_vji_ptr);

        /*External influence*/
        void vfct(typeRNum *out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u,
                          const std::vector<int>& x_index_xji_ptr, const std::vector<int>& u_index_uji_ptr, const std::vector<int>& u_index_vji_ptr);

        /*Partial derivate of the external influence with respect to states multiplied with Lagrangian multipliers*/
        void dvdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *vec,
                              const std::vector<int>& x_index_xji_ptr, const std::vector<int>& u_index_uji_ptr, const std::vector<int>& u_index_vji_ptr);

	    /*Partial derivate of the external influence with respect to controls multiplied with Lagrangian multipliers*/
        void dvdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *vec,
                              const std::vector<int>& x_index_xji_ptr, const std::vector<int>& u_index_uji_ptr, const std::vector<int>& u_index_vji_ptr);

    protected:
        Agent *agent_;
        Neighbor *neighbor_;
        AgentState desired_state_;
    };

}