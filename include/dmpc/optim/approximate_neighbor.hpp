/* This file is part of GRAMPC-D - (https://github.com/DanielBurk/GRAMPC-D.git)
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

#ifndef APPROXIMATE_NEIGHBOR_HPP
#define APPROXIMATE_NEIGHBOR_HPP

#include "dmpc/optim/problem_description_local_neighbor_approximation.hpp"

namespace dmpc {

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

typedef std::shared_ptr<ApproximateNeighbor> ApproximateNeighborPtr;

}

#endif // APPROXIMATE_NEIGHBOR_HPP
