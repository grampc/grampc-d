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

#ifndef MULTIPLIER_STATE_HPP
#define MULTIPLIER_STATE_HPP

#include "dmpc/util/types.hpp"

namespace dmpc
{

/**
 * @brief Lagrangian multiplier mu = [mu_x, mu_u]
 */
struct MultiplierState
{
public:

    int i_;

    typeRNum t0_;
    std::vector<typeRNum> t_;

    std::vector<typeRNum> mu_x_;
    std::vector<typeRNum> mu_u_;
    std::vector<typeRNum> mu_v_;
};

typedef std::shared_ptr<MultiplierState> MultiplierStatePtr;

}

#endif // MULTIPLIER_STATE_HPP
