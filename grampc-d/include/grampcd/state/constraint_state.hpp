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

#include "grampcd/util/types.hpp"

#include "cereal/archives/binary.hpp"
#include "cereal/types/vector.hpp"

namespace grampcd
{

    struct ConstraintState
    {
    public:

        /*Agent id*/
        int i_;

        /*Initial time step*/
        typeRNum t0_;

        /*Time vector*/
        std::vector<typeRNum> t_;

        /*Vector that contains the Lagrange Multipliers of the equality constraints*/
        std::vector<typeRNum> mu_g_;
        /*Vector that contains the Lagrange Multipliers of the inequality constraints*/
        std::vector<typeRNum> mu_h_;
        /*Vector that contains the penalties of the equality constraints*/
        std::vector<typeRNum> c_g_;
        /*Vector that contains the penalties of the inequality constraints*/
        std::vector<typeRNum> c_h_;

        template<class Archive>
        void serialize(Archive& ar)
        {
            ar(i_, t0_, t_, mu_g_, mu_h_, c_g_, c_h_);
        }
    };
}