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
#pragma once

#include "grampcd/util/types.hpp"

#include "cereal/archives/binary.hpp"
#include "cereal/types/vector.hpp"

namespace grampcd
{

    struct AgentState
    {
    public:

        /*Agent id*/
        int i_;

        /*Initial time step*/
        typeRNum t0_;
        /*Time vector*/
        std::vector<typeRNum> t_;

        /*Vector that contains states*/
        std::vector<typeRNum> x_;
        /*Vector that contains controls*/
        std::vector<typeRNum> u_;
        /*Vector that contains external influence*/
		std::vector<typeRNum> v_;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(i_, t0_, t_, x_, u_, v_);
		}
    };
}