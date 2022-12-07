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

    struct SensiState
    {
    public:

        int i_;  // id of the agent j who uses dJ_j / du_i

        /*Initial time step*/
        typeRNum t0_;
        /*Time vector*/
        std::vector<typeRNum> t_;

        /*Vector that contains the sensitivity regarding u (dH_j/du_i)*/
        std::vector<typeRNum> psi_u_;
        /*Vector that contains the sensitivity regarding x (dH_j/dx_i)*/
        std::vector<typeRNum> psi_x_;
        /*Vector that contains the sensitivity of the terminal cost (dVji/dx_i)*/
        std::vector<typeRNum> psi_V_;

        /*Vector that contains the second-order sensitivity regarding u (d^2H_j/du_i^2)*/
        std::vector<typeRNum> psi_uu_;
        /*Vector that contains the second-order sensitivity regarding x (d^2H_j/dx_i^2)*/
        std::vector<typeRNum> psi_xx_;
        /*Vector that contains the mixed second-order sensitivity regarding xu (d^2H_j/dx_i^2)*/
        std::vector<typeRNum> psi_xu_;
        /*Vector that contains the second-order sensitivity regarding the terminal cost (d^2Vji/dx_i)*/
        std::vector<typeRNum> psi_VV_;


		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(i_, t0_, t_, psi_x_, psi_u_,psi_V_, psi_uu_, psi_xx_, psi_xu_,psi_VV_);
		}
    };
}