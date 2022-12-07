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

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(i_, t0_, t_, mu_x_, mu_u_, mu_v_);
		}
	};

}
