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
	 * @brief Coupling variable $ z = [z_x, z_u] $.
	 */
	struct CouplingState
	{
	public:

		int i_;

		typeRNum t0_;
		std::vector<typeRNum> t_;

		std::vector<typeRNum> z_x_;
		std::vector<typeRNum> z_u_;
		std::vector<typeRNum> z_v_;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(i_, t0_, t_, z_x_, z_u_, z_v_);
		}
	};

}