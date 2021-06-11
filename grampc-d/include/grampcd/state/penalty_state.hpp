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

namespace grampcd
{

	/**
	 * @brief Lagrangian multiplier $ mu = [mu_x, mu_u] $.
	 */
	struct PenaltyState
	{
	public:

		int i_;

		typeRNum t0_;
		std::vector<typeRNum> t_;

		std::vector<typeRNum> rho_x_;
		std::vector<typeRNum> rho_u_;
		std::vector<typeRNum> rho_v_;
	};

}