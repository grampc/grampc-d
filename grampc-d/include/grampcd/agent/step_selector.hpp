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
#include "grampcd/util/class_forwarding.hpp"
#include "grampcd/optim/optim_util.hpp"

namespace grampcd
{
	class StepSelector
	{
	public:

		~StepSelector();

		/*Executes an ADMM Step*/
		virtual void execute_admmStep(const ADMMStep& step) = 0;
		/*returns the current ADMM Iterations*/
		virtual const unsigned int get_admmIter() = 0;
		
		// ADMM Iterations
		unsigned int admmIter_ = 0;

	protected:
		

	};

}