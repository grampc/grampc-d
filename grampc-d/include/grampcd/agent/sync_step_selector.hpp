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
#include"grampcd/agent/step_selector.hpp"


namespace grampcd
{
	class SyncStepSelector: public StepSelector
	{
	public:

		SyncStepSelector(SolverLocalPtr& local_solver);

		/*Executes an ADMM Step*/
		void execute_admmStep(const ADMMStep& step) override;
		/*returns the current ADMM Iterations*/
		const unsigned int get_admmIter() override;


	private:
		// local Solver
		SolverLocalPtr local_solver_;
	};

}
