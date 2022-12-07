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
#include"grampcd/agent/step_selector.hpp"


namespace grampcd
{
	class SyncStepSelector: public StepSelector
	{
	public:

		SyncStepSelector(SolverLocalPtr& local_solver, Agent* agent);

		/*Executes an ADMM Step*/
		void execute_algStep(const AlgStep& step) override;
		/*returns the current ADMM Iterations*/
		const unsigned int get_algIter() override;


	private:
		// local Solver
		SolverLocalPtr local_solver_;
		// agent 
		const Agent* agent_;
	};

}
