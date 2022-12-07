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
	class AsyncStepSelector : public StepSelector
	{
	public:

		AsyncStepSelector(SolverLocalPtr& local_solver, Agent* agent);

		/*executes an ADMM Step*/
		void execute_algStep(const AlgStep& step) override;
		/*returns the current ADMM iterations*/
		const unsigned int get_algIter() override;
		/*sets the flag to stop an algorithm (ADMM or Sensi)*/
		void set_flagStopAlg( bool flag);
		/*returns the flag if the algorithm has stopped*/
		const bool get_flagStopAlg();
		/*checks the delays of the neighbors*/
		const bool check_delays(const AlgStep& step);
		/*checks how to continue the ADMM algorithm*/
		void check_continuation();
		/*evaluates the criterion for the delay*/
		const int check_criterionRegardingDelay();
		/*checks if the admm iterations are finished*/
		void check_algIterations();
		/*checks if convergence has been achieved*/
		void check_convergence();

		// steps of ADMM algorithm
		AlgStep previous_step_;
		AlgStep current_step_;
		AlgStep next_step_;


	private:
		// local Solver
		SolverLocalPtr local_solver_;
		// agent
		const Agent* agent_;

		// flag to stop Algorithm (ADMM or Sensi)
	    bool flagStopAlg_;

	};

}
