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
	class AsyncStepSelector : public StepSelector
	{
	public:

		AsyncStepSelector(SolverLocalPtr& local_solver, Agent* agent);

		/*executes an ADMM Step*/
		void execute_admmStep(const ADMMStep& step) override;
		/*returns the current ADMM iterations*/
		const unsigned int get_admmIter() override;
		/*sets the flag to stop the ADMM algorithm*/
		void set_flagStopAdmm( bool flag);
		/*returns the flag if the algorithm has stopped*/
		const bool get_flagStopAdmm();
		/*checks the delays of the neighbors*/
		const bool check_delays(const ADMMStep& step);
		/*checks how to continue the ADMM algorithm*/
		void check_continuation();
		/*evaluates the criterion for the delay*/
		const int check_criterionRegardingDelay();
		/*checks if the admm iterations are finished*/
		void check_admmIterations();
		/*checks if convergence has been achieved*/
		void check_convergence();

		// steps of ADMM algorithm
		ADMMStep previous_step_;
		ADMMStep current_step_;
		ADMMStep next_step_;


	private:
		// local Solver
		SolverLocalPtr local_solver_;
		// agent
		const Agent* agent_;

		// flag to stop admm
	    bool flagStopAdmm_;

	};

}
