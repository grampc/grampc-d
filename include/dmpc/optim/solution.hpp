#ifndef SOLUTION_HPP
#define SOLUTION_HPP

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

#include "dmpc/state/agent_state.hpp"
#include <iostream>
#include <algorithm>

namespace dmpc
{
	class Solution
	{
	public:
		Solution();

		Solution(dmpc::AgentState agent_state, dmpc::AgentState predicted_agent_state);

		/*Update the solution with a simulated state.*/
		void update_state(const std::vector<typeRNum>& x, const typeRNum t0, const int id, const typeRNum cost);
		/*Update the solution with a predicted state.*/
		void update_predicted_state(const dmpc::AgentState& state, const std::vector< typeRNum >& predicted_cost);
		/*Update the solution with debug cost.*/
		void update_debug_cost(const typeRNum cost);

		/*Simulated agent state*/
		dmpc::AgentState agentState_;
		/*Predicted agent state*/
		dmpc::AgentState predicted_agentState_;

		/*Simulated cost*/
		std::vector< typeRNum > cost_;
		/*Predicted cost*/
		std::vector< typeRNum > predicted_cost_;
		/*Debug cost*/
		std::vector< typeRNum > debug_cost_;

		/*Maximum number of data points*/
		unsigned int maximum_number_of_data_points_ = 0;
	};

	typedef std::shared_ptr<Solution> SolutionPtr;

	std::ostream& operator<<(std::ostream& stream, const Solution& solution);
}

#endif // SOLUTION_HPP
