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

#include "grampcd/state/agent_state.hpp"

#include "cereal/archives/binary.hpp"
#include "cereal/types/vector.hpp"

#include <iostream>

namespace grampcd
{
	class Solution
	{
	public:
		Solution();

		Solution(AgentState agent_state, AgentState predicted_agent_state);

		/*Update the solution with a simulated state.*/
		void update_state(const std::vector<typeRNum>& x, const typeRNum t0, const int id, const typeRNum cost);
		/*Update the solution with a predicted state.*/
		void update_predicted_state(const AgentState& state, const std::vector< typeRNum >& predicted_cost);
		/*Update the solution with debug cost.*/
		void update_debug_cost(const typeRNum cost);

		/*Simulated agent state*/
		AgentState agentState_;
		/*Predicted agent state*/
		AgentState predicted_agentState_;

		/*Simulated cost*/
		std::vector< typeRNum > cost_;
		/*Predicted cost*/
		std::vector< typeRNum > predicted_cost_;
		/*Debug cost*/
		std::vector< typeRNum > debug_cost_;

		/*Maximum number of data points*/
		unsigned int maximum_number_of_data_points_ = 0;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(agentState_, predicted_agentState_, cost_, predicted_cost_, debug_cost_, maximum_number_of_data_points_);
		}
	};

	std::ostream& operator<<(std::ostream& stream, const Solution& solution);
}