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

#include "grampcd/optim/solution.hpp"

#include <algorithm>

namespace grampcd
{
	Solution::Solution() {}

	Solution::Solution(AgentState agent_state, AgentState predicted_agent_state)
	{
		agentState_ = agent_state;
		predicted_agentState_ = predicted_agent_state;
	}

	void Solution::update_state(const std::vector<typeRNum>& x, const typeRNum t0, const int id, const typeRNum cost)
	{

		// if predicted agent state is not set, this agent entered during runtime.
		if (predicted_agentState_.t_.empty())
			return;

		const auto Nu = static_cast<unsigned int>(predicted_agentState_.u_.size() / predicted_agentState_.t_.size());
		const auto Nx = static_cast<unsigned int>(x.size());

		agentState_.i_ = id;

		// check if maximum number of stored data points should be considered AND
		// maximum number of data points is already reached
		if (maximum_number_of_data_points_ > 0 && 
			agentState_.t_.size() == maximum_number_of_data_points_)
		{
			// shift time
			for (unsigned int i = 1; i < agentState_.t_.size(); ++i)
				agentState_.t_[i - 1] = agentState_.t_[i];
			// set time
			agentState_.t_.back() = t0;

			// shift controls
			for (unsigned int i = Nu; i < agentState_.u_.size(); ++i)
				agentState_.u_[i - Nu] = agentState_.u_[i];
			// set controls
			// As the predicted controls are used for the simulation, they can be 
			// directly used here.
			for (unsigned int i = 0; i < Nu; ++i)
				agentState_.u_[agentState_.u_.size() - Nu + i] = predicted_agentState_.u_[i];

			// shift states
			for (unsigned int i = Nx; i < agentState_.x_.size(); ++i)
				agentState_.x_[i - Nu] = agentState_.x_[i];
			// set states
			for (unsigned int i = 0; i < Nx; ++i)
				agentState_.x_[agentState_.x_.size() - Nx + i] = x[i];

			// shift cost
			for (unsigned int i = 1; i < cost_.size(); ++i)
				cost_[i - 1] = cost_[i];
			// set cost
			cost_.back() = cost;
		}
		else
		{
			// time
			agentState_.t_.push_back(t0);

			// state
			for (unsigned int i = 0; i < Nx; ++i)
				agentState_.x_.push_back(x[i]);

			// control
			for (unsigned int i = 0; i < Nu; ++i)
				agentState_.u_.push_back(predicted_agentState_.u_[i]);

			// cost
			cost_.push_back(cost);
		}
	}

	void Solution::update_predicted_state(const AgentState& state, const std::vector< typeRNum >& predicted_cost)
	{
		predicted_agentState_ = state;
		predicted_cost_ = predicted_cost;
	}

	void Solution::update_debug_cost(const typeRNum cost)
	{
		debug_cost_.push_back(cost);
	}

	std::ostream& operator<<(std::ostream& stream, const Solution& solution)
	{
		unsigned int Nx = 0;
		unsigned int Nu = 0;

		if (solution.agentState_.t_.size() > 0)
		{
			Nx = static_cast<int>(solution.agentState_.x_.size()) / static_cast<int>(solution.agentState_.t_.size());
			Nu = static_cast<int>(solution.agentState_.u_.size()) / static_cast<int>(solution.agentState_.t_.size());
		}

		size_t max_rows = std::max(solution.agentState_.t_.size(), solution.debug_cost_.size());

		stream << "AgentState_t\t";
		for (unsigned int i = 0; i < Nx; ++i)
			stream << "AgentState_x" << i << "\t";
		for (unsigned int i = 0; i < Nu; ++i)
			stream << "AgentState_u" << i << "\t";
		for (unsigned int i = 0; i < Nx; ++i)
			stream << "AgentState_v" << i << "\t";
		stream << "Cost\t"
			<< "Debug_Cost" << std::endl;

		for (unsigned int k = 0; k < max_rows; ++k)
		{
			/*
			Print agent state
			*/

			// print time
			if (solution.agentState_.t_.size() > k)
				stream << solution.agentState_.t_[k] << "\t";
			else
				stream << "0 \t ";

			// print states
			for (unsigned int i = 0; i < Nx; ++i)
			{
				if (solution.agentState_.x_.size() > Nx * k + i)
					stream << solution.agentState_.x_[Nx * k + i] << "\t";
				else
					stream << "0 \t ";
			}

			for (unsigned int i = 0; i < Nu; ++i)
			{
				// print control
				if (solution.agentState_.u_.size() > Nu * k + i)
					stream << solution.agentState_.u_[Nu * k + i] << "\t";
				else
					stream << "0 \t ";
			}

			for (unsigned int i = 0; i < Nx; ++i)
			{
				// print external influence
				if (solution.agentState_.v_.size() > Nx * k + i)
					stream << solution.agentState_.v_[Nx * k + i] << "\t";
				else
					stream << "0 \t ";
			}

			// print cost
			if (solution.cost_.size() > k)
				stream << solution.cost_[k] << "\t";
			else
				stream << "0 \t ";

			/*
			// print debug cost
			*/

			if (solution.debug_cost_.size() > k)
				stream << solution.debug_cost_[k] << "\t";
			else
				stream << "0";

			// newline
			if (k < max_rows)
				stream << "\n";
			else
				stream << std::endl;
		}

		return stream;
	}
}