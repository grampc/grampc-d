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

#ifndef CENTRALIZED_SOLVER_HPP
#define CENTRALIZED_SOLVER_HPP

#include "dmpc/info/optimization_info.hpp"
#include "dmpc/optim/problem_description_central.hpp"

namespace dmpc
{

/**
 * @brief Solver for the centralized MPC problem.
 */
class SolverCentral
{
public:
    SolverCentral(const std::vector<AgentPtr>& agents,
                      const OptimizationInfo& info);

    /*Solve the centralized OCP.*/
    void solve();

    /*Returns a vector with pointer to agents.*/
    const std::vector<AgentPtr>& get_agents() const;

private:
    OptimizationInfo optimizationInfo_;
	std::vector<AgentPtr> agents_;
    ProblemDescriptionCentral problem_description_;
    std::shared_ptr<grampc::Grampc> solver_;
};

typedef std::shared_ptr<SolverCentral> SolverCentralPtr;

}

#endif // CENTRALIZED_SOLVER_HPP
