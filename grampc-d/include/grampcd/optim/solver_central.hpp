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

#include "grampcd/info/optimization_info.hpp"
#include "grampcd/optim/problem_description_central.hpp"

namespace grampcd
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

        const bool is_converged() const;

    private:
        OptimizationInfo optimizationInfo_;
	    std::vector<AgentPtr> agents_;
        ProblemDescriptionCentral problem_description_;
        std::shared_ptr<grampc::Grampc> solver_;
    };

}