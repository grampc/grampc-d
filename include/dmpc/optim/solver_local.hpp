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

#ifndef LOCAL_SOLVER_HPP
#define LOCAL_SOLVER_HPP

#include "dmpc/info/optimization_info.hpp"
#include "dmpc/optim/problem_description_local_default.hpp"
#include "dmpc/optim/problem_description_local_neighbor_approximation.hpp"

namespace dmpc
{

class SolverLocal
{
public:
    SolverLocal(Agent* agent, const OptimizationInfo& info, const LoggingPtr& log);

    void update_agentStates();
    void update_couplingStates();
    void update_multiplierStates();

    const bool is_converged() const;

    void initialize_ADMM();

    const std::vector<int>& get_x_index_xji() const;
    const std::vector<int>& get_u_index_uji() const;
    const std::vector<int>& get_u_index_xji() const;
    const std::vector<int>& get_u_index_vji() const;
    const int get_x_index_xji(int j) const;
    const int get_u_index_xji(int j) const;
    const int get_u_index_vji(int j) const;
    const int get_u_index_uji(int j) const;

protected:
	ProblemDescriptionLocalDefault default_problem_description_;
    ProblemDescriptionLocalNeighborApproximation neighbor_approximation_problem_description_;

	SolverPtr solver_;
    LoggingPtr log_;

	int ADMM_iter_ = 0;
	double ADMM_PrimalResiduum = 0.0;
	double ADMM_DualResiduum = 0.0;

    typeRNum adaptPenaltyParameter( typeRNum primal_residuum, typeRNum dual_residuum, typeRNum penalty );
	void penaltyParameterAdaption();

	Agent* agent_;
	OptimizationInfo info_;
};

typedef std::shared_ptr<SolverLocal> SolverLocalPtr;

}

#endif // LOCAL_SOLVER_HPP
