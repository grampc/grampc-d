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
#include "grampcd/optim/problem_description_local_default.hpp"
#include "grampcd/optim/problem_description_local_neighbor_approximation.hpp"
#include "grampcd/comm/communication_interface.hpp"
#include "grampcd/optim/solver_local.hpp"

namespace grampcd
{

    class SolverLocalADMM : public SolverLocal
    {
    public:
        SolverLocalADMM(Agent* agent, const OptimizationInfo& info, const LoggingPtr& log, const CommunicationInterfacePtr& communication_interface);


        /***********************************************
      Generic functions
      ************************************************/
      /*updates the agent States*/
        void update_agentStates() override;
        /*sends the agent states to receiving neighbors*/
        void send_agentStates() override;
        /*checks if the algorithm is converged*/
        const bool is_converged() override;
        /*prints the debug costs to the solution file*/
        void print_debugCost() override;
        /*Sends the convergence flag*/
        void send_convergenceFlag() override;
        /*sends a flag if the agents finishes its ADMM iterations*/
        void send_stoppedAlgFlag() override;

        /***********************************************
      ADMM functions
      ************************************************/
        /*initializes the ADMM algorithm*/
        void initialize_ADMM() override;
        /*updates the coupling states*/
        void update_couplingStates() override;
        /*updates the multiplier states*/
        void update_multiplierStates() override;
        /*sends the coupling states*/
        void send_couplingStates() override;
        /*sends the multiplier states*/
        void send_multiplierStates() override;
        /*sends the number of neighbors */
        void send_numberofNeighbors() override;
          

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
	    double ADMM_PrimalResiduum = 0.0;
	    double ADMM_DualResiduum = 0.0;

        typeRNum adaptPenaltyParameter( typeRNum primal_residuum, typeRNum dual_residuum, typeRNum penalty );
	    void penaltyParameterAdaption();

	    Agent* agent_;
	    OptimizationInfo info_;
        CommunicationInterfacePtr communication_interface_;
    };

}