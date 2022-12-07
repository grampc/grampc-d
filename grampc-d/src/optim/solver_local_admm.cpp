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

#include "grampcd/agent/agent.hpp"
#include "grampcd/agent/neighbor.hpp"

#include "grampcd/optim/solver_local_admm.hpp"
#include "grampcd/optim/optim_util.hpp"
#include "grampcd/optim/approximate_neighbor.hpp"

#include "grampcd/model/agent_model.hpp"

#include "grampcd/util/logging.hpp"

#include "grampcd/agent/step_selector.hpp"
#include "grampcd/agent/async_step_selector.hpp"

#include "grampcd/optim/solution.hpp"
#include <cmath>
#include <algorithm>

namespace grampcd
{

    SolverLocalADMM::SolverLocalADMM(Agent* agent, const OptimizationInfo& info, const LoggingPtr& log, const CommunicationInterfacePtr& communication_interface)
        : agent_(agent),
          default_problem_description_(agent),
          neighbor_approximation_problem_description_(agent, info),
          solver_(new grampc::Grampc(&default_problem_description_)),
          info_(info),
          log_(log),
          communication_interface_(communication_interface)
    {
	    if (info.APPROX_ApproximateDynamics_)
		    solver_.reset(new grampc::Grampc(&neighbor_approximation_problem_description_));
	    else
		    solver_.reset(new grampc::Grampc(&default_problem_description_));

        configureSolver(solver_, info);
    }

    void SolverLocalADMM::initialize_ADMM()
    {
        // needed for neighbor approximation
        agent_->set_neighbors_initial_states();
        // steps needed for asynchronous execution
        send_numberofNeighbors();    
        agent_->reset_stopAdmmflag_of_neighbors();
        agent_->initialize_allNeighborDelays();
    
    }

    void SolverLocalADMM::update_agentStates()
    {
        const unsigned int Nx = solver_->getParameters()->Nx;
        const unsigned int Nu = solver_->getParameters()->Nu;
        const unsigned int Nxi = agent_->get_Nxi();
        const unsigned int Nui = agent_->get_Nui();
        const unsigned int Nhor = solver_->getOptions()->Nhor;

        // set initial state x0 and initial control
        {
            const AgentState& state = agent_->get_agentState();

            // set initial time
            solver_->setparam_real("t0", state.t0_);

            // set initial state
            for(unsigned int k = 0; k < Nxi; ++k)
                solver_->getParameters()->x0[k] = state.x_[k];

            // set initial control trajectory
            for(unsigned int i = 0; i < Nhor; ++i)
            {
                for(unsigned int k = 0; k < Nui; ++k)
                    solver_->getWorkspace()->u[i * Nu + k] = state.u_[i * Nui + k];
            }
        }

        // consider local copies
        for(const auto& neighbor : agent_->get_neighbors())
        {
            if( neighbor->is_sendingNeighbor() || neighbor->is_approximating() )
            {
                const AgentState& state = neighbor->get_localCopies();
                const unsigned int j = neighbor->get_id();
                const unsigned int Nxj = neighbor->get_Nxj();
                const unsigned int Nuj = neighbor->get_Nuj();

                // set trajectory of local copies as initial control trajectory
                for(unsigned int i = 0; i < Nhor; ++i)
                {
                    for(unsigned int k = 0; k < Nuj; ++k)
                        solver_->getWorkspace()->u[i * Nu + get_u_index_uji(j) + k] = state.u_[i * Nuj + k];

                    if(agent_->is_approximatingDynamics())
                    {
                        for(unsigned int k = 0; k < Nxj; ++k)
                            solver_->getWorkspace()->u[i * Nu + get_u_index_vji(j) + k] = state.v_[i * Nxj + k];
                    }
                    else
                    {
                        for(unsigned int k = 0; k < Nxj; ++k)
                            solver_->getWorkspace()->u[i * Nu + get_u_index_xji(j) + k] = state.x_[i * Nxj + k];
                    }
                }

                // set initial state
                if(agent_->is_approximatingDynamics())
                {
                    for(unsigned int k = 0; k < Nxj; ++k)
                        solver_->getParameters()->x0[k + get_x_index_xji(j)] = state.x_[k];
                }
            }
        }

        // construct control limits
        std::vector<typeRNum> umin(solver_->getParameters()->Nu, -INF);
        std::vector<typeRNum> umax(solver_->getParameters()->Nu, INF);

        agent_->get_agentModel()->get_controlLimits(umin, umax);

        if( info_.APPROX_ApproximateConstraints_ )
        {
            for( const auto& neighbor : agent_->get_neighbors() )
            {
                const unsigned int j = neighbor->get_id();
                const unsigned int Nuj = neighbor->get_Nuj();
                const unsigned int u_index_uji = get_u_index_uji(j);

                std::vector<typeRNum> umin_neighbor(Nuj, -INF);
                std::vector<typeRNum> umax_neighbor(Nuj, INF);

                neighbor->get_agentModel()->get_controlLimits( umin_neighbor, umax_neighbor );

                std::copy(umin_neighbor.begin(), umin_neighbor.end(), &umin[ u_index_uji ]);
                std::copy(umax_neighbor.begin(), umax_neighbor.end(), &umax[ u_index_uji ]);
            }
        }

        solver_->setparam_real_vector("umin", &umin[0]);
        solver_->setparam_real_vector("umax", &umax[0]);

        // solve
        solver_->run();

        // update x_i and u_i
        AgentState state = agent_->get_agentState();
        for (unsigned int i = 0; i < Nhor; ++i)
        {
            // update state trajectory x_i
            for (unsigned int k = 0; k < Nxi; ++k)
                state.x_[i * Nxi + k] = solver_->getWorkspace()->x[i * Nx + k];

            // update control trajectory u_i
            for (unsigned int k = 0; k < Nui; ++k)
                state.u_[i * Nui + k] = solver_->getWorkspace()->u[i * Nu + k];
        }
        agent_->set_agentState(state);

        // update predicted state and control trajectories of local copies
        for(const auto& neighbor : agent_->get_neighbors())
        {
            if(neighbor->is_sendingNeighbor() || agent_->is_approximating())
            {
                const unsigned int j = neighbor->get_id();
                const unsigned int Nxj = neighbor->get_Nxj();
                const unsigned int Nuj = neighbor->get_Nuj();

                state = neighbor->get_localCopies();
                for(unsigned int i = 0; i < Nhor; ++i)
                {
                    // update u_{ji}
                    for(unsigned int k = 0; k < Nuj; ++k)
                        state.u_[i * Nuj + k] = solver_->getWorkspace()->u[i * Nu + get_u_index_uji(j) + k];

                    if(agent_->is_approximatingDynamics())
                    {
                        for(unsigned int k = 0; k < Nxj; ++k)
                        {
                            // update x_{ji} as state trajectory
                            state.x_[i * Nxj + k] = solver_->getWorkspace()->x[i * Nx + get_x_index_xji(j) + k];

                            // update v_{ji}
                            state.v_[i * Nxj + k] = solver_->getWorkspace()->u[i * Nu + get_u_index_vji(j) + k];
                        }
                    }
                    else
                    {
                        // update x_{ji} as control trajectory
                        for(unsigned int k = 0; k < Nxj; ++k)
                            state.x_[i * Nxj + k] = solver_->getWorkspace()->u[i * Nu + get_u_index_xji(j) + k];
                    }
                }

                neighbor->set_localCopies(state);
            }
        }

        // update predicted state for external influence
        if(agent_->is_approximatingDynamics())
        {
            const auto& t = agent_->get_agentState().t_;

            for(const auto& neighbor : agent_->get_neighbors())
            {
                state = neighbor->get_externalInfluence_agentState();

                for(unsigned int i=0; i < Nhor; ++i)
                {
                    std::vector<typeRNum> vij( Nxi, 0.0 );

                    // evaluate v
                    neighbor->get_neighborApproximation()->vfct( &vij[0], t[i], solver_->getWorkspace()->x + i*Nx, solver_->getWorkspace()->u + i*Nu,
                        get_x_index_xji(), get_u_index_uji(), get_u_index_vji());

                    // update trajectory for external influence
                    for(unsigned int k = 0; k < Nxi; ++k)
                        state.v_[i * Nxi + k] = vij[k];
                }

                neighbor->set_externalInfluence_agentState(state);
            }
        }
    }

    void SolverLocalADMM::update_couplingStates()
    {
        // consider coupling state zx_i and zu_i
        {
            CouplingState coupling_state = agent_->get_couplingState();
            std::vector<typeRNum> normalization_factor_x(coupling_state.z_x_.size(), 0.0);
            std::vector<typeRNum> normalization_factor_u(coupling_state.z_u_.size(), 0.0);

            // consistency constraints for agent
            {
                const AgentState& state = agent_->get_agentState();
                const MultiplierState& multiplier = agent_->get_multiplierState();
                const PenaltyState& penalty = agent_->get_penaltyState();

                for(unsigned int k = 0; k < coupling_state.z_x_.size(); ++k)
                {
                    coupling_state.z_x_[k] = penalty.rho_x_[k] * state.x_[k] - multiplier.mu_x_[k];
                    normalization_factor_x[k] = penalty.rho_x_[k];
                }

                for(unsigned int k = 0; k < coupling_state.z_u_.size(); ++k)
                {
                    coupling_state.z_u_[k] = penalty.rho_u_[k] * state.u_[k] - multiplier.mu_u_[k];
                    normalization_factor_u[k] = penalty.rho_u_[k];
                }
            }

            // consistency constraints for receiving neighbors
            for(const NeighborPtr& neighbor : agent_->get_neighbors())
            {
                if( neighbor->is_receivingNeighbor() || neighbor->is_approximating() )
                {
                    const AgentState& neigh_local_copies = neighbor->get_neighbors_localCopies();
                    const MultiplierState& multiplier = neighbor->get_neighbors_coupled_multiplierState();
                    const PenaltyState& penalty = neighbor->get_neighbors_coupled_penaltyState();

                    // consistency constraints (zxj - x_{ji})
                    for(unsigned int k = 0; k < coupling_state.z_x_.size(); ++k)
                    {
                        coupling_state.z_x_[k] += penalty.rho_x_[k] * neigh_local_copies.x_[k] - multiplier.mu_x_[k];
                        normalization_factor_x[k] += penalty.rho_x_[k];
                    }

                    // consistency constraints (zu_j - u_{ji})
                    for(int k = 0; k < (int)coupling_state.z_u_.size(); ++k)
                    {
                        coupling_state.z_u_[k] += penalty.rho_u_[k] * neigh_local_copies.u_[k] - multiplier.mu_u_[k];
                        normalization_factor_u[k] += penalty.rho_u_[k];
                    }
                }
            }

            //normalize
            for(unsigned int k = 0; k < coupling_state.z_x_.size(); ++k)
                coupling_state.z_x_[k] /= normalization_factor_x[k];

            for(unsigned int k = 0; k < coupling_state.z_u_.size(); ++k)
                coupling_state.z_u_[k] /= normalization_factor_u[k];

            // set new coupling state
            agent_->set_couplingState(coupling_state);
        }

        if(agent_->is_approximatingDynamics())
        {
            for(const NeighborPtr& neighbor: agent_->get_neighbors())
            {
                CouplingState coupling = neighbor->get_externalInfluence_couplingState();
                std::vector<typeRNum> normalization_factor_v(coupling.z_v_.size(), 0.0);

                // consistency constraint (z_{v, ij} - v_{ij})
                {
				    const AgentState& state = neighbor->get_externalInfluence_agentState();
				    const MultiplierState& multiplier = neighbor->get_externalInfluence_multiplierState();
				    const PenaltyState& penalty = neighbor->get_externalInfluence_penaltyState();

				    for (unsigned int k = 0; k < coupling.z_v_.size(); ++k)
				    {
					    coupling.z_v_[k] = penalty.rho_v_[k] * state.v_[k] - multiplier.mu_v_[k];
					    normalization_factor_v[k] = penalty.rho_v_[k];
				    }
                }

                // consistency constraints (zv_{ji} - v_{ji})
                {
                    const AgentState& state = neighbor->get_neighbors_localCopies();
                    const MultiplierState& multiplier = neighbor->get_neighbors_coupled_multiplierState();
                    const PenaltyState& penalty = neighbor->get_neighbors_coupled_penaltyState();

                    for(unsigned int k = 0; k < coupling.z_v_.size(); ++k)
                    {
                        coupling.z_v_[k] += penalty.rho_v_[k] * state.v_[k] - multiplier.mu_v_[k];
                        normalization_factor_v[k] += penalty.rho_v_[k];
                    }
                }

                // normalize
                for(unsigned int k = 0; k < coupling.z_v_.size(); ++k)
                    coupling.z_v_[k] /= normalization_factor_v[k];

                // set new coupling
                neighbor->set_externalInfluence_couplingState(coupling);
            }
        }
    }

    void SolverLocalADMM::update_multiplierStates()
    {
        // adapt penalty parameter
        if( info_.ADMM_AdaptPenaltyParameter_ )
            penaltyParameterAdaption();

        MultiplierState multiplier;
        unsigned int k;

        {
            multiplier = agent_->get_multiplierState();
            const CouplingState& coupling = agent_->get_couplingState();
            const AgentState& state = agent_->get_agentState();
            const PenaltyState& penalty = agent_->get_penaltyState();

            // consistency constraints (zx_i - x_i)
            for(k = 0; k < multiplier.mu_x_.size(); ++k)
                multiplier.mu_x_[k] += penalty.rho_x_[k] * (coupling.z_x_[k] - state.x_[k]);

            // consistency constraints (zu_i - u_i)
            for(k = 0; k < multiplier.mu_u_.size(); ++k)
                multiplier.mu_u_[k] += penalty.rho_u_[k] * (coupling.z_u_[k] - state.u_[k]);

            // set new multiplier state
            agent_->set_multiplierState(multiplier);
        }

        // consistency constraints for neighbors
        for(const auto& neighbor : agent_->get_neighbors())
        {
            {
                multiplier = neighbor->get_coupled_multiplierState();
                const AgentState& local_copies = neighbor->get_localCopies();
                const CouplingState& coupling = neighbor->get_neighbors_couplingState();
                const PenaltyState& penalty = neighbor->get_coupled_penaltyState();

                // consistency constraints (zx_j - x_{ji})
                for (k = 0; k < multiplier.mu_x_.size(); ++k)
                    multiplier.mu_x_[k] += penalty.rho_x_[k] * (coupling.z_x_[k] - local_copies.x_[k]);

                // consistency constraints (zu_j - u_{ji})
                for (k = 0; k < multiplier.mu_u_.size(); ++k)
                    multiplier.mu_u_[k] += penalty.rho_u_[k] * (coupling.z_u_[k] - local_copies.u_[k]);

                // consistency constraints (zv_{ji} - v_{ji})
                const CouplingState& neighbors_externalInfluence_couplingState = neighbor->get_neighbors_externalInfluence_couplingState();
                for (k = 0; k < multiplier.mu_v_.size(); ++k)
                    multiplier.mu_v_[k] += penalty.rho_v_[k] * (neighbors_externalInfluence_couplingState.z_v_[k] - local_copies.v_[k]);

                neighbor->set_coupled_multiplierState(multiplier);
            }

            {
                multiplier = neighbor->get_externalInfluence_multiplierState();
                const AgentState& state = neighbor->get_externalInfluence_agentState();
                const CouplingState& coupling = neighbor->get_externalInfluence_couplingState();
                const PenaltyState& penalty = neighbor->get_externalInfluence_penaltyState();

                // consistency constraints (zv_{ij} - v_{ij})
                for (k = 0; k < multiplier.mu_v_.size(); ++k)
                    multiplier.mu_v_[k] += penalty.rho_v_[k] * (coupling.z_v_[k] - state.v_[k]);

                neighbor->set_externalInfluence_multiplierState(multiplier);
            }
        }
    }

    void SolverLocalADMM::send_agentStates()
    {
        for (const auto& neighbor : agent_->get_neighbors())
        {
            if (neighbor->is_sendingNeighbor() || neighbor->is_approximating())
            {
                if (neighbor->is_approximatingDynamics())
                {
                    AgentState local_copies = neighbor->get_localCopies();
                    local_copies.x_.clear();

                    communication_interface_->send_localCopies(local_copies, agent_->get_id(), neighbor->get_id());
                }
                else
                    communication_interface_->send_localCopies(neighbor->get_localCopies(), agent_->get_id(), neighbor->get_id());
            }
        }
    }

    void SolverLocalADMM::send_couplingStates()
    {

        for (const auto& neighbor : agent_->get_neighbors())
        {
            // send coupling state
            if ((neighbor->is_receivingNeighbor() || neighbor->is_approximating()) && !neighbor->is_approximatingDynamics())
                communication_interface_->send_couplingState(agent_->get_couplingState(), agent_->get_id(), neighbor->get_id());

            // send coupling state and ext_influence_coupling_state
            if (neighbor->is_approximatingDynamics())
                communication_interface_->send_couplingState(agent_->get_couplingState(), neighbor->get_externalInfluence_couplingState(), agent_->get_id(), neighbor->get_id());
        }


    }

    void SolverLocalADMM::send_multiplierStates()
    {
        for (const auto& neighbor : agent_->get_neighbors())
        {
            if (neighbor->is_sendingNeighbor() || neighbor->is_approximating())
                communication_interface_->send_multiplierState(neighbor->get_coupled_multiplierState(),
                    neighbor->get_coupled_penaltyState(), agent_->get_id(), neighbor->get_id());
        }
    }

    void SolverLocalADMM::send_convergenceFlag()
    {
        communication_interface_->send_convergenceFlag(is_converged(), agent_->get_id());
    }

    void SolverLocalADMM::send_numberofNeighbors()
    {
        for (const auto& neighbor : agent_->get_neighbors())
            communication_interface_->send_numberOfNeighbors(static_cast<int>(agent_->get_neighbors().size()), agent_->get_id(), neighbor->get_id());
    }

    void SolverLocalADMM::send_stoppedAlgFlag()
    {
        // get and cast step selector 
        AsyncStepSelectorPtr stepselector =  std::static_pointer_cast<AsyncStepSelector>(agent_->get_stepSelector());

        // send flag to neighbors
        for (const auto& neighbor : agent_->get_neighbors())
        {
            communication_interface_->send_stoppedAlgFlag(stepselector->get_flagStopAlg(), agent_->get_id(), neighbor->get_id());
        }

        // send flag to Coordinator 
        communication_interface_->send_stoppedAlgFlag(stepselector->get_flagStopAlg(), agent_->get_id());
    }

    void SolverLocalADMM::print_debugCost()
    {
        agent_->get_solution()->update_debug_cost(agent_->get_predicted_cost());
    }

    void SolverLocalADMM::penaltyParameterAdaption()
    {
        typeRNum primal_residuum = 0.0;
        typeRNum dual_residuum = 0.0;
        PenaltyState penalty;
        unsigned int k = 0;

        {
            const AgentState& state = agent_->get_agentState();
            const CouplingState& coupling = agent_->get_couplingState();
            const CouplingState& previous_coupling = agent_->get_previous_couplingState();
            penalty = agent_->get_penaltyState();

            //*************************************
            // consistency constraints (zx_i - x_i)
            //*************************************
            for (k = 0; k < penalty.rho_x_.size(); ++k)
            {
                primal_residuum = std::abs(coupling.z_x_[k] - state.x_[k]);
                dual_residuum = std::abs(penalty.rho_x_[k] * (coupling.z_x_[k] - previous_coupling.z_x_[k]));

                penalty.rho_x_[k] = adaptPenaltyParameter(primal_residuum, dual_residuum, penalty.rho_x_[k]);
            }

            //*************************************
            // consistency constraints (zu_i - u_i)
            //*************************************
            for (k = 0; k < penalty.rho_u_.size(); ++k)
            {
                primal_residuum = std::abs(coupling.z_u_[k] - state.u_[k]);
                dual_residuum = std::abs(penalty.rho_u_[k] * (coupling.z_u_[k] - previous_coupling.z_u_[k]));

                penalty.rho_u_[k] = adaptPenaltyParameter(primal_residuum, dual_residuum, penalty.rho_u_[k]);
            }

            agent_->set_penaltyState(penalty);
        }

        for( const auto& neighbor : agent_->get_neighbors() )
        {
            {
                penalty = neighbor->get_coupled_penaltyState();
                const AgentState& local_copies = neighbor->get_localCopies();
                const CouplingState& coupling = neighbor->get_neighbors_couplingState();
                const CouplingState& previous_coupling = neighbor->get_previous_neighbors_couplingState();

                //*************************************
                // consistency constraints (zx_j - x_{ji})
                //*************************************
                for(k = 0; k < penalty.rho_x_.size(); ++k)
                {
                    primal_residuum = std::abs( coupling.z_x_[k] - local_copies.x_[k] );
                    dual_residuum = std::abs( penalty.rho_x_[k] * (coupling.z_x_[k] - previous_coupling.z_x_[k]) );

                    penalty.rho_x_[k] = adaptPenaltyParameter(primal_residuum, dual_residuum, penalty.rho_x_[k]);
                }

                //*************************************
                // consistency constraints (zu_j - u_{ji})
                //*************************************
                for(k = 0; k < penalty.rho_u_.size(); ++k)
                {
                    primal_residuum = std::abs( coupling.z_u_[k] - local_copies.u_[k] );
                    dual_residuum = std::abs(penalty.rho_u_[k] * (coupling.z_u_[k] - previous_coupling.z_u_[k]));

                    penalty.rho_u_[k] = adaptPenaltyParameter(primal_residuum, dual_residuum, penalty.rho_u_[k]);
                }

                //*************************************
                // consistency constraints (zv_{ji} - v_{ji})
                //*************************************
                const CouplingState& coupling_externalInfluence = neighbor->get_neighbors_externalInfluence_couplingState();
                const CouplingState& coupling_previous_externalInfluence = neighbor->get_previous_neighbors_externalInfluence_couplingState();

                for(k = 0; k < penalty.rho_v_.size(); ++k)
                {
                    primal_residuum = std::abs(coupling_externalInfluence.z_v_[k] - local_copies.v_[k] );
                    dual_residuum = std::abs(penalty.rho_v_[k] * (coupling_externalInfluence.z_v_[k] - coupling_previous_externalInfluence.z_v_[k]));

                    penalty.rho_v_[k] = adaptPenaltyParameter(primal_residuum, dual_residuum, penalty.rho_v_[k]);
                }

                neighbor->set_coupled_penaltyState(penalty);
            }

            {
                penalty = neighbor->get_externalInfluence_penaltyState();
                const CouplingState& coupling = neighbor->get_externalInfluence_couplingState();
                const CouplingState& previous_coupling = neighbor->get_previous_externalInfluence_couplingState();
                const AgentState& state = neighbor->get_externalInfluence_agentState();

                //*************************************
                // consistency constraints (zv_{ij} - v_{ij})
                //*************************************
                for(k = 0; k < penalty.rho_v_.size(); ++k)
                {
                    primal_residuum = std::abs( coupling.z_v_[k] - state.v_[k] );
                    dual_residuum = std::abs(penalty.rho_v_[k] * (coupling.z_v_[k] - previous_coupling.z_v_[k]));

                    penalty.rho_v_[k] = adaptPenaltyParameter(primal_residuum, dual_residuum, penalty.rho_v_[k]);
                }

               neighbor->set_externInfluence_penaltyState(penalty);
            }
        }
    }

    typeRNum SolverLocalADMM::adaptPenaltyParameter( typeRNum primal_residuum, typeRNum dual_residuum, typeRNum penalty )
    {
	    typeRNum factor;
	    if (dual_residuum > 1e-10)
            factor = primal_residuum / dual_residuum;
	    else
            factor = 1;

	    factor = std::min(factor, info_.ADMM_PenaltyIncreaseFactor_);
	    factor = std::max(factor, info_.ADMM_PenaltyDecreaseFactor_);

	    typeRNum adaptedPenaltyParameter = factor * penalty;
	    adaptedPenaltyParameter = std::min(adaptedPenaltyParameter, info_.ADMM_PenaltyMax_);
	    adaptedPenaltyParameter = std::max(adaptedPenaltyParameter, info_.ADMM_PenaltyMin_);

	    return adaptedPenaltyParameter;
    }

    const bool SolverLocalADMM::is_converged() 
    {
	    typeRNum primal_residuum = 0;
	    size_t cnt = 0;
        unsigned int k = 0;

	    {
		    const AgentState& state = agent_->get_agentState();
		    const CouplingState& coupling = agent_->get_couplingState();

		    for ( k = 0; k < coupling.z_x_.size(); ++k)
                primal_residuum += std::pow(coupling.z_x_[k] - state.x_[k], 2);
		    cnt += coupling.z_x_.size();

            for (k = 0; k < coupling.z_u_.size(); ++k)
                primal_residuum += std::pow(coupling.z_u_[k] - state.u_[k], 2);
		    cnt += coupling.z_u_.size();
	    }
	    for (const auto& neighbor : agent_->get_neighbors())
	    {
		    if (neighbor->is_sendingNeighbor() || neighbor->is_approximating())
		    {
			    const AgentState& local_copies = neighbor->get_localCopies();
			    const CouplingState& neighbors_coupling = neighbor->get_neighbors_couplingState();

				for (k = 0; k < neighbors_coupling.z_x_.size(); ++k)
					primal_residuum += std::pow(neighbors_coupling.z_x_[k] - local_copies.x_[k], 2);
			    cnt += neighbors_coupling.z_x_.size();

				for (k = 0; k < neighbors_coupling.z_u_.size(); ++k)
					primal_residuum += std::pow(neighbors_coupling.z_u_[k] - local_copies.u_[k], 2);
			    cnt += neighbors_coupling.z_u_.size();
		    }
		    if (neighbor->is_approximatingDynamics())
		    {
			    const CouplingState& ext_infl_coupling = neighbor->get_externalInfluence_couplingState();
			    const AgentState& ext_infl_state = neighbor->get_externalInfluence_agentState();

				for (k = 0; k < ext_infl_coupling.z_v_.size(); ++k)
					primal_residuum += std::pow(ext_infl_coupling.z_v_[k] - ext_infl_state.v_[k], 2);
			    cnt += ext_infl_coupling.z_v_.size();

			    const CouplingState& neigh_ext_infl_coupling = neighbor->get_neighbors_externalInfluence_couplingState();
			    const AgentState& local_copies = neighbor->get_localCopies();

				for (k = 0; k < neigh_ext_infl_coupling.z_v_.size(); ++k)
					primal_residuum += std::pow(neigh_ext_infl_coupling.z_v_[k] - local_copies.v_[k], 2);
			    cnt += neigh_ext_infl_coupling.z_v_.size();
		    }
	    }
	    primal_residuum = primal_residuum / cnt;

	    return primal_residuum < std::pow(info_.ADMM_ConvergenceTolerance_, 2);
    }

    const std::vector<int>& SolverLocalADMM::get_u_index_uji() const
    {
        if(agent_->is_approximatingDynamics())
            return neighbor_approximation_problem_description_.get_u_index_uji();
        else
            return default_problem_description_.get_u_index_uji();
    }

    const std::vector<int>& SolverLocalADMM::get_u_index_vji() const
    {
        return neighbor_approximation_problem_description_.get_u_index_vji();
    }

    const std::vector<int>& SolverLocalADMM::get_u_index_xji() const
    {
        return default_problem_description_.get_u_index_xji();
    }

    const std::vector<int>& SolverLocalADMM::get_x_index_xji() const
    {
        return neighbor_approximation_problem_description_.get_x_index_xji();
    }

    const int SolverLocalADMM::get_u_index_uji(int j) const
    {
        if(agent_->is_approximatingDynamics())
            return neighbor_approximation_problem_description_.get_u_index_uji(j);
        else
            return default_problem_description_.get_u_index_uji(j);
    }

    const int SolverLocalADMM::get_u_index_vji(int j) const
    {
        if(agent_->is_approximatingDynamics())
            return neighbor_approximation_problem_description_.get_u_index_vji(j);
        else
        {
            log_->print(DebugType::Error) << "[SolverLocalADMM::get_u_index_vji]: "
                << "If dynamics are not approximated, there is no external influence." << std::endl;
            return 0;
        }
    }

    const int SolverLocalADMM::get_u_index_xji(int j) const
    {
        return default_problem_description_.get_u_index_xji(j);
    }

    const int SolverLocalADMM::get_x_index_xji(int j) const
    {
        if(agent_->is_approximatingDynamics())
		    return neighbor_approximation_problem_description_.get_x_index_xji(j);
	    else
	    {
		    log_->print(DebugType::Error) << "[SolverLocalADMM::get_x_index_xji]: "
			    << "If dynamics are not approximated, there is no extended state." << std::endl;
		    return 0;
	    }
    }
   
}
