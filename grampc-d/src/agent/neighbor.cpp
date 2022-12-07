/* This file is part of GRAMPC-D - (https://github.com/DanielBurk/GRAMPC-D.git)
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

#include "grampcd/optim/optim_util.hpp"
#include "grampcd/optim/approximate_neighbor.hpp"
#include "grampcd/optim/optim_util.hpp"

#include "grampcd/model/agent_model.hpp"
#include "grampcd/model/coupling_model.hpp"

#include "grampcd/info/optimization_info.hpp"

#include "grampcd/util/logging.hpp"

namespace grampcd
{

    Neighbor::Neighbor(const int id, const AgentModelPtr& agent_model, const LoggingPtr& log)
        :   id_(id),
            sending_neighbor_(false),
            receiving_neighbor_(false),
            agentModel_(agent_model),
            couplingModel_(nullptr),
            copied_couplingModel_(nullptr),
            approximate_neighbor_(nullptr),
            log_(log),
            agent_id_(-1),
            flag_StoppedAlg_(false)
    {
        sending_neighbor_ = false;
        receiving_neighbor_ = false;

        is_approximating_ = false;
        is_approximatingCost_ = false;
        is_approximatingDynamics_ = false;
        is_approximatingConstraints_ = false;

        delay_agentState_ = 0;
        delay_couplingState_ = 0;
        delay_multiplierState_ = 0;
        delay_sensiState_ = 0;
        delay_sensiAgentState_ = 0;
    }

    const int Neighbor::get_id() const
    {
        return id_;
    }

    const AgentModelPtr& Neighbor::get_agentModel() const
    {
        return agentModel_;
    }

    const CouplingModelPtr& Neighbor::get_couplingModel() const
    {
        return couplingModel_;
    }

    void Neighbor::set_couplingModel(const CouplingModelPtr& model)
    {
        couplingModel_ = model;
    }

    const CouplingModelPtr& Neighbor::get_copied_couplingModel() const
    {
        return copied_couplingModel_;
    }

    void Neighbor::set_copied_couplingModel(const CouplingModelPtr& model)
    {
        copied_couplingModel_ = model;
    }

    const ApproximateNeighborPtr Neighbor::get_neighborApproximation() const
    {
        return approximate_neighbor_;
    }

    const unsigned int Neighbor::get_numberOfNeighbors() const
    {
        return numberOfNeighbors_;
    }

    void Neighbor::set_numberOfNeighbors(const unsigned int number)
    {
        numberOfNeighbors_ = number;
    }

    const AgentState& Neighbor::get_externalInfluence_agentState() const
    {
        return  externalInfluence_agentState_;
    }

    const CouplingState& Neighbor::get_externalInfluence_couplingState() const
    {
        return externalInfluence_couplingState_;
    }

    const MultiplierState& Neighbor::get_externalInfluence_multiplierState() const
    {
        return externalInfluence_multiplierState_;
    }

    const PenaltyState& Neighbor::get_externalInfluence_penaltyState() const
    {
        return externalInfluence_penaltyState_;
    }

    const CouplingState& Neighbor::get_neighbors_externalInfluence_couplingState() const
    {
        return neighbors_externalInfluence_couplingState_;
    }

    const MultiplierState& Neighbor::get_neighbors_externalInfluence_multiplierState() const
    {
        return neighbors_externalInfluence_multiplierState_;
    }

    const PenaltyState& Neighbor::get_neighbors_externalInfluence_penaltyState() const
    {
        return neighbors_externalInfluence_penaltyState_;
    }

    const CouplingState& Neighbor::get_previous_externalInfluence_couplingState() const
    {
        return previous_externalInfluence_couplingState_;
    }

    const MultiplierState& Neighbor::get_previous_externalInfluence_multiplierState() const
    {
        return previous_externalInfluence_multiplierState_;
    }

    const AgentState& Neighbor::get_localCopies() const
    {
        return local_copies_;
    }

    const MultiplierState& Neighbor::get_coupled_multiplierState() const
    {
        return coupled_multiplierState_;
    }

    const PenaltyState& Neighbor::get_coupled_penaltyState() const
    {
        return coupled_penaltyState_;
    }

    const SensiState& Neighbor::get_sensiState() const
    {
        return sensiState_;
    }

    const ConstraintState& Neighbor::get_coupled_constraintState() const
    {
        return coupled_constraintState_;
    }

    const AgentState& Neighbor::get_neighbors_localCopies() const
    {
        return neighbors_localCopies_;
    }

    const AgentState& Neighbor::get_neighbors_agentState() const
    {
        return neighbors_agentState_;
    }

    const AgentState& Neighbor::get_neighbors_desiredAgentState() const
    {
        return neighbors_desiredAgentState_;
    }

    const CouplingState& Neighbor::get_neighbors_couplingState() const
    {
        return neighbors_couplingState_;
    }

    const MultiplierState& Neighbor::get_neighbors_coupled_multiplierState() const
    {
        return neighbors_coupled_multiplierState_;
    }

    const PenaltyState& Neighbor::get_neighbors_coupled_penaltyState() const
    {
        return neighbors_coupled_penaltyState_;
    }

    const ConstraintState& Neighbor::get_neighbors_coupled_constraintState() const
    {
        return neighbors_coupled_constraintState_;
    }

    const CouplingState& Neighbor::get_previous_neighbors_externalInfluence_couplingState() const
    {
        return previous_neighbors_externalInfluence_couplingState_;
    }

    const CouplingState& Neighbor::get_previous_neighbors_couplingState() const
    {
        return previous_neighbors_couplingState_;
    }

    const AgentState& Neighbor::get_desiredAgentState() const
    {
        return neighbors_desiredAgentState_;
    }

    const CouplingState& Neighbor::get_previous_couplingState() const
    {
        return previous_couplingState_;
    }

    const MultiplierState& Neighbor::get_previous_multiplierState() const
    {
        return previous_multiplierState_;
    }

    const bool Neighbor::is_approximating() const
    {
        return is_approximating_;
    }

    const bool Neighbor::is_approximatingCost() const
    {
        return is_approximatingCost_;
    }

    const bool Neighbor::is_approximatingConstraints() const
    {
        return is_approximatingConstraints_;
    }

    const bool Neighbor::is_approximatingDynamics() const
    {
        return is_approximatingDynamics_;
    }

    void Neighbor::set_externalInfluence_agentState( const AgentState& state )
    {
        if( compare_stateDimensions( externalInfluence_agentState_, state ) )
            externalInfluence_agentState_ = state;
        else
            log_->print(DebugType::Error) << "[Neighbor::set_externalInfluence_agentState] "
            << "Failed to set external influence agent state, as dimensions don't fit." << std::endl;
    }

    void Neighbor::set_externalInfluence_couplingState( const CouplingState& coupling )
    {
        if( compare_stateDimensions( externalInfluence_couplingState_, coupling ) )
        {
            previous_externalInfluence_couplingState_ = externalInfluence_couplingState_;
            externalInfluence_couplingState_ = coupling;
        }
        else
            log_->print(DebugType::Error) << "[Neighbor::set_externalInfluence_couplingState] "
            << "Failed to set external influence coupling state, as dimensions don't fit." << std::endl;
    }

    void Neighbor::set_externalInfluence_multiplierState(const MultiplierState &multiplier)
    {
        if( compare_stateDimensions( externalInfluence_multiplierState_, multiplier ) )
        {
            previous_externalInfluence_multiplierState_ = externalInfluence_multiplierState_;
            externalInfluence_multiplierState_ = multiplier;
        }
        else
            log_->print(DebugType::Error) << "[Neighbor::set_externalInfluence_multiplierState] "
            << "Failed to set external influence multiplier state, as dimensions don't fit." << std::endl;
    }

    void Neighbor::set_externInfluence_penaltyState( const PenaltyState& penalty )
    {
        if( compare_stateDimensions( externalInfluence_penaltyState_, penalty ) )
            externalInfluence_penaltyState_ = penalty;
        else
            log_->print(DebugType::Error) << "[Neighbor::set_externInfluence_penaltyState] "
            << "Failed to set external influence penalty state, as dimensions don't fit." << std::endl;
    }

    void Neighbor::set_neighbors_externalInfluence_penaltyState(const PenaltyState &penalty)
    {
        if( compare_stateDimensions( neighbors_externalInfluence_penaltyState_, penalty ) )
        {
            previous_neighbors_externalInfluence_couplingState_ = neighbors_externalInfluence_couplingState_;
            neighbors_externalInfluence_penaltyState_ = penalty;
        }
        else
            log_->print(DebugType::Error) << "[Neighbor::set_neighbors_externalInfluence_penaltyState] "
            << "Failed to set external influence penalty state, as dimensions don't fit." << std::endl;
    }

    void Neighbor::set_neighbors_externalInfluence_couplingState(const CouplingState &coupling)
    {
        if( compare_stateDimensions( neighbors_externalInfluence_couplingState_, coupling ) )
        {
            previous_neighbors_externalInfluence_couplingState_ = neighbors_externalInfluence_couplingState_;
            neighbors_externalInfluence_couplingState_ = coupling;
        }
        else
            log_->print(DebugType::Error) << "[Neighbor::set_neighbors_externalInfluence_couplingState] "
            << "Failed to set external influence coupling state, as dimensions don't fit." << std::endl;
    }

    void Neighbor::set_neighbors_externalInfluence_multiplierState(const MultiplierState &multiplier)
    {
        if( compare_stateDimensions( neighbors_externalInfluence_multiplierState_, multiplier ) )
            neighbors_externalInfluence_multiplierState_ = multiplier;
        else
            log_->print(DebugType::Error) << "[Neighbor::set_neighbors_externalInfluence_multiplierState] "
            << "Failed to set external influence multiplier state, as dimensions don't fit." << std::endl;
    }

    void Neighbor::set_localCopies(const AgentState &state)
    {
        if( compare_stateDimensions( local_copies_, state ) )
            local_copies_ = state;
        else
            log_->print(DebugType::Error) << "[Neighbor::set_localCopies] "
            << "Failed to set local copies, as dimensions don't fit." << std::endl;
    }

    void Neighbor::set_coupled_multiplierState(const MultiplierState &multiplier)
    {
        if( compare_stateDimensions( coupled_multiplierState_, multiplier ) )
        {
            previous_multiplierState_ = coupled_multiplierState_;
            coupled_multiplierState_ = multiplier;
        }
        else
            log_->print(DebugType::Error) << "[Neighbor::set_coupled_multiplierState] "
            << "Failed to set coupled multiplier state, as dimensions don't fit." << std::endl;
    }

    void Neighbor::set_coupled_penaltyState(const PenaltyState &penalty)
    {
        if( compare_stateDimensions( coupled_penaltyState_, penalty ) )
            coupled_penaltyState_ = penalty;
        else
            log_->print(DebugType::Error) << "[Neighbor::set_coupled_penaltyState] "
            << "Failed to set coupled penalty state, as dimensions don't fit." << std::endl;
    }

    void Neighbor::set_sensiState(const SensiState& state)
    {
        if (compare_stateDimensions(sensiState_, state))
            sensiState_ = state;
        else
            log_->print(DebugType::Error) << "[Neighbor::set_sensiState] "
            << "Failed to set sensi state, as dimensions don't fit." << std::endl;
    }

    void Neighbor::set_coupled_constraintState(const ConstraintState& state)
    {
        if (compare_stateDimensions(coupled_constraintState_, state))
            coupled_constraintState_ = state;
        else
            log_->print(DebugType::Error) << "[Neighbor::set_coupled_constraintState] "
            << "Failed to set coupled constraint state, as dimensions don't fit." << std::endl;
    }

    void Neighbor::set_neighbors_localCopies( const AgentState& state )
    {
        if( compare_stateDimensions( neighbors_localCopies_, state ) )
            neighbors_localCopies_ = state;
        else
            log_->print(DebugType::Error) << "[Neighbor::set_neighbors_localCopies] "
            << "Failed to set neighbors local copies, as dimensions don't fit." << std::endl;
    }

    void Neighbor::set_neighbors_agentState(const AgentState& state)
    {
        if (compare_stateDimensions(neighbors_agentState_, state))
            neighbors_agentState_ = state;
        else
            log_->print(DebugType::Error) << "[Neighbor::set_neighbors_agentState] "
            << "Failed to set neighbors agent states, as dimensions don't fit." << std::endl;
    }

    void Neighbor::set_neighbors_couplingState( const CouplingState& coupling )
    {
        if( compare_stateDimensions( neighbors_couplingState_, coupling ) )
        {
            previous_neighbors_couplingState_ = neighbors_couplingState_;
            neighbors_couplingState_ = coupling;
        }
        else
            log_->print(DebugType::Error) << "[Neighbor::set_neighbors_couplingState] "
            << "Failed to set neighbors coupling state, as dimensions don't fit." << std::endl;
    }

    void Neighbor::set_neighbors_coupled_multiplierState( const MultiplierState& multiplier )
    {
        if( compare_stateDimensions( neighbors_coupled_multiplierState_, multiplier ) )
            neighbors_coupled_multiplierState_ = multiplier;
        else
            log_->print(DebugType::Error) << "[Neighbor::set_neighbors_coupled_multiplierState] "
            << "Failed to set neighbors coupled multiplier state, as dimensions don't fit." << std::endl;
    }

    void Neighbor::set_neighbors_coupled_penaltyState( const PenaltyState& penalty )
    {
        if( compare_stateDimensions( neighbors_coupled_penaltyState_, penalty ) )
            neighbors_coupled_penaltyState_ = penalty;
        else
            log_->print(DebugType::Error) << "[Neighbor::set_neighbors_coupled_penaltyState] "
            << "Failed to set neighbors coupled penalty state, as dimensions don't fit." << std::endl;
    }

    void Neighbor::set_neighbors_coupled_constraintState(const ConstraintState& state)
    {
        if (compare_stateDimensions(neighbors_coupled_constraintState_, state))
            neighbors_coupled_constraintState_ = state;
        else
            log_->print(DebugType::Error) << "[Neighbor::set_neighbors_coupled_constraintState] "
            << "Failed to set neighbors coupled constraint state, as dimensions don't fit." << std::endl;
    }

    void Neighbor::set_neighbors_desiredAgentState( const AgentState& state )
    {
        if( compare_stateDimensions( neighbors_desiredAgentState_, state ) )
            neighbors_desiredAgentState_ = state;
        else
            log_->print(DebugType::Error) << "[Neighbor::set_neighbors_desiredAgentState] "
            << "Failed to neighbors desired agent state, as dimensions don't fit." << std::endl;
    }

    void Neighbor::initialize_neighbor(const std::vector<double>& t, const OptimizationInfo& optimization_info, Agent* agent)
    {
        is_approximatingConstraints_ = optimization_info.APPROX_ApproximateConstraints_;
        is_approximatingCost_ = optimization_info.APPROX_ApproximateCost_;
        is_approximatingDynamics_ = optimization_info.APPROX_ApproximateDynamics_;

        is_approximating_ = is_approximatingConstraints_ || is_approximatingCost_ || is_approximatingDynamics_;

        agent_id_ = agent->get_id();

        const unsigned int Nhor = optimization_info.COMMON_Nhor_;

        const auto Nui = get_Nui();
        const auto Nxi = get_Nxi();
        const auto Nuj = get_Nuj();
        const auto Nxj = get_Nxj();

        if( ( is_sendingNeighbor() || is_approximating() ) && ! is_approximatingDynamics() )
        {
            local_copies_.t_ = t;
            local_copies_.i_ = this->id_;
            local_copies_.x_.resize(Nxj*Nhor, 0.0);
            local_copies_.u_.resize(Nuj*Nhor, 0.0);
            local_copies_.lambda_.resize(Nxj * Nhor, 0.0);

            coupled_multiplierState_.t_ = t;
            coupled_multiplierState_.i_ = this->id_;
            coupled_multiplierState_.mu_x_.resize(Nxj*Nhor, 0.0);
            coupled_multiplierState_.mu_u_.resize(Nuj*Nhor, 0.0);

            coupled_penaltyState_.t_ = t;
            coupled_penaltyState_.i_ = this->id_;
            coupled_penaltyState_.rho_x_.resize(Nxj*Nhor, optimization_info.ADMM_PenaltyInit_);
            coupled_penaltyState_.rho_u_.resize(Nuj*Nhor, optimization_info.ADMM_PenaltyInit_);
                        
            neighbors_couplingState_.t_ = t;
            neighbors_couplingState_.i_ = this->id_;
            neighbors_couplingState_.z_x_.resize(Nxj*Nhor, 0.0);
            neighbors_couplingState_.z_u_.resize(Nuj*Nhor, 0.0);
                     
            previous_neighbors_couplingState_.t_ = t;
            previous_neighbors_couplingState_.i_ = this->id_;
            previous_neighbors_couplingState_.z_x_.resize(Nxj*Nhor, 0.0);
            previous_neighbors_couplingState_.z_u_.resize(Nuj*Nhor, 0.0);
        }

        if( ( is_receivingNeighbor() || is_approximating() ) && ! is_approximatingDynamics() )
        {
            neighbors_localCopies_.t_ = t;
            neighbors_localCopies_.i_ = agent->get_id();
            neighbors_localCopies_.x_.resize(Nxi*Nhor, 0.0);
            neighbors_localCopies_.u_.resize(Nui*Nhor, 0.0);
            neighbors_localCopies_.lambda_.resize(Nxi * Nhor, 0.0);

            neighbors_coupled_multiplierState_.t_ = t;
            neighbors_coupled_multiplierState_.i_ = agent->get_id();
            neighbors_coupled_multiplierState_.mu_x_.resize(Nxi*Nhor, 0.0);
            neighbors_coupled_multiplierState_.mu_u_.resize(Nui*Nhor, 0.0);

            neighbors_coupled_penaltyState_.t_ = t;
            neighbors_coupled_penaltyState_.i_ = agent->get_id();
            neighbors_coupled_penaltyState_.rho_x_.resize(Nxi*Nhor, optimization_info.ADMM_PenaltyInit_);
            neighbors_coupled_penaltyState_.rho_u_.resize(Nui*Nhor, optimization_info.ADMM_PenaltyInit_);

        }

        if( is_approximatingDynamics() )
        {
            local_copies_.t_ = t;
            local_copies_.i_ = this->id_;
            local_copies_.x_.resize(Nxj*Nhor, 0.0);
            local_copies_.u_.resize(Nuj*Nhor, 0.0);
            local_copies_.v_.resize(Nxj*Nhor, 0.0);
            local_copies_.lambda_.resize(Nxj * Nhor, 0.0);

            coupled_multiplierState_.t_ = t;
            coupled_multiplierState_.i_ = this->id_;
            coupled_multiplierState_.mu_u_.resize(Nuj*Nhor, 0.0);
            coupled_multiplierState_.mu_v_.resize(Nxj*Nhor, 0.0);

            coupled_penaltyState_.t_ = t;
            coupled_penaltyState_.i_ = this->id_;
            coupled_penaltyState_.rho_u_.resize(Nuj*Nhor, optimization_info.ADMM_PenaltyInit_);
            coupled_penaltyState_.rho_v_.resize(Nxj*Nhor, optimization_info.ADMM_PenaltyInit_);

            externalInfluence_agentState_.t_ = t;
            externalInfluence_agentState_.i_ = agent->get_id();
            externalInfluence_agentState_.v_.resize(Nxi*Nhor, 0.0);
            externalInfluence_agentState_.lambda_.resize(Nxi * Nhor, 0.0);

            externalInfluence_couplingState_.t_ = t;
            externalInfluence_agentState_.i_ = agent->get_id();
            externalInfluence_couplingState_.z_v_.resize(Nxi*Nhor, 0.0);

            externalInfluence_multiplierState_.t_ = t;
            externalInfluence_multiplierState_.i_ = agent->get_id();
            externalInfluence_multiplierState_.mu_v_.resize(Nxi*Nhor, 0.0);

            externalInfluence_penaltyState_.t_ = t;
            externalInfluence_penaltyState_.i_ = agent->get_id();
            externalInfluence_penaltyState_.rho_v_.resize(Nxi*Nhor, optimization_info.ADMM_PenaltyInit_);

            neighbors_localCopies_.t_ = t;
            neighbors_localCopies_.i_ = agent->get_id();
            neighbors_localCopies_.u_.resize(Nui*Nhor, 0.0);
            neighbors_localCopies_.v_.resize(Nxi*Nhor, 0.0);
            neighbors_localCopies_.lambda_.resize(Nxi * Nhor, 0.0);

            neighbors_coupled_multiplierState_.t_ = t;
            neighbors_coupled_multiplierState_.i_ = agent->get_id();
            neighbors_coupled_multiplierState_.mu_u_.resize(Nui*Nhor, 0.0);
            neighbors_coupled_multiplierState_.mu_v_.resize(Nxi*Nhor, 0.0);

            neighbors_coupled_penaltyState_.t_ = t;
            neighbors_coupled_penaltyState_.i_ = agent->get_id();
            neighbors_coupled_penaltyState_.rho_u_.resize(Nui*Nhor, optimization_info.ADMM_PenaltyInit_);
            neighbors_coupled_penaltyState_.rho_v_.resize(Nxi*Nhor, optimization_info.ADMM_PenaltyInit_);

            neighbors_couplingState_.t_ = t;
            neighbors_couplingState_.i_ = this->id_;
            neighbors_couplingState_.z_u_.resize(Nuj*Nhor, 0.0);

            previous_neighbors_couplingState_.t_ = t;
            previous_neighbors_couplingState_.i_ = this->id_;
            previous_neighbors_couplingState_.z_u_.resize(Nuj*Nhor, 0.0);

            neighbors_externalInfluence_couplingState_.t_ = t;
            neighbors_externalInfluence_couplingState_.i_ = this->get_id();
            neighbors_externalInfluence_couplingState_.z_v_.resize(Nxj*Nhor, 0.0);

            previous_neighbors_externalInfluence_couplingState_.t_ = t;
            previous_neighbors_externalInfluence_couplingState_.i_ = this->get_id();
            previous_neighbors_externalInfluence_couplingState_.z_v_.resize(Nxj*Nhor, 0.0);

            previous_externalInfluence_couplingState_.t_ = t;
            previous_externalInfluence_couplingState_.i_ = this->get_id();
            previous_externalInfluence_couplingState_.z_v_.resize(Nxj*Nhor, 0.0);

            neighbors_externalInfluence_multiplierState_.t_ = t;
            neighbors_externalInfluence_multiplierState_.i_ = this->get_id();
            neighbors_externalInfluence_multiplierState_.mu_v_.resize(Nxj*Nhor, 0.0);

            neighbors_externalInfluence_penaltyState_.t_ = t;
            neighbors_externalInfluence_penaltyState_.i_ = this->get_id();
            neighbors_externalInfluence_penaltyState_.rho_v_.resize(Nxj*Nhor, optimization_info.ADMM_PenaltyInit_);
        }

        if( is_approximatingCost() )
        {
            neighbors_desiredAgentState_.t_ = t;
            neighbors_desiredAgentState_.i_ = this->get_id();
            neighbors_desiredAgentState_.x_.resize(Nxj*Nhor, 0.0);
            neighbors_desiredAgentState_.u_.resize(Nuj*Nhor, 0.0);
            neighbors_desiredAgentState_.lambda_.resize(Nxj* Nhor, 0.0);
        }
       
        if( is_approximating() )
            approximate_neighbor_.reset(new grampcd::ApproximateNeighbor(agent, this));

        // for sensitivity-based algorithm
        neighbors_agentState_.t_ = t;
        neighbors_agentState_.i_ = agent->get_id();
        neighbors_agentState_.x_.resize(Nxj* Nhor, 0.0);
        neighbors_agentState_.u_.resize(Nuj* Nhor, 0.0);
        neighbors_agentState_.lambda_.resize(Nxj* Nhor, 0.0);

        if (is_sendingNeighbor())
        {
            const auto Ngij = get_couplingModel()->get_Ngij();
            const auto Nhij = get_couplingModel()->get_Nhij();

            coupled_constraintState_.t_ = t;
            coupled_constraintState_.i_ = agent->get_id();
            coupled_constraintState_.mu_g_.resize(Ngij* Nhor, 0.0);
            coupled_constraintState_.mu_h_.resize(Nhij* Nhor, 0.0);
            coupled_constraintState_.c_g_.resize(Ngij* Nhor, 0.0);
            coupled_constraintState_.c_h_.resize(Nhij* Nhor, 0.0);

        }
        if (is_receivingNeighbor())
        {
            const auto Ngji = get_copied_couplingModel()->get_Ngij();
            const auto Nhji = get_copied_couplingModel()->get_Nhij();

            sensiState_.t_ = t;
            sensiState_.i_ = agent->get_id();
            sensiState_.psi_x_.resize(Nxi * Nhor, 0.0);
            sensiState_.psi_u_.resize(Nui * Nhor, 0.0);
            sensiState_.psi_V_.resize(Nxi, 0.0);

            if (optimization_info.SENSI_higherOrder_)
            {
                sensiState_.psi_xx_.resize(Nxi * Nxi * Nhor, 0.0);
                sensiState_.psi_uu_.resize(Nui * Nui * Nhor, 0.0);
                sensiState_.psi_xu_.resize(Nxi * Nui * Nhor, 0.0);
                sensiState_.psi_VV_.resize(Nxi * Nxi, 0.0);
            }

            neighbors_coupled_constraintState_.t_ = t;
            neighbors_coupled_constraintState_.i_ = agent->get_id();
            neighbors_coupled_constraintState_.mu_g_.resize(Ngji* Nhor, 0.0);
            neighbors_coupled_constraintState_.mu_h_.resize(Nhji* Nhor, 0.0);
            neighbors_coupled_constraintState_.c_g_.resize(Ngji* Nhor, 0.0);
            neighbors_coupled_constraintState_.c_h_.resize(Nhji* Nhor, 0.0);

        }
    }

    void Neighbor::reset_neighborStates()
    {
        const auto& t = local_copies_.t_;

        const unsigned int i = agent_id_;
        const unsigned int j = get_id();

        resetState(local_copies_, j, t);
        resetState(coupled_multiplierState_, j, t);
        resetState(coupled_penaltyState_, j, t);

        resetState(neighbors_couplingState_, j, t);
        resetState(neighbors_localCopies_, i, t);
        resetState(neighbors_coupled_multiplierState_, i, t);
        resetState(neighbors_coupled_penaltyState_, i, t);

        resetState(externalInfluence_agentState_, j, t);
        resetState(externalInfluence_couplingState_, j, t);
        resetState(externalInfluence_multiplierState_, j, t);
        resetState(externalInfluence_penaltyState_, j, t);

        resetState(neighbors_externalInfluence_penaltyState_, j, t);
        resetState(neighbors_externalInfluence_couplingState_, j, t);
        resetState(neighbors_externalInfluence_multiplierState_, j, t);

        resetState(neighbors_agentState_, j, t);
        resetState(sensiState_, j, t);
        resetState(coupled_constraintState_, j, t);
        resetState(neighbors_coupled_constraintState_, j, t);

        resetState(previous_externalInfluence_couplingState_, j, t);
        resetState(previous_externalInfluence_multiplierState_, j, t);
        resetState(previous_couplingState_, j, t);
        resetState(previous_multiplierState_, j, t);

        resetState(neighbors_desiredAgentState_, j, t);
    }

    const bool Neighbor::is_sendingNeighbor() const
    {
        return sending_neighbor_;
    }

    const bool Neighbor::is_receivingNeighbor() const
    {
        return receiving_neighbor_;
    }

    void Neighbor::defineAs_sendingNeighbor(const bool isSending)
    {
        sending_neighbor_ = isSending;
    }

    void Neighbor::defineAs_receivingNeighbor(const bool isReceiving)
    {
        receiving_neighbor_ = isReceiving;
    }

    void Neighbor::shift_states(const typeRNum dt, const typeRNum t0)
    {
        shiftState(local_copies_, dt, t0);
        shiftState(coupled_multiplierState_, dt, t0);
        shiftState(coupled_penaltyState_, dt, t0);

        shiftState(neighbors_localCopies_, dt, t0);
        shiftState(neighbors_couplingState_, dt, t0);
        shiftState(neighbors_coupled_multiplierState_, dt, t0);
        shiftState(neighbors_coupled_penaltyState_, dt, t0);

        shiftState(externalInfluence_agentState_, dt, t0);
        shiftState(externalInfluence_couplingState_, dt, t0);
        shiftState(externalInfluence_multiplierState_, dt, t0);
        shiftState(externalInfluence_penaltyState_, dt, t0);

        shiftState(neighbors_externalInfluence_couplingState_, dt, t0);
        shiftState(neighbors_externalInfluence_multiplierState_, dt, t0);
        shiftState(neighbors_externalInfluence_penaltyState_, dt, t0);

        shiftState(neighbors_agentState_, dt, t0);
        shiftState(sensiState_, dt, t0);
        shiftState(coupled_constraintState_, dt, t0);
        shiftState(neighbors_coupled_constraintState_, dt, t0);

        shiftState(previous_couplingState_, dt, t0);
        shiftState(previous_multiplierState_, dt, t0);

        shiftState(previous_externalInfluence_couplingState_, dt, t0);
        shiftState(previous_externalInfluence_multiplierState_, dt, t0);
        shiftState(previous_neighbors_externalInfluence_couplingState_, dt, t0);
        shiftState(previous_neighbors_couplingState_, dt, t0);
    }

    const unsigned int Neighbor::get_Nui() const
    {
        if( is_sendingNeighbor() )
            return couplingModel_->get_Nui();
        else
            return copied_couplingModel_->get_Nuj();
    }

    const unsigned int Neighbor::get_Nxi() const
    {
        if( is_sendingNeighbor() )
            return couplingModel_->get_Nxi();
        else
            return copied_couplingModel_->get_Nxj();
    }

    const unsigned int Neighbor::get_Nuj() const
    {
        if( is_sendingNeighbor() )
            return couplingModel_->get_Nuj();
        else
            return copied_couplingModel_->get_Nui();
    }

    const unsigned int Neighbor::get_Nxj() const
    {
        if( is_sendingNeighbor() )
            return couplingModel_->get_Nxj();
        else
            return copied_couplingModel_->get_Nxi();
    }

    void Neighbor::increase_delays(const AlgStep& step)
    {
        switch (step)
        {
        case(AlgStep::ADMM_UPDATE_AGENT_STATE):
            ++delay_agentState_;
            break;

        case(AlgStep::ADMM_UPDATE_COUPLING_STATE):
            ++delay_couplingState_;
            break;

        case(AlgStep::ADMM_UPDATE_MULTIPLIER_STATE):
            ++delay_multiplierState_;
            break;

        case (AlgStep::SENSI_UPDATE_AGENT_STATE):
            ++delay_sensiAgentState_;
            break;

        default:
            break;
        }
    }

    void Neighbor::reset_delays(const AlgStep& step)
    {
        switch (step)
        {
        case(AlgStep::ADMM_UPDATE_AGENT_STATE):
            delay_agentState_ = 0;
            break;

        case(AlgStep::ADMM_UPDATE_COUPLING_STATE):
            delay_couplingState_ = 0;
            break;

        case(AlgStep::ADMM_UPDATE_MULTIPLIER_STATE):
            delay_multiplierState_ = 0;
            break;

        case (AlgStep::SENSI_UPDATE_AGENT_STATE):
            delay_sensiAgentState_ = 0;
            break;

        default:
            break;
        }
    }

    int Neighbor::get_delays(const AlgStep& step)
    {
        switch (step)
        {
        case(AlgStep::ADMM_UPDATE_AGENT_STATE):
            return delay_agentState_;
            break;

        case(AlgStep::ADMM_UPDATE_COUPLING_STATE):
            return delay_couplingState_;
            break;

        case(AlgStep::ADMM_UPDATE_MULTIPLIER_STATE):
            return delay_multiplierState_;
            break;

        case (AlgStep::SENSI_UPDATE_AGENT_STATE):
            return delay_sensiAgentState_;
            break;

        default:
            return 0;
            break;
        }
    }

    void Neighbor::initialize_delays()
    {
        delay_agentState_ = 65535;
        delay_couplingState_ = 65535;
        delay_multiplierState_ = 0;
        delay_sensiAgentState_ = 0;   
    }

}
