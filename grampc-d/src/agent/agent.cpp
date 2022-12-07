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

#include "grampcd/util/data_conversion.hpp"
#include "grampcd/util/logging.hpp"

#include "grampcd/model/model_factory.hpp"
#include "grampcd/model/agent_model.hpp"
#include "grampcd/model/coupling_model.hpp"

#include "grampcd/comm/communication_interface.hpp"

#include "grampcd/optim/optim_util.hpp"
#include "grampcd/optim/solver_local.hpp"
#include "grampcd/optim/solver_local_admm.hpp"
#include "grampcd/optim/solver_local_sensi.hpp"
#include "grampcd/optim/solution.hpp"
#include "grampcd/agent/step_selector.hpp"
#include "grampcd/agent/sync_step_selector.hpp"
#include "grampcd/agent/async_step_selector.hpp"



namespace grampcd
{

    Agent::Agent(const CommunicationInterfacePtr& communication_interface,
                 const ModelFactoryPtr& model_factory,
                 const AgentInfo& agent_info,
         const LoggingPtr& log)
        :   
        communication_interface_(communication_interface),
        model_factory_(model_factory),
        model_(model_factory->create_agentModel(agent_info)),
        info_(agent_info),
        local_solver_(nullptr),
        solution_(new Solution()),
        log_(log),
        step_selector_(nullptr)
    {
        x_des_ = std::vector<typeRNum>(model_->get_Nxi(), 0.0);
        u_des_ = std::vector<typeRNum>(model_->get_Nui(), 0.0);
    }

    const int Agent::get_id() const
    {
        return info_.id_;
    }

    const AgentInfo& Agent::get_agentInfo() const
    {
        return info_;
    }

    const AgentModelPtr& Agent::get_agentModel() const
    {
        return model_;
    }

    const std::shared_ptr< std::map<int, CouplingModelPtr> > Agent::get_couplingModels() const
    {
        std::shared_ptr< std::map<int, CouplingModelPtr> > couplings(new std::map<int, CouplingModelPtr>);

        // consider coupling to each sending neighbor
        for( const auto& neighbor : sending_neighbors_ )
            couplings->insert(std::make_pair(neighbor->get_id(), neighbor->get_couplingModel()));

        return couplings;
    }

    const std::vector<NeighborPtr>& Agent::get_neighbors() const
    {
        return neighbors_;
    }

    const bool Agent::add_neighbor(const grampcd::CouplingInfo &coupling_info, const grampcd::AgentInfo &neighbor_info)
    {
        // either add sending or receiving neighbor
        if( coupling_info.agent_id_ == get_id() )
            return add_sendingNeighbor(coupling_info, neighbor_info);
        else
            return add_receivingNeighbor(coupling_info, neighbor_info);
    }

    const bool Agent::add_sendingNeighbor(const grampcd::CouplingInfo &coupling_info, const grampcd::AgentInfo &neighbor_info)
    {
        // check if neighbor already exists
        NeighborPtr neighbor = nullptr;
        for(const auto& search_neighbor : neighbors_)
        {
            if( search_neighbor->get_id() == coupling_info.neighbor_id_ )
            {
                if( search_neighbor->is_sendingNeighbor() )
                {
                    log_->print(DebugType::Error) << "[Agent::add_sendingNeighbor] Agent " << get_id() << ": "
                        << "Failed to add neighbor, as agent " << coupling_info.agent_id_
                        << " already has sending neighbor " << coupling_info.neighbor_id_ << "." << std::endl;
                    return false;
                }
                else
                {
                    neighbor = search_neighbor;
                    break;
                }
            }
        }

        // if neighbor does not exist, create new neighbor
        if( neighbor == nullptr )
        {
            const auto agent_model = model_factory_->create_agentModel(neighbor_info);
            neighbor = NeighborPtr(new Neighbor(coupling_info.neighbor_id_, agent_model, log_));
            neighbors_.push_back(neighbor);
        }

        // create coupling model
        const auto coupling_model = model_factory_->create_couplingModel(coupling_info);
        neighbor->set_couplingModel( coupling_model );

        // set flag is_sendingNeighbor and add neighbor to the list sending_neighbors_
        neighbor->defineAs_sendingNeighbor( true );
        sending_neighbors_.push_back(neighbor);

        // check if coupling is correct initialized
        if( (model_->get_Nxi() != neighbor->get_couplingModel()->get_Nxi()) ||
            (model_->get_Nui() != neighbor->get_couplingModel()->get_Nui()) )
        {
            log_->print(DebugType::Error) << "[Agent::add_sendingNeighbor] Agent " << get_id() << ": "
                << "Invalid coupling model between agent " << get_id() << " and neighbor " << neighbor->get_id() 
                << ". " << std::endl;

            //delete neighbor
            neighbors_.pop_back();
            sending_neighbors_.erase( sending_neighbors_.begin() + sending_neighbors_.size() - 1 );

            return false;
        }

        // initialize neighbor
        neighbor->reset_neighborStates();
        neighbor->initialize_neighbor(agentState_.t_, optimizationInfo_, this);

        return true;
    }

    const bool Agent::add_receivingNeighbor(const grampcd::CouplingInfo &coupling_info, const grampcd::AgentInfo &neighbor_info)
    {
        // check if neighbor already exists
        NeighborPtr neighbor = nullptr;
        for( NeighborPtr search_neighbor : neighbors_ )
        {
            if( search_neighbor->get_id() == coupling_info.agent_id_ )
            {
                if( search_neighbor->is_receivingNeighbor() )
                {
                    log_->print(DebugType::Error) << "[Agent::add_receivingNeighbor] Agent " << get_id() << ": "
                        << "Failed to add neighbor, as agent " << coupling_info.neighbor_id_
                        << " already has receiving neighbor " << coupling_info.agent_id_
                        << std::endl;

                    return false;
                }
                else
                {
                    neighbor = search_neighbor;
                    break;
                }
            }
        }

        // if neighbor does not exist, create new neighbor
        if( neighbor == nullptr )
        {
            const auto agent_model = model_factory_->create_agentModel(neighbor_info);
            neighbor = NeighborPtr(new Neighbor(coupling_info.agent_id_, agent_model, log_));
            neighbors_.push_back(neighbor);
        }

        // create coupling model
        const CouplingModelPtr coupling_model = model_factory_->create_couplingModel(coupling_info);
        neighbor->set_copied_couplingModel( coupling_model );
        neighbor->defineAs_receivingNeighbor( true );

        receiving_neighbors_.push_back(neighbor);

        // initialize neighbor
        neighbor->reset_neighborStates();
        neighbor->initialize_neighbor(agentState_.t_, optimizationInfo_, this);

        return true;
    }

    void Agent::remove_neighbor(const grampcd::CouplingInfo &coupling_info)
    {
        NeighborPtr neighbor = nullptr;
        // check if neighbor exists
        for( const auto& searchNeighbor : neighbors_ )
        {
            if( searchNeighbor->get_id() == coupling_info.neighbor_id_ ||
                   searchNeighbor->get_id() == coupling_info.agent_id_ )
            {
                neighbor = searchNeighbor;
                break;
            }
        }

        // check if neighbor is found
        if (neighbor == nullptr)
        {
            log_->print(DebugType::Error) << "[Agent::remove_neighbor] Agent " << get_id() << ": "
                << "Neighbor not found." << std::endl;
            return;
        }

        // check if receiving or sending neighbor should be deleted
        if (coupling_info.agent_id_ == neighbor->get_id())
        {
            neighbor->defineAs_receivingNeighbor(false);
            DataConversion::erase_element_from_vector(receiving_neighbors_, neighbor);
        }
        else
        {
            neighbor->defineAs_sendingNeighbor(false);
            DataConversion::erase_element_from_vector(sending_neighbors_, neighbor);
        }

        // check if neighbor is still a neighbor
        const bool is_neighbor = DataConversion::is_element_in_vector(sending_neighbors_, neighbor) ||
            DataConversion::is_element_in_vector(receiving_neighbors_, neighbor);

        // remove neighbor from list neighbors
        if (!is_neighbor)
            DataConversion::erase_element_from_vector(neighbors_, neighbor);
    }

    /*************************************************************************
     state functions
     *************************************************************************/

    void Agent::initialize(const OptimizationInfo& optimization_info)
    {
        optimizationInfo_ = optimization_info;
        initial_penalty_ = optimizationInfo_.ADMM_PenaltyInit_;

        is_approximatingCost_ = optimizationInfo_.APPROX_ApproximateCost_;
        is_approximatingDynamics_ = optimization_info.APPROX_ApproximateDynamics_;
        is_approximatingConstraints_ = optimization_info.APPROX_ApproximateConstraints_;

        is_approximating_ = is_approximatingCost_ || is_approximatingDynamics_ || is_approximatingConstraints_;

        const auto Nhor = optimizationInfo_.COMMON_Nhor_;
        const typeRNum Thor = optimizationInfo_.COMMON_Thor_;
        const auto Nxi = get_agentModel()->get_Nxi();
        const auto Nui = get_agentModel()->get_Nui();
        const auto Nhi = get_agentModel()->get_Nhi();
        const int id = info_.id_;

        // create uniform time discretization
        std::vector<typeRNum> t = std::vector<typeRNum>(Nhor, 0.0);
        for(unsigned int i = 0; i < Nhor; ++i)
            t[i] = i * Thor / (Nhor - 1.0);

        // reset states
        resetState(agentState_, id, t);
        resetState(couplingState_, id, t);
        resetState(multiplierState_, id, t);
        resetState(desired_agentState_, id, t);

        resetState(previous_couplingState_, id, t);
        resetState(previous_multiplierState_, id, t);
        resetState(previous_agentState_, id, t);

        resetState(penaltyState_, id, t);

        // resize states
        agentState_.x_.resize( Nhor * Nxi, 0.0 );
        agentState_.u_.resize( Nhor * Nui, 0.0 );
        agentState_.lambda_.resize(Nhor * Nxi, 0.0);
        couplingState_.z_u_.resize( Nhor * Nui, 0.0 );
        multiplierState_.mu_u_.resize( Nhor * Nui, 0.0 );
	    penaltyState_.rho_u_.resize(Nhor * Nui, initial_penalty_);
        desired_agentState_.x_.resize(Nhor * Nxi, 0.0);
        desired_agentState_.u_.resize(Nhor * Nui, 0.0);
        desired_agentState_.lambda_.resize(Nhor * Nxi, 0.0);

        if(!is_approximatingDynamics_)
        {
            couplingState_.z_x_.resize( Nhor * Nxi, 0.0 );
            multiplierState_.mu_x_.resize( Nhor * Nxi, 0.0 );
            penaltyState_.rho_x_.resize( Nhor * Nxi, initial_penalty_ );
        }

        // initialize states
        for(unsigned int i = 0; i < Nhor; ++i)
        {
            for(unsigned int j = 0; j < model_->get_Nxi(); ++j)
            {
			    agentState_.x_[i * model_->get_Nxi() + j] = x_init_[j];
                desired_agentState_.x_[i * model_->get_Nxi() + j] = x_des_[j];
            }
        }

        // initialize controls
        for(unsigned int i = 0; i < Nhor; ++i)
        {
            for(unsigned int j = 0; j < model_->get_Nui(); ++j)
            {
			    agentState_.u_[i * model_->get_Nui() + j] = u_init_[j];
			    desired_agentState_.u_[i * model_->get_Nui() + j] = u_des_[j];
            }
        }

        // initialize adjoint states (to initial state)
        std::vector<typeRNum> adj_init;
        adj_init.resize(model_->get_Nxi(), 0.0);
        model_->dVdx(&adj_init[0], agentState_.t_[Nhor - 1], &x_init_[0], &x_des_[0]);
        for (unsigned int i = 0; i < Nhor; ++i)
        {
            for (unsigned int j = 0; j < model_->get_Nxi(); ++j)
            {
                agentState_.lambda_[i * model_->get_Nxi() + j] = adj_init[j];
            }
        }

        // initialize coupling states
        couplingState_.z_u_ = agentState_.u_;
        if( !is_approximating_ )
            couplingState_.z_x_ = agentState_.x_;

        // initialize previous states
        previous_couplingState_ = couplingState_;
        previous_multiplierState_ = multiplierState_;
        previous_agentState_ = agentState_; 

        //initialize neighbors
        for( const auto& neighbor : neighbors_ )
        {
            neighbor->reset_neighborStates();
            neighbor->initialize_neighbor(agentState_.t_, optimizationInfo_, this);
        }
    }

    void Agent::shift_states(typeRNum dt, typeRNum t0)
    {
        // shift own states
        shiftState(agentState_, dt, t0);
        shiftState(couplingState_, dt, t0);
        shiftState(multiplierState_, dt, t0);
        shiftState(penaltyState_, dt, t0);
        shiftState(previous_couplingState_, dt, t0);
        shiftState(previous_multiplierState_, dt, t0);
        shiftState(previous_agentState_, dt, t0);

        // shift states of neighbor
        for( const auto& neighbor : neighbors_ )
            neighbor->shift_states(dt, t0);
    }

    const AgentState& Agent::get_desiredAgentState() const
    {
        return desired_agentState_;
    }

    void Agent::set_desiredAgentState(const AgentState& state)
    {
        if( compare_stateDimensions( desired_agentState_, state ) )
            desired_agentState_ = agentState_;
        else
            log_->print(DebugType::Error) << "[Agent::set_desiredAgentState] Agent " << get_id() << ": "
            << "Failed to set desired agent state, as dimensions don't fit." << std::endl;
    }

    void Agent::set_desiredAgentState(const std::vector<typeRNum>& x_des)
    {
        x_des_ = x_des;
        u_des_ = std::vector<typeRNum>(model_->get_Nui(), 0.0);
    }

    void Agent::set_desiredAgentState(const std::vector<typeRNum>& x_des, const std::vector<typeRNum>& u_des)
    {
	    x_des_ = x_des;
        u_des_ = u_des;
    }

    const AgentState& Agent::get_agentState() const
    {
        return agentState_;
    }

    const AgentState& Agent::get_previous_agentState() const
    {
        return previous_agentState_;
    }

    void Agent::set_agentState(const AgentState& state)
    {
        if(!compare_stateDimensions( agentState_, state ))
        {
            log_->print(DebugType::Error) << "[Agent::set_agentState] Agent " << get_id() << ": "
                << "Failed to set agent state, as dimensions don't fit." << std::endl;

            return;
	    }

        // set agent state
        previous_agentState_ = agentState_;
	    agentState_ = state;

        // evaluate predicted cost
	    std::vector< typeRNum > predicted_cost(optimizationInfo_.COMMON_Nhor_, 0.0);
        for (unsigned int i = 0; i < predicted_cost.size(); ++i)
        {
            //check if the number of inputs of the agent is greater than 0
            ctypeRNum* u = agentState_.u_.size() > 0 ? &agentState_.u_[i] : nullptr;
            model_->lfct(&predicted_cost[i], agentState_.t_[i], &agentState_.x_[i], u, &desired_agentState_.x_[i]);
        }

        // update solution
	    solution_->update_predicted_state(state, predicted_cost);
    }

    const CouplingState& Agent::get_couplingState() const
    {
        return couplingState_;
    }

    const CouplingState& Agent::get_previous_couplingState() const
    {
        return previous_couplingState_;
    }

    void Agent::set_couplingState(const CouplingState& state)
    {
        if(compare_stateDimensions( couplingState_, state ))
        {
            previous_couplingState_ = couplingState_;
            couplingState_ = state;
        }
        else
            log_->print(DebugType::Error) << "[Agent::set_couplingState] Agent " << get_id() << ": "
            << "Failed to set coupling state, as dimensions don't fit." << std::endl;
    }

    const MultiplierState& Agent::get_multiplierState() const
    {
        return multiplierState_;
    }

    const MultiplierState& Agent::get_previous_multiplierState() const
    {
        return previous_multiplierState_;
    }

    void Agent::set_multiplierState(const MultiplierState& state)
    {
        if(compare_stateDimensions( state, multiplierState_ ))
        {
            previous_multiplierState_ = multiplierState_;
            multiplierState_ = state;
        }
        else
            log_->print(DebugType::Error) << "[Agent::set_multiplierState] Agent " << get_id() << ": "
            << "Failed to set multiplier state, as dimensions don't fit." << std::endl;
    }

    const PenaltyState& Agent::get_penaltyState() const
    {
        return penaltyState_;
    }

    void Agent::set_penaltyState(const PenaltyState &penalty)
    {
        if(compare_stateDimensions( penaltyState_, penalty ))
            penaltyState_ = penalty;
        else
            log_->print(DebugType::Error) << "[Agent::get_penaltyState] Agent " << get_id() << ": "
            << "Failed to set penalty state, as dimensions don't fit." << std::endl;
    }

    void Agent::fromCommunication_registered_coupling(const CouplingInfo& coupling_info, const AgentInfo& neighbor_info)
    {
        add_neighbor(coupling_info, neighbor_info);

        // reinitialize after adding neighbor
        fromCommunication_configured_optimization(optimizationInfo_);

        log_->print(DebugType::Message) << "[Agent::fromCommunication_registered_coupling] Agent " << get_id() << ": "
            << "Registered coupling between agent "
            << coupling_info.agent_id_ << " (receiving) and "
            << coupling_info.neighbor_id_ << " (sending)." << std::endl;
    }

    void Agent::fromCommunication_deregistered_coupling(const CouplingInfo& coupling_info)
    {
        remove_neighbor(coupling_info);

        log_->print(DebugType::Message) << "[Agent::fromCommunication_deregistered_coupling] Agent " << get_id() << ": "
            << " Deregistered coupling with agent "
            << coupling_info.agent_id_ << " (receiving) and "
            << coupling_info.neighbor_id_ << " (sending)." << std::endl;
    }

    void Agent::fromCommunication_received_localCopies(const AgentState& state, int from)
    {
        if (state.i_ != get_id())
        {
			log_->print(DebugType::Warning) << "[Agent::fromCommunication_received_localCopies] Agent " << get_id() << ": "
				<< "Agent " << get_id() << " received state of agent " << state.i_ << " from agent "
				<< from << "." << std::endl;

            return;
        }
        
		const auto neighbor = DataConversion::get_element_from_vector(neighbors_, from);

        if (neighbor == nullptr)
        {
			log_->print(DebugType::Warning) << "[Agent::fromCommunication_received_localCopies] Agent " << get_id() << ": "
				<< "Could not find neighbor with id " << from << "." << std::endl;

            return;
        }
        
        neighbor->set_neighbors_localCopies(state);

        if (optimizationInfo_.ASYNC_Active_)
        {
            neighbor->reset_delays(AlgStep::ADMM_UPDATE_AGENT_STATE);
            step_selector_->execute_algStep(AlgStep::ADMM_UPDATE_COUPLING_STATE);    
        }
    }

    void Agent::fromCommunication_received_agentState(const AgentState& state, const ConstraintState& constr_state, int from)
    {
        const auto neighbor = DataConversion::get_element_from_vector(neighbors_, from);

        if (state.i_ != neighbor->get_id())
        {
            log_->print(DebugType::Warning) << "[Agent::fromCommunication_received_agentState] Agent " << get_id() << ": "
                << "Agent " << get_id() << " received state of agent " << state.i_ << " from agent "
                << from << "." << std::endl;

            return;
        }

        if (neighbor == nullptr)
        {
            log_->print(DebugType::Warning) << "[Agent::fromCommunication_received_agentState] Agent " << get_id() << ": "
                << "Could not find neighbor with id " << from << "." << std::endl;

            return;
        }

        // set neighbors state
        neighbor->set_neighbors_agentState(state);
        // only set constraint state if neighbor is receiving neighbor
        if (neighbor->is_receivingNeighbor())
            neighbor->set_neighbors_coupled_constraintState(constr_state);

        if (optimizationInfo_.ASYNC_Active_)
        {
            neighbor->reset_delays(AlgStep::SENSI_UPDATE_AGENT_STATE);
            step_selector_->execute_algStep(AlgStep::SENSI_UPDATE_AGENT_STATE);
        }
    }

    void Agent::fromCommunication_received_desiredAgentState(const AgentState& state, int from)
    {
        if (state.i_ != from)
        {
			log_->print(DebugType::Warning) << "[Agent::fromCommunication_received_desiredAgentState] Agent " << get_id() << ": "
				<< "Agent " << get_id() << " received state of agent " << state.i_ << " from agent " << from << "." << std::endl;

            return;
        }
            
		const auto neighbor = DataConversion::get_element_from_vector(neighbors_, from);

        if (neighbor == nullptr)
        {
			log_->print(DebugType::Warning) << "[Agent::fromCommunication_received_desiredAgentState] Agent " << get_id() << ": "
				<< "Could not find neighbor with id " << from << "." << std::endl;

            return;
        }
			
        neighbor->set_neighbors_desiredAgentState(state);
    }

    void Agent::fromCommunication_received_couplingState(const CouplingState& coupling, int from)
    {
        // received neighbors coupling state
        if (coupling.i_ != from)
        {
			log_->print(DebugType::Warning) << "[Agent::fromCommunication_received_couplingState] Agent " << get_id() << ": "
				<< "Agent " << get_id() << " received coupling state of agent " << coupling.i_
				<< " from agent " << from << "." << std::endl;

            return;
        }
        
		const auto neighbor = DataConversion::get_element_from_vector(neighbors_, from);

		if (neighbor == nullptr)
		{
			log_->print(DebugType::Warning) << "[Agent::fromCommunication_received_couplingState] Agent " << get_id() << ": "
				<< "Could not find agent with id " << from << "." << std::endl;
			return;
		}

        if (!compare_stateDimensions(neighbor->get_neighbors_couplingState(), coupling))
        {
			log_->print(DebugType::Error) << "[Agent::fromCommunication_received_couplingState] Agent " << get_id() << ": "
				<< "Failed to receive coupling, as dimensions don't fit." << std::endl;

            return;
        }
			
        neighbor->set_neighbors_couplingState(coupling);

        // asynchronous execution
        if (optimizationInfo_.ASYNC_Active_)
        {
            neighbor->reset_delays(AlgStep::ADMM_UPDATE_COUPLING_STATE);
            step_selector_->execute_algStep(AlgStep::ADMM_UPDATE_MULTIPLIER_STATE);
        }
    }

    void Agent::fromCommunication_received_couplingState(const CouplingState& coupling, const CouplingState& ext_influence_coupling, int from)
    {
        // received neighbors coupling state
        if (coupling.i_ != from)
        {
			log_->print(DebugType::Warning) << "[Agent::fromCommunication_received_couplingState] Agent " << get_id() << ": "
				<< "Agent " << get_id() << " received coupling state of agent "
				<< coupling.i_ << " from agent " << from << "." << std::endl;

            return;
        }

		const auto neighbor = DataConversion::get_element_from_vector(neighbors_, from);

		if (neighbor == nullptr)
		{
			log_->print(DebugType::Warning) << "[Agent::fromCommunication_received_couplingState] Agent " << get_id() << ": "
				<< "Could not find agent with id " << from << "." << std::endl;
			return;
		}

        if (!compare_stateDimensions(neighbor->get_neighbors_couplingState(), coupling) ||
            !compare_stateDimensions(neighbor->get_neighbors_externalInfluence_couplingState(), ext_influence_coupling))
        {
			log_->print(DebugType::Error) << "[Agent::fromCommunication_received_couplingState] Agent " << get_id() << ": "
				<< "Failed to receive coupling, as dimensions don't fit." << std::endl;

            return;
        }		
			
        neighbor->set_neighbors_externalInfluence_couplingState(ext_influence_coupling);			
		neighbor->set_neighbors_couplingState(coupling);

        // asynchronous execution
        if (optimizationInfo_.ASYNC_Active_)
        {
            neighbor->reset_delays(AlgStep::ADMM_UPDATE_COUPLING_STATE);
            step_selector_->execute_algStep(AlgStep::ADMM_UPDATE_MULTIPLIER_STATE);
        }
    }

    void Agent::fromCommunication_received_multiplierState(const MultiplierState& multiplier, PenaltyState penalty, int from)
    {
        // received neighbors coupled multiplier state
        if (multiplier.i_ != get_id())
        {
			log_->print(DebugType::Warning) << "[Agent::fromCommunication_received_multiplierState] Agent " << get_id() << ": "
				<< "Agent " << get_id() << " received multiplier state of agent " << multiplier.i_
				<< " from agent " << from << "." << std::endl;

            return;
        }
        
		const auto neighbor = DataConversion::get_element_from_vector(neighbors_, from);

		if (neighbor == nullptr)
		{
			log_->print(DebugType::Warning) << "[Agent::fromCommunication_received_multiplierState] Agent " << get_id() << ": "
				<< "Could not find agent with id " << from << "." << std::endl;
			return;
		}

        if (!compare_stateDimensions(neighbor->get_neighbors_coupled_multiplierState(), multiplier) ||
            !compare_stateDimensions(neighbor->get_neighbors_coupled_penaltyState(), penalty))
        {
			log_->print(DebugType::Error) << "[Agent::fromCommunication_received_multiplierState] Agent " << get_id() << ": "
				<< "Failed to receive multiplier state, as dimensions don't fit." << std::endl;

            return;
        }

		neighbor->set_neighbors_coupled_multiplierState(multiplier);
		neighbor->set_neighbors_coupled_penaltyState(penalty);

        // asynchronous execution
        if (optimizationInfo_.ASYNC_Active_)
        {
            neighbor->reset_delays(AlgStep::ADMM_UPDATE_MULTIPLIER_STATE);
            step_selector_->execute_algStep(AlgStep::ADMM_UPDATE_AGENT_STATE);
        }
    }

    void Agent::fromCommunication_configured_optimization(const OptimizationInfo& info)
    {
        initialize(info);

        // create local solver
        if (optimizationInfo_.COMMON_Solver_ == "ADMM")
            local_solver_.reset(new SolverLocalADMM(this, info, log_, communication_interface_));
        else if (optimizationInfo_.COMMON_Solver_ == "Sensi")
            local_solver_.reset(new SolverLocalSensi(this, info, log_, communication_interface_));
        else
            log_->print(DebugType::Error) << "[Agent::fromCommunication_configured_optimization] Agent" << get_id() << ": "
            << "is not initialized with a proper solver" << std::endl;

        // create step selector
        if (optimizationInfo_.ASYNC_Active_)
            step_selector_.reset(new AsyncStepSelector(local_solver_, this));
        else
            step_selector_.reset(new SyncStepSelector(local_solver_, this));

    }

    const OptimizationInfo& Agent::get_optimizationInfo() const
    {
        return optimizationInfo_;
    }

    void Agent::fromCommunication_received_numberOfNeighbors(const unsigned int number, const int from)
    {
        const auto neighbor = DataConversion::get_element_from_vector(neighbors_, from);

        if (neighbor == nullptr)
        {
			log_->print(DebugType::Warning) << "[Agent::onNumerOfNeighborsReceived] Agent " << get_id() << ": "
				<< "Received number of neighbors from unknown neighbor " << from << "." << std::endl;

            return;
        }
        
        neighbor->set_numberOfNeighbors(number);
    }

    void Agent::fromCommunication_trigger_step(const AlgStep& step)
    {
       step_selector_->execute_algStep(step);
    }

    void Agent::fromCommunication_received_flagToStopAdmm(const bool flag)
    {
        // remember flag stop algorithm
        std::static_pointer_cast<AsyncStepSelector>(step_selector_)->set_flagStopAlg(flag);

        // send acknowledge flag
        communication_interface_->send_stoppedAlgFlag(flag,this->get_id());
    }

    void Agent::fromCommunication_received_flagStoppedAdmm(const bool flag, int from)
    {

        const auto neighbor = DataConversion::get_element_from_vector(neighbors_, from);

        if (neighbor == nullptr)
        {
            log_->print(DebugType::Warning) << "[Agent::fromCommunication_received_admmIter] Agent " << get_id() << ": "
                << "Received number of neighbors from unknown neighbor " << from << "." << std::endl;

            return;
        }
        neighbor->flag_StoppedAlg_ = true;

    }
    

    void Agent::set_neighbors_initial_states()
    {
        for (const auto& neighbor : get_neighbors())
        {
            if (optimizationInfo_.APPROX_ApproximateDynamics_)
            {
                auto state = neighbor->get_localCopies();
                const auto trueState = communication_interface_->get_agentState_from_agent(neighbor->get_id());

                if (trueState == nullptr)
                {
                    log_->print(DebugType::Warning) << "[Agent::set_neighbors_initial_states] Agent " << get_id() << ": "
                        << "Neighbor with id " << neighbor->get_id()
                        << " did not respond. Initialization failed." << std::endl;

                    return;
                }

                for (unsigned int j = 0; j < neighbor->get_Nxj(); ++j)
                    state.x_[j] = trueState->x_[j];

                neighbor->set_localCopies(state);
            }

            if (optimizationInfo_.APPROX_ApproximateCost_)
            {
                const auto desiredAgentState = communication_interface_->get_desiredAgentState_from_agent(neighbor->get_id());

                if (desiredAgentState == nullptr)
                {
                    log_->print(DebugType::Warning) << "[Agent::set_neighbors_initial_states] Agent " << get_id() << ": "
                        << "Neighbor with id " << neighbor->get_id()
                        << " did not respond. Initialization failed." << std::endl;

                    return;
                }

                neighbor->set_neighbors_desiredAgentState(*desiredAgentState);
            }
        }    
    }

    void Agent::set_updatedState(const std::vector<typeRNum>& new_state, const typeRNum dt, const typeRNum t0, const typeRNum cost)
    {
        // shift all states
        shift_states(dt, t0);    

        // set simulated state
        for( int k = 0; k < new_state.size(); ++k )
        {
            agentState_.x_[k] = new_state[k];
            x_init_[k] = new_state[k];
        }

        // update solution
        solution_->update_state(new_state, t0, get_id(), cost);
    }

    const std::vector<NeighborPtr>& Agent::get_sendingNeighbors() const
    {
        return sending_neighbors_;
    }

    const std::vector<NeighborPtr>& Agent::get_receivingNeighbors() const
    {
        return receiving_neighbors_;
    }

    const unsigned int Agent::get_Nxi() const
    {
        return model_->get_Nxi();
    }

    const unsigned int Agent::get_Nui() const
    {
        return model_->get_Nui();
    }

    const bool Agent::is_approximating() const
    {
        return is_approximating_;
    }

    const bool Agent::is_approximatingCost() const
    {
        return is_approximatingCost_;
    }

    const bool Agent::is_approximatingConstraints() const
    {
        return is_approximatingConstraints_;
    }

    const bool Agent::is_approximatingDynamics() const
    {
        return is_approximatingDynamics_;
    }

    const SolutionPtr& Agent::get_solution() const
    {
        return solution_;
    }

    void Agent::set_solution(const SolutionPtr& solution)
    {
        solution_ = solution;
    }

    void Agent::reset_solution()
    {
        solution_.reset(new Solution());
    }

    void Agent::set_initialState(const std::vector<typeRNum>& x_init)
    {
        if (x_init_.size() != x_init.size())
        {
            log_->print(DebugType::Error) << "[Agent::set_initialState]: "
                << "Dimension of new initial state differs from "
                << "dimension of old initial state." << std::endl;

            return;
        }

        x_init_ = x_init;
    }

    void Agent::set_initialState(const std::vector<typeRNum>& x_init, const std::vector<typeRNum>& u_init)
    {
        x_init_ = x_init;
        u_init_ = u_init;
    }

    const typeRNum Agent::get_predicted_cost() const
    {

	    // evaluate Nhor
	    const unsigned int Nhor = static_cast<unsigned int>(agentState_.t_.size());

	    const unsigned int Nxi = get_Nxi();
	    const unsigned int Nui = get_Nui();

	    // prepare cost
	    typeRNum cost = 0.0;

	    // consider horizon
        for (unsigned int i = 0; i < Nhor; ++i)
        {
            get_agentModel()->lfct(&cost, agentState_.t_[i], &agentState_.x_[i * Nxi], &agentState_.u_[i * Nui], &desired_agentState_.x_[0]);

            for (const auto& neighbor : get_sendingNeighbors())
            {
                const auto Nxj = neighbor->get_Nxj();
                const auto Nuj = neighbor->get_Nuj();
                if (optimizationInfo_.COMMON_Solver_ == "ADMM")
                {
                    const auto& local_copies = neighbor->get_localCopies();

                    neighbor->get_couplingModel()->lfct
                    (
                        &cost, agentState_.t_[i],
                        &agentState_.x_[i * Nxi],
                        &agentState_.u_[i * Nui],
                        &local_copies.x_[i * Nxj],
                        &local_copies.u_[i * Nuj]
                    );
                }
                else if (optimizationInfo_.COMMON_Solver_ == "Sensi")
                {
                    const auto& neighbors_agentState = neighbor->get_neighbors_agentState();
                    neighbor->get_couplingModel()->lfct
                    (
                        &cost, agentState_.t_[i],
                        &agentState_.x_[i * Nxi],
                        &agentState_.u_[i * Nui],
                        &neighbors_agentState.x_[i * Nxj],
                        &neighbors_agentState.u_[i * Nuj]
                    );

                }              
            }
        }

	    // multiply sum with timestep
	    cost *= agentState_.t_[1] - agentState_.t_[0];

	    // consider terminal cost
	    get_agentModel()->Vfct(&cost, agentState_.t_[Nhor - 1], &agentState_.x_[Nxi * (Nhor - 1)], &desired_agentState_.x_[0]);

		for (const auto& neighbor : get_sendingNeighbors())
		{
			const auto Nxj = neighbor->get_Nxj();
			const auto Nuj = neighbor->get_Nuj();
            if (optimizationInfo_.COMMON_Solver_ == "ADMM")
            {
                const auto& local_copies = neighbor->get_localCopies();
                neighbor->get_couplingModel()->Vfct
                (
                    &cost, agentState_.t_[Nhor - 1],
                    &agentState_.x_[(Nhor - 1) * Nxi],
                    &local_copies.x_[(Nhor - 1) * Nxj]
                );
            }
            else if (optimizationInfo_.COMMON_Solver_ == "Sensi")
            {
                const auto& neighbors_agentState = neighbor->get_neighbors_agentState();
                neighbor->get_couplingModel()->Vfct
                (
                    &cost, agentState_.t_[Nhor - 1],
                    &agentState_.x_[(Nhor - 1) * Nxi],
                    &neighbors_agentState.x_[(Nhor - 1) * Nxj]
                );
            }
		}
	    return cost;
    }

    void Agent::increase_all_delays(const AlgStep& step) const
    {
        for (auto& neighbor : sending_neighbors_)
            neighbor->increase_delays(step);
    }

    void Agent::initialize_allNeighborDelays() const
    {
        for (const auto& neighbor : neighbors_)
            neighbor->initialize_delays();
    }

    const int Agent::get_delay_sending_neighbors(const AlgStep& step) const 
    {
        int max_delay = 0;

        for (auto& neighbor : sending_neighbors_)
        {
            int previous_max_delay = max_delay;

            max_delay = std::max(max_delay, neighbor->get_delays(step));

            // ignore delay of finished neighbors
            if(neighbor->flag_StoppedAlg_)
                max_delay = previous_max_delay;
        }
        return max_delay;
    }

    const int Agent::get_delay_receiving_neighbors(const AlgStep& step) const
    {
        int max_delay = 0;

        for (auto& neighbor : receiving_neighbors_)
        {
            int previous_max_delay = max_delay;
            
            max_delay = std::max(max_delay, neighbor->get_delays(step));

            // ignore delay of finished neighbors
            if (neighbor->flag_StoppedAlg_)
                max_delay = previous_max_delay;
        }
        return max_delay;
    }

    const StepSelectorPtr& Agent::get_stepSelector() const
    {
        return step_selector_;
    }

    void Agent::reset_stopAdmmflag_of_neighbors()
    {
        for (auto& neighbor : neighbors_)
            neighbor->flag_StoppedAlg_ = false;
    }

}

