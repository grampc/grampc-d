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
#include "grampcd/optim/solution.hpp"

namespace grampcd
{

    Agent::Agent(const CommunicationInterfacePtr& communication_interface,
                 const ModelFactoryPtr& model_factory,
                 const AgentInfo& agent_info,
        const LoggingPtr& log)
        :   
        communication_interface_(communication_interface),
        model_factory_(model_factory),
        model_(nullptr),
        info_(agent_info),
        local_solver_(nullptr),
        solution_(new Solution()),
        log_(log)
    {
        model_ = model_factory->create_agentModel(agent_info);
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

        resetState(penaltyState_, id, t);

        // resize states
        agentState_.x_.resize( Nhor * Nxi, 0.0 );
        agentState_.u_.resize( Nhor * Nui, 0.0 );
        couplingState_.z_u_.resize( Nhor * Nui, 0.0 );
        multiplierState_.mu_u_.resize( Nhor * Nui, 0.0 );
	    penaltyState_.rho_u_.resize(Nhor * Nui, initial_penalty_);
        desired_agentState_.x_.resize(Nhor * Nxi, 0.0);
        desired_agentState_.u_.resize(Nhor * Nui, 0.0);

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

        // initialize coupling states
        couplingState_.z_u_ = agentState_.u_;
        if( !is_approximating_ )
            couplingState_.z_x_ = agentState_.x_;

        // initialize previous states
        previous_couplingState_ = couplingState_;
        previous_multiplierState_ = multiplierState_;

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

    void Agent::set_agentState(const AgentState& state)
    {
        if(!compare_stateDimensions( agentState_, state ))
        {
            log_->print(DebugType::Error) << "[Agent::set_agentState] Agent " << get_id() << ": "
                << "Failed to set agent state, as dimensions don't fit." << std::endl;

            return;
	    }

        // set agent state
	    agentState_ = state;

        // evaluate predicted cost
	    std::vector< typeRNum > predicted_cost(optimizationInfo_.COMMON_Nhor_, 0.0);
	    for (unsigned int i = 0; i < predicted_cost.size(); ++i)
            model_->lfct(&predicted_cost[i], agentState_.t_[i], &agentState_.x_[i], &agentState_.u_[i], &desired_agentState_.x_[i]);

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

    void Agent::fromCommunication_received_agentState(const AgentState& state, int from)
    {
        if (state.i_ != get_id())
        {
			log_->print(DebugType::Warning) << "[Agent::fromCommunication_received_agentState] Agent " << get_id() << ": "
				<< "Agent " << get_id() << " received state of agent " << state.i_ << " from agent "
				<< from << "." << std::endl;

            return;
        }
        
		const auto neighbor = DataConversion::get_element_from_vector(neighbors_, from);

        if (neighbor == nullptr)
        {
			log_->print(DebugType::Warning) << "[Agent::fromCommunication_received_agentState] Agent " << get_id() << ": "
				<< "Could not find neighbor with id " << from << "." << std::endl;

            return;
        }
        
        neighbor->set_neighbors_localCopies(state);
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
    }

    void Agent::fromCommunication_configured_optimization(const OptimizationInfo& info)
    {
        initialize(info);

        // create local solver
        local_solver_.reset(new SolverLocal(this, info, log_));
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

    void Agent::fromCommunication_trigger_step(const ADMMStep& step)
    {
        switch(step)
        {

        case ADMMStep::UPDATE_AGENT_STATE:
            local_solver_->update_agentStates();
            break;

        case ADMMStep::SEND_AGENT_STATE:
            // send local copies to neighbors
            for( const auto& neighbor : neighbors_ )
            {
                if( neighbor->is_sendingNeighbor() || neighbor->is_approximating() )
                {
                    if (neighbor->is_approximatingDynamics())
                    {
                        AgentState local_copies = neighbor->get_localCopies();
                        local_copies.x_.clear();

                        communication_interface_->send_agentState(local_copies, get_id(), neighbor->get_id());
                    }
                    else
                        communication_interface_->send_agentState(neighbor->get_localCopies(), get_id(), neighbor->get_id());
                }
            }
            break;

        case ADMMStep::UPDATE_COUPLING_STATE:
            local_solver_->update_couplingStates();
            break;

        case ADMMStep::SEND_COUPLING_STATE:
            // send coupling state to neighbors
            {
                for( const auto& neighbor : neighbors_ )
                {
                    // send coupling state
                    if( ( neighbor->is_receivingNeighbor() || neighbor->is_approximating() ) && ! neighbor->is_approximatingDynamics() )
                        communication_interface_->send_couplingState(get_couplingState(), get_id(), neighbor->get_id());

                    // send coupling state and ext_influence_coupling_state
                    if( neighbor->is_approximatingDynamics() )
                        communication_interface_->send_couplingState(get_couplingState(), neighbor->get_externalInfluence_couplingState(), get_id(), neighbor->get_id());
                }
            }
            break;

        case ADMMStep::UPDATE_MULTIPLIER_STATE:
            local_solver_->update_multiplierStates();
            break;

        case ADMMStep::SEND_MULTIPLIER_STATE:
            // send multiplier state to all sending neighbors
            for( const auto& neighbor : neighbors_ )
            {
                if( neighbor->is_sendingNeighbor() || neighbor->is_approximating() )
                    communication_interface_->send_multiplierState(neighbor->get_coupled_multiplierState(),
                        neighbor->get_coupled_penaltyState(), get_id(), neighbor->get_id());
            }
            break;

        case ADMMStep::SEND_CONVERGENCE_FLAG:
            communication_interface_->send_convergenceFlag(local_solver_->is_converged(), get_id());
            break;

        case ADMMStep::INITIALIZE:
            // send number of neighbors to each neighbor
            for( const auto& neighbor : get_neighbors() )
                communication_interface_->send_numberOfNeighbors(static_cast<int>(neighbors_.size()), get_id(), neighbor->get_id());

            local_solver_->initialize_ADMM();
            set_neighbors_initial_states();
            break;

        case ADMMStep::PRINT:
            solution_->update_debug_cost(get_predicted_cost());
            break;
        }
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
        }

	    // multiply sum with timestep
	    cost *= agentState_.t_[1] - agentState_.t_[0];

	    // consider terminal cost
	    get_agentModel()->Vfct(&cost, agentState_.t_[Nhor - 1], &agentState_.x_[Nxi * (Nhor - 1)], &desired_agentState_.x_[0]);

		for (const auto& neighbor : get_sendingNeighbors())
		{
			const auto Nxj = neighbor->get_Nxj();
			const auto Nuj = neighbor->get_Nuj();
			const auto& local_copies = neighbor->get_localCopies();

			neighbor->get_couplingModel()->Vfct
			(
				&cost, agentState_.t_[Nhor - 1],
				&agentState_.x_[(Nhor - 1) * Nxi],
				&local_copies.x_[(Nhor - 1) * Nxj]
			);
		}

	    return cost;
    }

}
