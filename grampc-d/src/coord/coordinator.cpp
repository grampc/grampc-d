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

#include "grampcd/coord/coordinator.hpp"

#include "grampcd/optim/optim_util.hpp"

#include "grampcd/util/logging.hpp"
#include "grampcd/util/data_conversion.hpp"

#include "grampcd/comm/communication_interface.hpp"

#include "grampcd/info/agent_info.hpp"
#include "grampcd/info/coupling_info.hpp"

namespace grampcd
{

    Coordinator::Coordinator(const CommunicationInterfacePtr& communication_interface, bool simulation, const LoggingPtr& log)
        : communication_interface_(communication_interface),
          simulation_(simulation),
        log_(log)
    {
    }

    const bool Coordinator::register_agent(const AgentInfo& agent_info)
    {
        // check if agent is already registered
        if(agents_.find(agent_info.id_) != agents_.end())
        {
            log_->print(DebugType::Error) << "Failed to register agent " << agent_info.id_ 
                << ", because it already exists." << std::endl;

		    return false;
        }

        // add agent to list
        agents_.insert(std::make_pair( agent_info.id_, std::make_shared<AgentInfo>(agent_info) ));

        // prepare entry in sending neighbors
	    sending_neighbors_.insert(std::make_pair(agent_info.id_, std::vector< std::shared_ptr< CouplingInfo> >()));

	    // prepare entry in receiving neighbors
	    receiving_neighbors_.insert(std::make_pair(agent_info.id_, std::vector< std::shared_ptr< CouplingInfo> >()));

        log_->print(DebugType::Message) << "[Coordinator::register_agent] "
            << "Agent " << agent_info.id_ << " registered. "
            << "Model name '" << agent_info.model_name_ << "'." << std::endl;

        return true;
    }

    const bool Coordinator::register_coupling(const CouplingInfo& coupling_info)
    {
        // check if agent is known
	    const bool agent_is_known = agents_.find(coupling_info.agent_id_) != agents_.end();

	    // check if neighbor is known
	    const bool neighbor_is_known = agents_.find(coupling_info.neighbor_id_) != agents_.end();
    
        if( !agent_is_known || !neighbor_is_known )
        {
            if (!agent_is_known && !neighbor_is_known)
                log_->print(DebugType::Warning) << "[Coordinator::register_coupling] Failed to register coupling as agents "
                << coupling_info.agent_id_ << " and " << coupling_info.neighbor_id_ << " are not known." << std::endl;
            else
            {
                log_->print(DebugType::Warning) << "[Coordinator::register_coupling] "
                    << "Failed to register coupling as agent ";
                if (!agent_is_known)
                    log_->print(DebugType::Warning) << coupling_info.agent_id_;
                else
                    log_->print(DebugType::Warning) << coupling_info.neighbor_id_;
                log_->print(DebugType::Warning) << " is not known." << std::endl;
            }
            return false;
        }

	    // check if coupling already exists 
        auto& sending_neighbors = sending_neighbors_.find(coupling_info.agent_id_)->second;
        auto& receiving_neighbors = receiving_neighbors_.find(coupling_info.neighbor_id_)->second;

        if(DataConversion::is_element_in_vector(sending_neighbors, coupling_info) ||
            DataConversion::is_element_in_vector(receiving_neighbors, coupling_info) )
        {
            log_->print(DebugType::Warning) << "[Coordinator::register_coupling] Failed to register coupling of agent " << coupling_info.agent_id_
                << " with neighbor " << coupling_info.neighbor_id_ << " as coupling is already registered." << std::endl;
            return false;
        }

        const CouplingInfoPtr info(std::make_shared<CouplingInfo>(coupling_info));

        // register coupling for sending neighbor
        sending_neighbors.push_back(info);

        // register coupling for receiving neighbor
        receiving_neighbors.push_back(info);

        log_->print(DebugType::Message) << "[Coordinator::register_coupling] "
            << "Coupling between agent " << coupling_info.agent_id_ << " (receiving) and agent " << coupling_info.neighbor_id_ << " (sending) registered. "
            << "Model name '" << coupling_info.model_name_ << "'." << std::endl;

        return true;
    }

    const bool Coordinator::deregister_agent(const AgentInfo &agent_info)
    {
        // check if agent is known
        if (agents_.find(agent_info.id_) == agents_.end())
        {
            log_->print(DebugType::Warning) << "[Coordinator::deregister_agent] Failed to deregister agent "
                << agent_info.id_ << " as agent is not known." << std::endl;

            return false;
        }

        const auto& sending_neighbors = sending_neighbors_.find(agent_info.id_)->second;
        const auto& receiving_neighbors = receiving_neighbors_.find(agent_info.id_)->second;

        // deregister couplings of sending neighbors
        while( sending_neighbors.size() > 0 )
            deregister_coupling(sending_neighbors[0]);

        // deregister couplings of receiving neigbors
        while( receiving_neighbors.size() > 0 )
            deregister_coupling(receiving_neighbors[0]);

        agents_.erase(agents_.find(agent_info.id_));

        return true;
    }

    const bool Coordinator::deregister_coupling(CouplingInfoPtr coupling_info)
    {
	    auto& sending_neighbors = sending_neighbors_.find(coupling_info->agent_id_)->second;
        auto& receiving_neighbors = receiving_neighbors_.find(coupling_info->neighbor_id_)->second;

        // check if coupling is registered in sending neighbors
        if (!DataConversion::is_element_in_vector(sending_neighbors, *coupling_info))
        {
            log_->print(DebugType::Error) << "[Coordinator::deregister_coupling] "
                << "Failed to deregister coupling, as agent " << coupling_info->agent_id_
                << " has no sending neighbor " << coupling_info->neighbor_id_ << "." << std::endl;
		    return false;
	    }

	    // check if coupling is registered in receiving neighbors
	    if (!DataConversion::is_element_in_vector(receiving_neighbors, *coupling_info))
	    {
            log_->print(DebugType::Error) << "[Coordinator::deregister_coupling] "
                << "Failed to deregister coupling, as agent " << coupling_info->neighbor_id_
                << " has no receiving neighbor " << coupling_info->agent_id_ << "." << std::endl;
		    return false;
	    }

        // inform both agents that coupling is deregistered
        communication_interface_->fromCommunication_deregistered_coupling(*coupling_info);

        // and delete coupling
        DataConversion::erase_element_from_vector(sending_neighbors, *coupling_info);
        DataConversion::erase_element_from_vector(receiving_neighbors, *coupling_info);

        log_->print(DebugType::Message) << "[Coordinator::deregister_coupling] "
            << "Coupling between agent " << coupling_info->agent_id_ << " (receiving) and agent " << coupling_info->neighbor_id_
            << " (sending) deregistered." << std::endl;

        return true;
    }

    void Coordinator::print_network() const
    {
        for( const auto& [id, info] : agents_ )
        {
            log_->print(DebugType::Message) << "Agent " << id << ": sending neighbors [";
            for (const auto& coupling : sending_neighbors_.find(id)->second)
                log_->print(DebugType::Message) << coupling->neighbor_id_ << ", ";

            log_->print(DebugType::Message) << "] receiving neighbors [";
            for (const auto& coupling : receiving_neighbors_.find(id)->second)
                log_->print(DebugType::Message) << coupling->agent_id_ << ", ";

            log_->print(DebugType::Message) << "]" << std::endl;
        }
    }

    const unsigned int Coordinator::get_numberOfAgents() const
    {
        return static_cast<unsigned int>(agents_.size());
    }

    /*************************************************************************
     coordination of alternating direction method of multipliers
     *************************************************************************/

    void Coordinator::initialize_ADMM(const OptimizationInfo& oi)
    {
        optimizationInfo_ = oi;

        // initialize local solvers
        communication_interface_->configure_optimization(oi);

        // send initial agent states
        communication_interface_->trigger_step(ADMMStep::SEND_AGENT_STATE);

        // send initial coupling states
        communication_interface_->trigger_step(ADMMStep::SEND_COUPLING_STATE);

        // send initial multiplier states
        communication_interface_->trigger_step(ADMMStep::SEND_MULTIPLIER_STATE);
    }

    void Coordinator::solve_ADMM(int outer_iterations, int inner_iterations)
    {
        communication_interface_->trigger_step( ADMMStep::INITIALIZE );

        for(int i = 0; i < outer_iterations; ++i)
        {
            for(int j = 0; j < inner_iterations; ++j)
            {
                // solve local minimization problem for agent states
                communication_interface_->trigger_step(ADMMStep::UPDATE_AGENT_STATE);

                // send updated agent states to receiving neighbors
                communication_interface_->trigger_step(ADMMStep::SEND_AGENT_STATE);

                // solve local minimization problem for coupling states
                communication_interface_->trigger_step(ADMMStep::UPDATE_COUPLING_STATE);

                // send updated coupling states to sending neighbors
                communication_interface_->trigger_step(ADMMStep::SEND_COUPLING_STATE);
            }

            // solve local maximization problem for multiplier states
            communication_interface_->trigger_step(ADMMStep::UPDATE_MULTIPLIER_STATE);

            // send updated multiplier states to receiving neighbors
            communication_interface_->trigger_step(ADMMStep::SEND_MULTIPLIER_STATE);

            if(optimizationInfo_.ADMM_DebugCost_)
                communication_interface_->trigger_step(ADMMStep::PRINT);

            // evaluate convergence
            ADMM_converged_ = true;
            communication_interface_->trigger_step(ADMMStep::SEND_CONVERGENCE_FLAG);
            if(ADMM_converged_)
                break;
        }
    }

    void Coordinator::fromCommunication_received_convergenceFlag(bool converged, int from)
    {
        ADMM_converged_ = ADMM_converged_ && converged;
    }

    void Coordinator::trigger_simulation(const std::string& Integrator, typeRNum dt) const
    {
        communication_interface_->trigger_simulation(Integrator, dt);
    }

    const std::map<unsigned int, AgentInfoPtr >& Coordinator::get_agentInfos()const
    {
        return agents_;
    }

    const std::map< unsigned int, std::vector< CouplingInfoPtr > >& Coordinator::get_sendingNeighbors() const
    {
        return sending_neighbors_;
    }

    const std::map< unsigned int, std::vector< CouplingInfoPtr > >& Coordinator::get_receivingNeighbors() const
    {
        return receiving_neighbors_;
    }

    const OptimizationInfo& Coordinator::get_optimizationInfo() const
    {
        return optimizationInfo_;
    }

}
