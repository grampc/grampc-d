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

#include "grampcd/comm/communication_interface_central.hpp"

#include "grampcd/model/agent_model.hpp"

#include "grampcd/agent/agent.hpp"

#include "grampcd/coord/coordinator.hpp"

#include "grampcd/optim/solution.hpp"

#include "grampcd/simulator/simulator.hpp"

#include "grampcd/util/logging.hpp"

#include <algorithm>

namespace grampcd
{

    CommunicationInterfaceCentral::CommunicationInterfaceCentral(const LoggingPtr& log, const int number_of_threads) 
		:
		work_(ioService_),
        log_(log)
	{
		// start threads
        if (number_of_threads > 1)
        {
            for (int i = 0; i < number_of_threads; ++i)
                threads_.push_back(std::thread([this]() {ioService_.run(); }));
        }

        log_->print(DebugType::Message) << "[CentralizedCommunicationInterface::CommunicationInterfaceCentral] "
            << "Communication interface is running." << std::endl;
    }

void CommunicationInterfaceCentral::set_coordinator(const CoordinatorPtr& coordinator)
{
    coordinator_ = coordinator;
}

void CommunicationInterfaceCentral::set_simulator(const SimulatorPtr& simulator)
{
    simulator_ = simulator;
}

const bool CommunicationInterfaceCentral::register_agent(const AgentPtr& agent)
{
    if(!coordinator_)
    {
        log_->print(DebugType::Error) << "[CentralizedCommunicationInterface::register_agent] "
            << "Missing coordinator, call setCoordinator(...) first." << std::endl;
        return false;
    }
    
    if(!coordinator_->register_agent(agent->get_agentInfo()))
    {
        log_->print(DebugType::Error) << "[CentralizedCommunicationInterface::registerState] "
            << "Could not register agent with id " << agent->get_id() << "." << std::endl;
        return false; 
    }

    // add agent to list
    if (agent->get_id() >= agents_.size())
        agents_.resize(agent->get_id() + 1);
    agents_[agent->get_id()] = agent;

    return true;
}

const bool CommunicationInterfaceCentral::deregister_agent(const AgentInfo& agent)
{
	if (!coordinator_)
	{
        log_->print(DebugType::Error) << "[CentralizedCommunicationInterface::deregister_agent] "
            << "Missing coordinator, call setCoordinator(...) first." << std::endl;
		return false;
	}

    if (agent.id_ > agents_.size())
    {
        log_->print(DebugType::Error) << "[CommunicationInterfaceCentral::deregister_agent] "
            << "Failed to deregister agent as agent is not known." << std::endl;
        return false;
    }

	if (!coordinator_->deregister_agent(agent))
	{
        log_->print(DebugType::Error) << "[CentralizedCommunicationInterface::registerState] "
            << "Could not deregister agent with id '" << agent.id_ << "'" << std::endl;
		return false;
	}

    // delete agent from list
    agents_[agent.id_] = nullptr;

    return true;
}

const bool CommunicationInterfaceCentral::register_coupling(const CouplingInfo& coupling)
{
    if(!coordinator_)
    {
        log_->print(DebugType::Error) << "[CentralizedCommunicationInterface::register_coupling] "
            << "Failed to register coupling as coordinator is missing." << std::endl;
        return false;
    }

    if(!coordinator_->register_coupling(coupling))
    {
        log_->print(DebugType::Error) << "[CentralizedCommunicationInterface::register_coupling] "
            << "Could not register coupling between " << coupling.agent_id_ << " and "
            << coupling.neighbor_id_ << "." << std::endl;
        return false; 
    }

    if (std::max(coupling.agent_id_, coupling.neighbor_id_) >= agents_.size())
    {
        log_->print(DebugType::Error) << "[CentralizedCommunicationInterface::register_coupling] "
            << "Could not register coupling between " << coupling.agent_id_ << " and "
            << coupling.neighbor_id_ << "." << std::endl;
        return false;
    }
    
    const AgentPtr& agent = agents_[coupling.agent_id_];
    const AgentPtr& neighbor = agents_[coupling.neighbor_id_];

    agent->fromCommunication_registered_coupling(coupling, neighbor->get_agentInfo());
    neighbor->fromCommunication_registered_coupling(coupling, agent->get_agentInfo());

    return true;
}

const bool CommunicationInterfaceCentral::fromCommunication_deregistered_coupling(const CouplingInfo &coupling)
{
    if (std::max(coupling.agent_id_, coupling.neighbor_id_) >= agents_.size())
    {
        log_->print(DebugType::Error) << "[CentralizedCommunicationInterface::register_coupling] "
            << "Could not deregister coupling between " << coupling.agent_id_ << " and "
            << coupling.neighbor_id_ << "." << std::endl;
        return false;
    }

    agents_[coupling.agent_id_]->fromCommunication_deregistered_coupling(coupling);
    agents_[coupling.neighbor_id_]->fromCommunication_deregistered_coupling(coupling);
    return true;
}

const bool CommunicationInterfaceCentral::deregister_coupling(const CouplingInfo& coupling)
{
    if(!coordinator_)
    {
        log_->print(DebugType::Error) << "[CentralizedCommunicationInterface::deregister_coupling] "
            << "Missing coordinator, call setCoordinator(...) first." << std::endl;
        return false;
    }

    return coordinator_->deregister_coupling(std::make_shared<CouplingInfo>(coupling));
}

const bool CommunicationInterfaceCentral::send_numberOfNeighbors(const int number, const int from, const int to)
{
    if (to >= agents_.size())
    {
        log_->print(DebugType::Warning) << "[CommunicationInterfaceCentral::send_numberOfNeighbors] "
            << "Could not find agent with id " << to << "." << std::endl;
        return false;
    }

    agents_[to]->fromCommunication_received_numberOfNeighbors(number, from);
    return true;    
}

const bool CommunicationInterfaceCentral::send_localCopies(const AgentState& state, const int from, const int to)
{
    if( to >= agents_.size() )
    {
        log_->print(DebugType::Warning) << "[CentralizedCommunicationInterface::send_localCopies] "
            << "Could not find agent with id " << to << "." << std::endl;
        return false;
    }

    agents_[to]->fromCommunication_received_localCopies(state, from);
    return true;
}

const bool CommunicationInterfaceCentral::send_agentState(const AgentState& state, const ConstraintState& constr_state, const int from, const int to)
{
    if (to >= agents_.size())
    {
        log_->print(DebugType::Warning) << "[CentralizedCommunicationInterface::send_agentState] "
            << "Could not find agent with id " << to << "." << std::endl;
        return false;
    }

    agents_[to]->fromCommunication_received_agentState(state, constr_state, from);
    return true;
}

const bool CommunicationInterfaceCentral::send_desiredAgentState(const AgentState &desired_state, const int from, const int to)
{
    if( to >= agents_.size() )
    {
        log_->print(DebugType::Warning) << "[CentralizedCommunicationInterface::send_desiredAgentState] "
            << "Could not find agent with id " << to << "." << std::endl;
        return false;
    }

    agents_[to]->fromCommunication_received_desiredAgentState(desired_state, from);
    return true;
}

const std::shared_ptr< AgentState > CommunicationInterfaceCentral::get_agentState_from_agent(const int agentId) const
{
    return std::make_shared<AgentState>(agents_[agentId]->get_agentState());
}

const std::shared_ptr< AgentState > CommunicationInterfaceCentral::get_desiredAgentState_from_agent(const int agentId) const
{
    return std::make_shared<AgentState>(agents_[agentId]->get_desiredAgentState());
}

const std::shared_ptr< std::map<unsigned int, AgentInfoPtr > > CommunicationInterfaceCentral::get_agentInfo_from_coordinator() const
{
    return std::make_shared< std::map<unsigned int, AgentInfoPtr > >(coordinator_->get_agentInfos());
}

const std::shared_ptr< std::map< unsigned int, std::vector< CouplingInfoPtr > > > CommunicationInterfaceCentral::get_sendingNeighbors_from_coordinator() const
{
    return std::make_shared< std::map< unsigned int, std::vector< CouplingInfoPtr > > >(coordinator_->get_sendingNeighbors());
}

const std::shared_ptr< std::map< unsigned int, std::vector< CouplingInfoPtr > > > CommunicationInterfaceCentral::get_receivingNeighbors_from_coordinator() const
{
    return std::make_shared< std::map< unsigned int, std::vector< CouplingInfoPtr > > >(coordinator_->get_receivingNeighbors());
}

const AgentModelPtr CommunicationInterfaceCentral::get_agentModel(const int agentId) const
{
    return agents_[ agentId ]->get_agentModel();
}

const std::shared_ptr< std::map<int, CouplingModelPtr> > CommunicationInterfaceCentral::get_couplingModels_from_agent(const int agentId ) const
{
    return agents_[agentId]->get_couplingModels();
}

const bool CommunicationInterfaceCentral::send_couplingState(const CouplingState& state, const int from, const int to)
{
    if( to >= agents_.size() )
    {
        log_->print(DebugType::Warning) << "[CentralizedCommunicationInterface::sendCouplingState] "
            << "Could not find agent with id " << to << "." << std::endl;
        return false;
    }

    agents_[to]->fromCommunication_received_couplingState(state, from);
    return true;
}

const bool CommunicationInterfaceCentral::send_couplingState(const CouplingState& state, const CouplingState& state2, const int from, const int to)
{
    if( to >= agents_.size() )
    {
        log_->print(DebugType::Warning) << "[CentralizedCommunicationInterface::sendCouplingState] "
            << "Could not find agent with id " << to << "." << std::endl;
        return false;
    }

    agents_[to]->fromCommunication_received_couplingState(state, state2, from);
    return true;
}

const bool CommunicationInterfaceCentral::send_multiplierState(const MultiplierState& state, PenaltyState penalty, const int from, const int to)
{
    if( to >= agents_.size() )
    {
        log_->print(DebugType::Warning) << "[CentralizedCommunicationInterface::sendMultiplierState] "
            << "Could not find agent with id " << to << "." << std::endl;
        return false;
    }

    agents_[to]->fromCommunication_received_multiplierState(state, penalty, from);
    return true;
}

const bool CommunicationInterfaceCentral::send_convergenceFlag(const bool converged, const int from)
{
    coordinator_->fromCommunication_received_convergenceFlag(converged, from);
    return true;
}

const bool  CommunicationInterfaceCentral::send_stoppedAlgFlag(const bool flag, const int from, const int to)
{
    if (to >= agents_.size())
    {
        log_->print(DebugType::Warning) << "[CentralizedCommunicationInterface::send_admmIter] "
            << "Could not find agent with id " << to << "." << std::endl;
        return false;
    }
  
    agents_[to]->fromCommunication_received_flagStoppedAdmm(flag, from);
    return true;
}

const bool  CommunicationInterfaceCentral::send_stoppedAlgFlag(const bool flag, const int from)
{
    coordinator_->fromCommunication_received_flagStoppedAlg(flag, from);
    return true;
}

const bool CommunicationInterfaceCentral::configure_optimization(const OptimizationInfo& info)
{
    // configure optimization for each agent
    for( const AgentPtr& agent : agents_ )
    {
        if( agent != nullptr )
            agent->fromCommunication_configured_optimization(info);
    }
    return true;
}

const bool CommunicationInterfaceCentral::trigger_step(const AlgStep& step)
{
    // if number of threads is chosen to be 1, simply trigger steps
    if (threads_.size() == 0)
    {
		for (const auto& agent : agents_)
		{
			if (agent != nullptr)
				agent->fromCommunication_trigger_step(step);
		}
		return true;
	}

    // if multiple threads are available, start trigger step
    // as asynchronous task
    std::unique_lock<std::mutex> guard(mutex_triggerStep_);
    counter_triggerStep_ = 0;
	for (const auto& agent : agents_)
	{
		if (agent != nullptr)
		{
			++counter_triggerStep_;
			asio::post([agent, step, this]()
			{
				agent->fromCommunication_trigger_step(step);
				std::unique_lock<std::mutex> guard(mutex_triggerStep_);
				--counter_triggerStep_;
				condVar_triggerStep_.notify_one();
			});
		}
	}
	condVar_triggerStep_.wait(guard, [this]() { return counter_triggerStep_ == 0; });

    return true;
}

void CommunicationInterfaceCentral::trigger_simulation(const std::string& Integrator, const typeRNum dt)
{
    simulator_->distributed_simulation(Integrator, dt);
}

void CommunicationInterfaceCentral::set_simulatedState_for_agent
(
    const int agentId, 
    const std::vector<typeRNum>& new_state, 
    const typeRNum dt, 
	const typeRNum t0,
	const typeRNum cost
)
{
    agents_[ agentId ]->set_updatedState( new_state, dt, t0, cost);
}

const SolutionPtr CommunicationInterfaceCentral::get_solution(const unsigned int agent_id) const
{
	for (auto agent : agents_)
	{
        if (agent->get_id() == agent_id)
            return agent->get_solution();
	}

	return nullptr;
}

const std::vector<SolutionPtr> CommunicationInterfaceCentral::get_solution(const std::string& agents) const
{
    if (agents != "all")
    {
        log_->print(DebugType::Error) << "[CommunicationInterfaceCentral::get_solution] "
            << "Unknown set of agents." << std::endl;

        return std::vector<SolutionPtr>();
    }

    if (!coordinator_)
    {
        log_->print(DebugType::Error) << "[CommunicationInterfaceCentral::get_solution] "
            << "Coordinator is required to sample all solutions." << std::endl;

		return std::vector<SolutionPtr>();
    }

    std::vector<SolutionPtr> solutions;
    for (auto agent : agents_)
    {
        if(agent != nullptr)
            solutions.push_back(agent->get_solution());
    }

    return solutions;
}

void CommunicationInterfaceCentral::reset_solution(const unsigned int agent_id)
{
    for (const auto& agent : agents_)
    {
        if (agent->get_id() == agent_id)
        {
            agent->reset_solution();
            return;
        }
    }
}

void CommunicationInterfaceCentral::reset_solution(const std::string& agents)
{
    if (agents == "all")
        for (const auto& agent : agents_)
            agent->reset_solution();
    else
        log_->print(DebugType::Error) << "[CommunicationInterfaceCentral::reset_solution] "
        << "Set of agents is not known." << std::endl;
}

void CommunicationInterfaceCentral::waitFor_connection(const int agents, const int couplings)
{
    log_->print(DebugType::Error) << "[CommunicationInterfaceCentral::waitForConnection] "
        << "This function is not implemented for the centralized communication interface." << std::endl;
}

void CommunicationInterfaceCentral::set_passive()
{
    log_->print(DebugType::Error) << "[CommunicationInterfaceCentral::set_passive] "
        << "This function is not implemented for the centralized communication interface." << std::endl;
}

void CommunicationInterfaceCentral::send_flag_to_agents(const int agent_id) const
{
    log_->print(DebugType::Error) << "[CommunicationInterfaceCentral::send_flag_to_agents] "
        << "This function is not implemented for the centralized communication interface." << std::endl;
}

void CommunicationInterfaceCentral::send_flag_to_agents(const std::vector<int>& agent_ids) const
{
    log_->print(DebugType::Error) << "[CommunicationInterfaceCentral::send_flag_to_agents] "
        << "This function is not implemented for the centralized communication interface." << std::endl;
}

void CommunicationInterfaceCentral::send_flag_to_agents(const std::string& agents) const
{
    log_->print(DebugType::Error) << "[CommunicationInterfaceCentral::send_flag_to_agents] "
        << "This function is not implemented for the centralized communication interface." << std::endl;
}

void CommunicationInterfaceCentral::waitFor_flag_from_coordinator()
{
    log_->print(DebugType::Error) << "[CommunicationInterfaceCentral::waitFor_flag_from_coordinator] "
        << "This function is not implemented for the centralized communication interface." << std::endl;
}

void CommunicationInterfaceCentral::cap_stored_data(const unsigned int data_points)
{
    for (auto agent : agents_)
    {
        if(agent != nullptr)
            agent->get_solution()->maximum_number_of_data_points_ = data_points;
    }
}

const unsigned int CommunicationInterfaceCentral::get_numberOfAgents() const
{
    if (!coordinator_)
    {
        log_->print(DebugType::Error) << "CommunicationInterfaceCentral::get_numberOfAgents"
            << "This method requires the coordinator." << std::endl;
        return 1;
    }

    return coordinator_->get_numberOfAgents();
}

}
