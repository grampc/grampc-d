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

#include "grampcd/comm/communication_interface_local.hpp"
#include "grampcd/comm/message.hpp"
#include "grampcd/comm/message_handler.hpp"

#include "grampcd/info/communication_data.hpp"

#include "general_model_factory.hpp"

#include "grampcd/simulator/simulator.hpp"

#include "grampcd/agent/agent.hpp"

#include "grampcd/coord/coordinator.hpp"

#include "grampcd/util/logging.hpp"
#include "grampcd/util/data_conversion.hpp"
#include "grampcd/util/constants.hpp"

#include "grampcd/optim/optim_util.hpp"

#include "asio.hpp"

#include <cereal/archives/binary.hpp>
#include <iostream>
#include <sstream>
#include <chrono>
#include <charconv>
#include <thread>

namespace grampcd
{
    struct EncapsulatedAsio
    {
        EncapsulatedAsio(const unsigned short port):
			acceptor_(ioService_, asio::ip::tcp::endpoint(asio::ip::tcp::v4(), port)),
			timer_waitForAck_(ioService_),
			timer_waitTrue_(ioService_)
        {}

		asio::io_service ioService_;
		asio::ip::tcp::acceptor acceptor_;
		asio::basic_waitable_timer<std::chrono::system_clock>  timer_waitForAck_;
		asio::basic_waitable_timer<std::chrono::system_clock>  timer_waitTrue_;
		std::vector<std::thread> threads_for_communication_;
    };

    // constructor for agents
    CommunicationInterfaceLocal::CommunicationInterfaceLocal(const LoggingPtr& log, const CommunicationInfo& comm_info_coordinator)
        :
		asio_(std::make_shared<EncapsulatedAsio>(0)),
		log_(log),
		message_handler_(new MessageHandler(this, log_))
    {
        comm_info_local_.agent_type_ = "agent";
        comm_info_local_.ip_ = asio_->acceptor_.local_endpoint().address().to_string();
        comm_info_local_.port_ = std::to_string(asio_->acceptor_.local_endpoint().port());

        // set coordinator
        std::unique_lock<std::shared_mutex> guard(mutex_coordinator_);
        comm_info_coordinator_ = comm_info_coordinator;

        // start communication as client
        start_client();

        log_->print(DebugType::Message) << "[CommunicationInterfaceLocal::CommunicationInterfaceLocal] "
            << "Communication interface is set up for an agent." << std::endl;
    }

    // constructor for coordinator
    CommunicationInterfaceLocal::CommunicationInterfaceLocal(const LoggingPtr& log, const unsigned short port)
		:
		asio_(std::make_shared<EncapsulatedAsio>(port)),
		log_(log),
		message_handler_(new MessageHandler(this, log_))
    {
        comm_info_local_.agent_type_ = "coordinator";
        comm_info_local_.ip_ = asio_->acceptor_.local_endpoint().address().to_string();
        comm_info_local_.port_ = std::to_string(asio_->acceptor_.local_endpoint().port());

        // start communication as server
        start_server();

        log_->print(DebugType::Message) << "[CommunicationInterfaceLocal::CommunicationInterfaceLocal] "
            << "Communication interface is set up for a coordinator." << std::endl;
    }

    CommunicationInterfaceLocal::~CommunicationInterfaceLocal()
    {
        // stop the ioService
        asio_->ioService_.stop();

        // join all threads
        for (auto& thread : asio_->threads_for_communication_)
            thread.join();
    }

    void CommunicationInterfaceLocal::handle_disconnect(const CommunicationDataPtr& comm_data)
    {
        // reset flag
        comm_data->is_connected_ = false;

        // reset timer
        comm_data->timer_ping_.cancel();
        comm_data->timer_polling_.cancel();
        comm_data->timer_check_connection_.cancel();

        // clear data array
        comm_data->data_.clear();

        // close socket
        close_shutdown_socket(comm_data);

        if (comm_info_local_.agent_type_ == "agent")
            handle_disconnect_as_agent(comm_data);
        else if (comm_info_local_.agent_type_ == "coordinator")
            handle_disconnect_as_coordinator(comm_data);
    }

    void CommunicationInterfaceLocal::handle_disconnect_as_agent(const CommunicationDataPtr& comm_data)
    {
        // An agent disconnected
        if(comm_data->communication_info_->agent_type_ == "agent")
        {
            log_->print(DebugType::Message) << "[CommunicationInterfaceLocal::handleDisconnectAsAgent]"
                << " Agent with id " << comm_data->communication_info_->id_
                << " disconnected." << std::endl;

            // remove corresponding communication data from the list        
            std::unique_lock<std::shared_mutex> guard(mutex_comm_data_vec_);
            DataConversion::erase_element_from_vector(comm_data_vec_, comm_data);
        }
        // The coordinator disconnected, so fully restart
        else if(comm_data->communication_info_->agent_type_ == "coordinator")
        {
            log_->print(DebugType::Message) << "[CommunicationInterfaceLocal::handleDisconnectAsAgent] "
                << "Coordinator disconnected." << std::endl;

            std::lock_guard<std::shared_mutex> guard(mutex_infos_);

            // move all active registration infos to pending
            while( agent_infos_active_.size() > 0 )
            {
                agent_infos_registering_.push_back(agent_infos_active_[0]);
                agent_infos_active_.erase(agent_infos_active_.begin());
            }

            // move all active coupling infos to pending and delete all corresponding neighbors
            std::unique_lock<std::shared_mutex> guard_mutex_basics(mutex_basics_);
            while(coupling_infos_active_.size() > 0)
            {
                const auto& info = coupling_infos_active_[0];

                agent_->fromCommunication_deregistered_coupling(info);

                // check if coupling should be re-registered
                const bool reregister = (!DataConversion::is_element_in_vector(coupling_infos_blocked_, info))
                    && (info.agent_id_ == agent_->get_id());

                // if it should not be re-registered, delete it form the list pending couplings
                if (!reregister)
                    DataConversion::erase_element_from_vector(coupling_infos_pending_, info);
                else if(reregister && !DataConversion::is_element_in_vector(coupling_infos_pending_, info))
                    coupling_infos_pending_.push_back(info);

                // delete coupling from the list deregistering couplings
                DataConversion::erase_element_from_vector(coupling_infos_deregistering_, info);

                // delete coupling from the list active couplings
                DataConversion::erase_element_from_vector(coupling_infos_active_, info);
            }

            // close all connections to agents
            std::unique_lock<std::shared_mutex> guard_comm_data_vec(mutex_comm_data_vec_);
            while(comm_data_vec_.size() > 1)
            {
                close_shutdown_socket(comm_data_vec_[1]);
                DataConversion::erase_element_from_vector(comm_data_vec_,comm_data_vec_[1]);
            }

            // try to reconnect
            log_->print(DebugType::Message) << "[CommunicationInterfaceLocal::handleDisconnectAsAgent] "
                << "Try to reconnect to coordinator." << std::endl;

            async_connect(comm_data, comm_data->communication_info_->ip_, comm_data->communication_info_->port_);
        }
    }

    void CommunicationInterfaceLocal::close_socket(const CommunicationDataPtr& comm_data) const
    {
        std::unique_lock<std::shared_mutex> guard_socket(comm_data->mutex_socket_);

        if (!comm_data->socket_.is_open())
            return;

        std::error_code ec;
        comm_data->socket_.close(ec);
        if (ec)
        {
		    log_->print(DebugType::Error) << "[CommunicationInterfaceLocal::close_socket] "
			    << "Failed to close socket with" << std::endl
			    << "Error code: " << ec.value() << std::endl
			    << "Message: " << ec.message() << std::endl;
        }
    }

    void CommunicationInterfaceLocal::close_shutdown_socket(const CommunicationDataPtr& comm_data) const
    {
        std::unique_lock<std::shared_mutex> guard_socket(comm_data->mutex_socket_);

        if (!comm_data->socket_.is_open())
            return;

        std::error_code ec;
        comm_data->socket_.shutdown(asio::ip::tcp::socket::shutdown_both, ec);
        if (ec)
        {
		    log_->print(DebugType::Error) << "[CommunicationInterfaceLocal::close_shutdown_socket] "
			    << "Failed to shutdown socket with" << std::endl
			    << "Error code: " << ec.value() << std::endl
			    << "Message: " << ec.message() << std::endl;
        }

        comm_data->socket_.close(ec);
        if (ec)
        {
		    log_->print(DebugType::Error) << "[CommunicationInterfaceLocal::close_shutdown_socket] "
			    << "Failed to close socket with" << std::endl
			    << "Error code: " << ec.value() << std::endl
			    << "Message: " << ec.message() << std::endl;
        }
    }

    void CommunicationInterfaceLocal::handle_disconnect_as_coordinator(const CommunicationDataPtr& comm_data)
    {
        // deregister agent
        std::unique_lock<std::shared_mutex> guard_coordinator(mutex_coordinator_);

        if(comm_data->agent_info_->id_ >= 0)
            if (coordinator_->deregister_agent(*comm_data->agent_info_))
            {
			    log_->print(DebugType::Message) << "[CommunicationInterfaceLocal::handleDisconnectAsCoordinator]"
				    << " Agent " << comm_data->agent_info_->id_ << " disconnected." << std::endl;
            }

        // delete CommunicationData from vector
        std::unique_lock<std::shared_mutex> guard_comm_data_vec(mutex_comm_data_vec_);
        DataConversion::erase_element_from_vector(comm_data_vec_, comm_data);
    }

    void CommunicationInterfaceLocal::acceptHandler(const std::error_code& ec, CommunicationDataPtr comm_data)
    {
        if(ec)
        {
            // reset flag
            comm_data->is_connected_ = false;

            log_->print(DebugType::Error) << "Failed to accept connection with " << comm_data->communication_info_->agent_type_
                << " with id " << comm_data->communication_info_->id_ << std::endl 
                << ". Error number: " << ec.value() << std::endl 
                << "Error message: " << ec.message() << std::endl;
            return;
        }
    
        // set flag
        comm_data->is_connected_ = true;

        // start sending a periodic ping to check, if the connection is still open
        start_ping(comm_data);

        //it's always an agent that is accepted
        comm_data->communication_info_->agent_type_ = "agent";

        // read data from socket
        async_read_some(comm_data);

        // create new socket
        std::unique_lock<std::shared_mutex> guard(mutex_comm_data_vec_);
        comm_data_vec_.push_back(std::shared_ptr<CommunicationData>(new CommunicationData(asio_->ioService_)));
        guard.unlock();

        // accept connection on new socket
        async_accept(comm_data_vec_.back());
    }

    void CommunicationInterfaceLocal::connectHandler(const std::error_code &ec, CommunicationDataPtr comm_data)
    {
        if(ec)
        {
            // reset flag
            comm_data->is_connected_ = false;

            log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::connectHandler] "
                << "Failed to connect with " << comm_data->communication_info_->agent_type_;
            if( comm_data->communication_info_->agent_type_ == "agent" )
                log_->print(DebugType::Warning) << " with id " << comm_data->communication_info_->id_;
            log_->print(DebugType::Warning) << "." << std::endl << "Error value: " << ec.value() << std::endl
                      << "Error message: " << ec.message() << std::endl;

            // close the socket
            close_socket(comm_data);

            // wait for a second
            std::this_thread::sleep_for(std::chrono::seconds(general_waiting_time_s_));

            // try to reconnect if comm_data is still in comm_data_vec
            async_connect(comm_data, comm_data->communication_info_->ip_, comm_data->communication_info_->port_);

            return;
        }
    
        // set flag
        comm_data->is_connected_ = true;

        // start the period ping to check if connection is still open
        start_ping(comm_data);

        // connected to coordinator
        if (comm_data->communication_info_->agent_type_ == "coordinator")
        {
            log_->print(DebugType::Message) << "[CommunicationInterfaceLocal::connectHandler] "
                << "Connected to coordinator." << std::endl;

            // start polling the registrations
            start_polling(comm_data_vec_.back());

            // accept new connection
            std::lock_guard<std::shared_mutex> guard(mutex_comm_data_vec_);
            comm_data_vec_.push_back(std::shared_ptr<CommunicationData>(new CommunicationData(asio_->ioService_)));
            async_accept(comm_data_vec_.back());
        }
        // connected to agent
        else
        {
			// send communication info
			const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_send_communication_info>(comm_info_local_));
			async_send(comm_data, serialize(message));

            log_->print(DebugType::Message) << "[CommunicationInterfaceLocal::connectHandler] Connected to agent with id "
                << std::to_string(comm_data->communication_info_->id_) << std::endl;
        }

        // read data from socket
        async_read_some(comm_data);
    }

    void CommunicationInterfaceLocal::start_server()
    {
        std::lock_guard<std::shared_mutex> guard(mutex_comm_data_vec_);

        //accept new connection
        comm_data_vec_.push_back(std::shared_ptr<CommunicationData>(new CommunicationData(asio_->ioService_)));
        async_accept(comm_data_vec_.back());

        // start threads
        for( int i = 0; i < number_of_threads_; ++i  )
            asio_->threads_for_communication_.push_back(std::thread([this]() {asio_->ioService_.run(); }));
    }

    void CommunicationInterfaceLocal::start_client()
    {
        std::lock_guard<std::shared_mutex> guard(mutex_comm_data_vec_);

        // create communication data
        comm_data_vec_.push_back(std::shared_ptr<CommunicationData>(new CommunicationData(asio_->ioService_)));

        // connect to coordinator
        async_connect(comm_data_vec_.back(), comm_info_coordinator_.ip_, comm_info_coordinator_.port_);
        comm_info_coordinator_.agent_type_ = "coordinator";

        // safe coordinator communication info in communication data
        comm_data_vec_.back()->communication_info_ = std::make_shared<CommunicationInfo>(comm_info_coordinator_);

        // start threads
        for( int i = 0; i < number_of_threads_; ++i  )
            asio_->threads_for_communication_.push_back(std::thread([this]() {asio_->ioService_.run(); }));
    }

    void CommunicationInterfaceLocal::poll_check_connection(const std::error_code &ec, CommunicationDataPtr comm_data)
    {
        if(ec || !disconnect_at_timeout_)
            return;

        if(comm_data->communication_info_->agent_type_ == "coordinator")
            log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::poll_check_connection] Coordinator is not responding." << std::endl;
        else if(comm_data->communication_info_->agent_type_ == "agent")
            log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::poll_check_connection] Agent with id " << comm_data->communication_info_->id_
                      << " is not responding." << std::endl;    

        if(comm_data->is_connected_) 
            handle_disconnect(comm_data);
    }

    void CommunicationInterfaceLocal::start_ping(const CommunicationDataPtr& comm_data)
    {
        // start polling ping
        comm_data->timer_ping_.expires_from_now(std::chrono::seconds(ping_period_s_));
        comm_data->timer_ping_.async_wait(std::bind(&CommunicationInterfaceLocal::poll_ping,
            this, std::placeholders::_1, comm_data) );
    }

    void CommunicationInterfaceLocal::poll_ping(const std::error_code &ec, CommunicationDataPtr comm_data)
    {
        if(ec)
            return;

        // do not ping, if connection is closed
        if( !comm_data->is_connected_ ) 
            return;

        // start timer to measure delay
        comm_data->timer_check_connection_.expires_from_now(std::chrono::milliseconds(ping_waiting_time_ms_));
        comm_data->timer_check_connection_.async_wait(std::bind(&CommunicationInterfaceLocal::poll_check_connection,
            this, std::placeholders::_1, comm_data) );

		// send ping
		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_get_ping>());
		async_send(comm_data, serialize(message));

        // restart timer for next ping
        start_ping(comm_data);
    }

    void CommunicationInterfaceLocal::start_polling(const CommunicationDataPtr& comm_data) const
    {
        comm_data->timer_polling_.expires_from_now(std::chrono::seconds(polling_period_s_));
        comm_data->timer_polling_.async_wait(std::bind(&CommunicationInterfaceLocal::poll_registration,
            this, std::placeholders::_1, comm_data) );
    }

    void CommunicationInterfaceLocal::poll_registration(const std::error_code &ec, CommunicationDataPtr comm_data) const
    {
        if(ec)
            return;

        std::shared_lock<std::shared_mutex> guard(mutex_infos_);

        // register pending Agent infos
        for (const auto& info : agent_infos_registering_)
        {
			const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_register_agent>(info));
			async_send(comm_data, serialize(message));
        }

        // deregister agent infos
        for (const auto& info : agent_infos_deregistering_)
        {
			const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_deregister_agent>(info));
			async_send(comm_data, serialize(message));
        }

        // register pending Coupling Infos
        for (const auto& info : coupling_infos_pending_)
        {
			const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_register_coupling>(info));
			async_send(comm_data, serialize(message));
        }

        // deregister Coupling Infos
        for (const auto& info : coupling_infos_deregistering_)
        {
			const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_deregister_coupling>(info));
			async_send(comm_data, serialize(message));
        }

        // restart timer
        start_polling(comm_data);
    }

    void CommunicationInterfaceLocal::set_agent(const AgentPtr &agent)
    {
        std::lock_guard<std::shared_mutex> guard(mutex_basics_);
        agent_ = agent;
    }

    void CommunicationInterfaceLocal::set_coordinator(const CoordinatorPtr& coordinator)
    {
        std::lock_guard<std::shared_mutex> guard(mutex_coordinator_);
        coordinator_ = coordinator;
    }

    void CommunicationInterfaceLocal::set_communicationInfo_for_coordinator(const CommunicationInfo& info)
    {
        std::lock_guard<std::shared_mutex> guard(mutex_basics_);
        comm_info_coordinator_ = info;
    }

    void CommunicationInterfaceLocal::set_simulator(const SimulatorPtr& simulator)
    {
        std::lock_guard<std::shared_mutex> guard(mutex_basics_);
        simulator_ = simulator;
    }

    const bool CommunicationInterfaceLocal::register_agent(const AgentPtr& agent)
    {
        std::lock_guard<std::shared_mutex> guard_basics(mutex_basics_);
        std::unique_lock<std::shared_mutex> guard_infos(mutex_infos_);

        comm_info_local_.id_ = agent->get_id();
        agent_ = agent;
        agent_infos_registering_.push_back( agent->get_agentInfo() );

        return true;
    }

    const bool CommunicationInterfaceLocal::deregister_agent(const AgentInfo& agent)
    {
	    std::unique_lock<std::shared_mutex> guard_infos(mutex_infos_);
        agent_infos_deregistering_.push_back(agent);

        return true;
    }

    const bool CommunicationInterfaceLocal::register_coupling(const CouplingInfo& coupling)
    {
        std::lock_guard<std::shared_mutex> guard(mutex_infos_);

        // send coordinator message to register coupling
        coupling_infos_pending_.push_back(coupling);

        // do not block coupling anymore
        DataConversion::erase_element_from_vector(coupling_infos_blocked_, coupling);

        return true;
    }

    const bool CommunicationInterfaceLocal::deregister_coupling(const CouplingInfo& coupling)
    {
        std::lock_guard<std::shared_mutex> guard(mutex_infos_);

        // send coordinator message to deregister coupling
        coupling_infos_deregistering_.push_back(coupling);

        // block the coupling from re-registering
        coupling_infos_blocked_.push_back(coupling);

        return true;
    }

    const bool CommunicationInterfaceLocal::send_numberOfNeighbors(const int number, const int from, const int to)
    {
        // get CommunicationData
        const auto comm_data = get_communicationData(to);

        // check if agent is known
        if(comm_data == nullptr)
        {
            log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::sendNumberOfNeighbors] "
                << "Could not find agent with id " << to << "." << std::endl;
            return false;
        }

		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_send_numberOfNeighbors>(number, from));
		async_send(comm_data, serialize(message));

        return true;
    }

    const bool CommunicationInterfaceLocal::send_localCopies(const AgentState& state, int from, int to)
    {
        // get CommunicationData
        const auto comm_data = get_communicationData(to);

        // check if agent is known
        if(comm_data == nullptr)
        {
            log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::sendAgentState] "
                << "Could not find agent with id " << to << "." << std::endl;
            return false;
        }

		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_send_local_copies>(state, from));
		async_send(comm_data, serialize(message));

        return true;
    }

    const bool CommunicationInterfaceLocal::send_agentState(const AgentState& state, const ConstraintState& constr_state, int from, int to)
    {
        // get CommunicationData
        const auto comm_data = get_communicationData(to);

        // check if agent is known
        if (comm_data == nullptr)
        {
            log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::send_agentState] "
                << "Could not find agent with id " << to << "." << std::endl;
            return false;
        }

        const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_send_agent_state>(state,constr_state, from));
        async_send(comm_data, serialize(message));

        return true;
    }

    const bool CommunicationInterfaceLocal::send_desiredAgentState(const AgentState &desired_state, int from, int to)
    {
        // get CommunicationData
        const auto comm_data = get_communicationData(to);

        // check if agent is known
        if(comm_data == nullptr)
        {
            log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::sendDesiredState] "
                << "Could not find agent with id " << to << "." << std::endl;
            return false;
        }

		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_send_desired_agent_state>(desired_state, from));
		async_send(comm_data, serialize(message));

        return true;
    }

    const std::shared_ptr< AgentState > CommunicationInterfaceLocal::get_agentState_from_agent(int agentId) const
    {
        const auto comm_data = get_communicationData(agentId);

        if (comm_data == nullptr)
        {
            log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::get_agentStateFromAgent] Agent with id " << agentId
                << " not found." << std::endl;
            return nullptr;
        }

	    // send request for true agent state
	    std::unique_lock<std::mutex> guard(comm_data->mutex_agentState_for_simulation_);
	    comm_data->flag_agentState_for_simulation_ = false;
		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_get_agentState_for_simulation>());
		async_send(comm_data, serialize(message));

	    // wait for response
	    comm_data->cond_var_agentState_for_simulation_.wait_for(guard, std::chrono::seconds(general_waiting_time_s_),
		    [comm_data] {return comm_data->flag_agentState_for_simulation_; });

        // check if agent did respond
        if (!comm_data->flag_agentState_for_simulation_)
            log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::get_agentStateFromAgent] "
            << "Agent did not respond. I return last known state." << std::endl;

        return comm_data->agentState_for_simulation_;
    }

    const std::shared_ptr< AgentState > CommunicationInterfaceLocal::get_desiredAgentState_from_agent(int agentId) const
    {
        const auto comm_data = get_communicationData(agentId);

        if( comm_data == nullptr )
        {
            log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::get_desiredAgentStateFromAgent] Agent with id "
                << agentId << " not found." << std::endl;
            return nullptr;
        }

        // send request for desired agentState
        std::unique_lock<std::mutex> guard(comm_data->mutex_desired_agentState_);
        comm_data->flag_desired_agentState_ = false;

		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_get_desired_agent_state_from_agent>());
		async_send(comm_data, serialize(message));

        // wait for response
        comm_data->cond_var_desired_agentState_.wait_for(guard, std::chrono::seconds(general_waiting_time_s_),
            [comm_data] {return comm_data->flag_desired_agentState_; });

        // check if agent did respond
        if (!comm_data->flag_desired_agentState_)
            log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::get_desiredAgentStateFromAgent] Agent with id " << agentId
            << " did not respond. I use last known state." << std::endl;

        return comm_data->desired_agentState_;

    }

    const std::shared_ptr< std::map<unsigned int, AgentInfoPtr > > CommunicationInterfaceLocal::get_agentInfo_from_coordinator() const
    {
        std::shared_lock<std::shared_mutex> guard(mutex_coordinator_);
        const auto& agentInfos = coordinator_->get_agentInfos();

        // check if map with agents should be adapted
        bool adapt_map = false;
	    for (const auto& [id, info] : agentInfos)
	    {
		    const auto comm_data = get_communicationData(id);

		    // do not consider agents, that are not known
            if (comm_data == nullptr)
                adapt_map = true;
		    // do not consider agents that are not configured
            else if (!comm_data->is_configured)
                adapt_map = true;
	    }

        // if the map is ok, return it
        if(!adapt_map)
            return std::make_shared< std::map< unsigned int, AgentInfoPtr > >(agentInfos);

        // otherwise safe a local copy and adapt it
        auto localCopy_agentInfos = agentInfos;
        for (auto iterator = localCopy_agentInfos.begin(); iterator != localCopy_agentInfos.end(); )
        {
            const auto id = iterator->first;

            const auto comm_data = get_communicationData(id);

            // do not consider agents, that are not known
            if (comm_data == nullptr)
                iterator = localCopy_agentInfos.erase(iterator);
            // do not consider agents that are not configured
            else if (!comm_data->is_configured)
                iterator = localCopy_agentInfos.erase(iterator);
            else
                ++iterator;
        }

        return std::make_shared< std::map< unsigned int, AgentInfoPtr > >(localCopy_agentInfos);
    }

    const std::shared_ptr< std::map< unsigned int, std::vector< CouplingInfoPtr > > > CommunicationInterfaceLocal::get_sendingNeighbors_from_coordinator() const
    {
        std::shared_lock<std::shared_mutex> guard(mutex_coordinator_);
        return std::make_shared< std::map< unsigned int, std::vector< CouplingInfoPtr > > > (coordinator_->get_sendingNeighbors());
    }

    const std::shared_ptr< std::map< unsigned int, std::vector< CouplingInfoPtr > > > CommunicationInterfaceLocal::get_receivingNeighbors_from_coordinator() const
    {
        std::shared_lock<std::shared_mutex> guard(mutex_coordinator_);
        return std::make_shared< std::map< unsigned int, std::vector< CouplingInfoPtr > > > (coordinator_->get_receivingNeighbors());
    }

    const AgentModelPtr CommunicationInterfaceLocal::get_agentModel(int agentId) const
    {
        // get corresponding communication data
        const auto comm_data = get_communicationData(agentId);

        // check if communication data is found
        if (comm_data == nullptr)
        {
            log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::get_agentModel] "
                << "Agent with id " << agentId << " not found." << std::endl;
		    return nullptr;
        }

        std::unique_lock<std::mutex> guard(comm_data->mutex_agentModel_);
        comm_data->flag_agentModel_ = false;

	    // send request for agent model
		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_get_agent_model_from_agent>());
		async_send(comm_data, serialize(message));

        // wait for response
        comm_data->cond_var_agentModel_.wait_for(guard, std::chrono::seconds(general_waiting_time_s_),
            [comm_data] { return comm_data->flag_agentModel_; });

        // check if agent did respond
	    if (!comm_data->flag_agentModel_)
            log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::get_agentModel] Agent with id " << agentId
            << " did not respond. I return last known model." << std::endl;

        return comm_data->agentModel_for_simulation_;
    }

    const std::shared_ptr< std::map<int, CouplingModelPtr> > CommunicationInterfaceLocal::get_couplingModels_from_agent( int agentId ) const
    {
        // get corresponding communication data
        const auto comm_data = get_communicationData(agentId);

        // check if communication data is found
        if (comm_data == nullptr)
        {
            log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::get_couplingModelsFromAgent] "
                << "Agent with id " << agentId << " not found." << std::endl;
            return nullptr;
        }

        std::unique_lock<std::mutex> guard(comm_data->mutex_couplingModels_);
        comm_data->flag_couplingModels_ = false;

        // send request for coupling models
		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_get_coupling_model_from_agent>());
		async_send(comm_data, serialize(message));

        // wait for response
        comm_data->cond_var_couplingModels_.wait_for(guard, std::chrono::seconds(general_waiting_time_s_),
            [comm_data] { return comm_data->flag_couplingModels_; });

        // check if agent did respond
	    if (!comm_data->flag_couplingModels_)
            log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::get_couplingModelsFromAgent] "
            << "Agent with id " << agentId << " did not respond." << std::endl;

        return comm_data->couplingModel_for_simulation_;
    }

    const bool CommunicationInterfaceLocal::send_couplingState(const CouplingState& state, int from, int to)
    {
        // get CommunicationData
        const auto comm_data = get_communicationData(to);

        // check if agent is known
        if(comm_data == nullptr)
        {
            log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::sendCouplingState] Agent with id " << to
                << " not found." << std::endl;

            return false;
        }

        // send coupling states
		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_send_coupling_state>(state, from));
		async_send(comm_data, serialize(message));

        return true;
    }

    const bool CommunicationInterfaceLocal::send_couplingState(const CouplingState& state, const CouplingState& state2, int from, int to)
    {
        // get CommunicationData
        const auto comm_data = get_communicationData(to);

        // check if agent is known
        if(comm_data == nullptr)
        {
            log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::sendCouplingState] Agent with id " << to
                << " not found." << std::endl;

            return false;
        }

        // send coupling states
		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_send_two_coupling_states>(state, state2, from));
		async_send(comm_data, serialize(message));

        return true;
    }

    const bool CommunicationInterfaceLocal::send_multiplierState(const MultiplierState& multiplier, PenaltyState penalty, int from, int to)
    {
        // get CommunicationData
        const auto comm_data = get_communicationData(to);

        // check if agent is known
        if(comm_data == nullptr)
        {
            log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::sendMultiplierState] Agent with id " << to
                << " not found." << std::endl;

            return false;
        }

        // send multiplier states
		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_send_multiplier_state>(multiplier, penalty, from));
		async_send(comm_data, serialize(message));

        return true;
    }

    const bool CommunicationInterfaceLocal::send_convergenceFlag(bool converged, int from)
    {
        // get CommunicationData
        const auto comm_data = get_communicationData("coordinator");

        // check if agent is known
        if (comm_data == nullptr)
        {
            log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::sendConvergenceFlag] "
                << "Coordinator not found." << std::endl;

            return false;
        }

        // send convergence flag
		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_send_convergenceFlag>(converged, from));
		async_send(comm_data, serialize(message));

        return true;
    }

    const bool  CommunicationInterfaceLocal::send_stoppedAlgFlag(const bool flag, const int from, const int to)
    {
         // get communication data 
        const auto comm_data = get_communicationData(to);

        // check if agent is known
        if (comm_data == nullptr)
        {
            log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::sendflagStoppedAdmm] Agent with id " << to
                << " not found." << std::endl;

            return false;
        }

        // send admm stop flag to agent
        const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_send_flag_stopped_admm>(flag,from));
        async_send(comm_data, serialize(message));

        return true;
    }
 
    const bool CommunicationInterfaceLocal::send_stoppedAlgFlag(const bool flag, const int from)
    {
        // get CommunicationData
        const auto comm_data = get_communicationData("coordinator");

        // check if agent is known
        if (comm_data == nullptr)
        {
            log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::sendConvergenceFlag] "
                << "Coordinator not found." << std::endl;

            return false;
        }

        // send admm stop flag to coordinator 
        const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_send_flag_stopped_admm_coordinator>(flag, from));
        async_send(comm_data, serialize(message));

        return true;

    }

    const bool CommunicationInterfaceLocal::send_flagToStopAdmm( const bool flag, const int to)
    {

        // get CommunicationData
        const auto comm_data = get_communicationData(to);

        // check if agent is known
        if (comm_data == nullptr)
        {
            log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::send_flagToStopAdmm]  Agent with id " << to
                << " not found." << std::endl;

            return false;
        }

        // send flag to stop admm
        const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_send_flag_to_stop_admm>(flag));
        async_send(comm_data, serialize(message));

        return true;

    }

    const bool CommunicationInterfaceLocal::configure_optimization(const OptimizationInfo& info)
    {
	    std::shared_lock<std::shared_mutex> guard(mutex_comm_data_vec_);
        std::unique_lock<std::mutex> guard_config_optim(mutex_config_optimizationInfo_);
        std::vector<CommunicationDataPtr> comm_data_to_wait_for;
        numberOfNotifications_config_optimizationInfo_ = 0;

        // configure optimization for each agent
		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_send_optimizationInfo>(info));
		const auto data = serialize(message);
        for(const auto& comm_data : comm_data_vec_)
        {
            if( comm_data->is_connected_ )
            {
                comm_data_to_wait_for.push_back(comm_data);
                ++numberOfNotifications_config_optimizationInfo_;
                async_send(comm_data, data);
            }
        }

        guard.unlock();

        // wait for response
        conditionVariable_config_optimizationInfo_.wait_for(guard_config_optim, std::chrono::seconds(general_waiting_time_s_),
            [this] { return numberOfNotifications_config_optimizationInfo_ == 0; });

        // check if every agent did configure in time
        if (numberOfNotifications_config_optimizationInfo_ == 0)
            return true;

        // check which agent did not configure in time
        for (const auto& comm_data : comm_data_to_wait_for)
        {
            if (!comm_data->is_configured)
                log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::configureOptimization] Agent "
                << comm_data->communication_info_->id_
                << " did not configure optimization in time." << std::endl;
	    }

	    return false;
    }

    void CommunicationInterfaceLocal::configureOptimization(const CommunicationDataPtr& comm_data)
    {
        if (!comm_data->is_connected_ || comm_data->communication_info_->id_ < 0 
            || comm_data->is_configured || comm_data->is_configuring)
            return;

        // set flag
        comm_data->is_configuring = true;

        // send request to configure
	    std::unique_lock<std::mutex> guard_config_optim(mutex_config_optimizationInfo_); 
		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_send_optimizationInfo>(coordinator_->get_optimizationInfo()));
		async_send(comm_data, serialize(message));

        // wait for response
	    conditionVariable_config_optimizationInfo_.wait_for(guard_config_optim, std::chrono::seconds(general_waiting_time_s_),
		    [comm_data] { return comm_data->is_configured; });

        // reset flag
        comm_data->is_configuring = false;

        // check if agent did configure
        if (!comm_data->is_configured)
        {
            log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::configureOptimization] Agent "
                << comm_data->communication_info_->id_ << " did not configure optimization in time." << std::endl; 
            return;
        }

        log_->print(DebugType::Message) << "[CommunicationInterfaceLocal::configureOptimization] Configured agent "
            << comm_data->communication_info_->id_ << " during runtime." << std::endl;    
    }

    void CommunicationInterfaceLocal::waitFor_connection(const int agents, const int couplings)
    {
        unsigned int number_of_agents = 0;
        unsigned int number_of_couplings = 0;

        while( (number_of_agents != agents) || (number_of_couplings != couplings) )
        {
            number_of_agents = 0;
            number_of_couplings = 0;

            std::shared_lock<std::shared_mutex> guard(mutex_comm_data_vec_);
            std::unique_lock<std::mutex> guard_waitForConnection(mutex_waitForConnection_);
            numberOfNotifications_waitForConnection_ = 0;

            // ask all agents for number of connected agents
			const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_get_number_of_active_couplings>());
			const auto data = serialize(message);
            std::vector<CommunicationDataPtr> comm_data_to_wait_for;
            for (const auto& comm_data : comm_data_vec_)
            {
                if( comm_data->is_connected_ )
                {
                    comm_data_to_wait_for.push_back(comm_data);
                    ++numberOfNotifications_waitForConnection_;
                    async_send(comm_data, data);
                }
            }

            guard.unlock();

            // wait for response
            conditionVariable_waitForConnection_.wait_for(guard_waitForConnection, std::chrono::seconds(general_waiting_time_s_),
                [this] { return numberOfNotifications_waitForConnection_ == 0; });

            // check if agents did respond
            if (numberOfNotifications_waitForConnection_ != 0)
                log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::waitFor_connection] "
                << "Some agents did not respond." << std::endl;

            std::shared_lock<std::shared_mutex> guard_comm_data_vec(mutex_comm_data_vec_);

            // count number of answers and number of couplings
            for (const auto& comm_data : comm_data_to_wait_for)
            {
			    ++number_of_agents;
			    number_of_couplings += comm_data->number_of_connected_agents_;
            }

		    // wait some time until next check
		    std::this_thread::sleep_for(std::chrono::seconds(general_waiting_time_s_));
        }
    }

    const bool CommunicationInterfaceLocal::trigger_step(const AlgStep& step)
    {
        // prepare protocol
		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_trigger_step>(step));
		const auto data = serialize(message);

        // send trigger to each agent
        std::shared_lock<std::shared_mutex> guard_commDataVec(mutex_comm_data_vec_);
        std::unique_lock<std::mutex> guard_triggerStep(mutex_triggerStep_);
        numberOfNotifications_triggerStep_ = 0;

        for(const auto& comm_data : comm_data_vec_)
        {
            if (comm_data->is_connected_ && !comm_data->is_configured)
                asio::post([this, comm_data]() { configureOptimization(comm_data); });
            else if( comm_data->is_connected_ && comm_data->is_configured)
		    {
			    ++numberOfNotifications_triggerStep_;
                async_send(comm_data, data);
            }
        }
        guard_commDataVec.unlock();

        // wait for response
        conditionVariable_triggerStep_.wait_for(guard_triggerStep, std::chrono::seconds(general_waiting_time_s_),
            [this] {return numberOfNotifications_triggerStep_ == 0; });

        if (numberOfNotifications_triggerStep_ != 0)
            log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::triggerStep] "
            << "Some agents did not execute triggered step in time." << std::endl;

        return true;
    }

    void CommunicationInterfaceLocal::trigger_simulation(const std::string& Integrator, typeRNum dt)
    {
        std::shared_lock<std::shared_mutex> guard(mutex_basics_);
        simulator_->distributed_simulation(Integrator, dt);
    }

    void CommunicationInterfaceLocal::set_simulatedState_for_agent
    (
        int agentId, 
        const std::vector<typeRNum>& new_state, 
        typeRNum dt, 
		typeRNum t0,
		const typeRNum cost
    )
    {
        const auto comm_data = get_communicationData(agentId);

        std::unique_lock<std::mutex> guard_triggerStep(mutex_triggerStep_);
        numberOfNotifications_triggerStep_ = 0;
        ++numberOfNotifications_triggerStep_;

		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_send_simulated_state>(new_state, dt, t0, cost));
		async_send(comm_data, serialize(message));

        // wait for response
        conditionVariable_triggerStep_.wait_for(guard_triggerStep, std::chrono::seconds(general_waiting_time_s_),
            [this] {return numberOfNotifications_triggerStep_ == 0; });

        if (numberOfNotifications_triggerStep_ != 0)
            log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::set_simulatedState_for_agent] "
            << "Some agents did not execute triggered step "
            << "in time." << std::endl;
    }

    void CommunicationInterfaceLocal::writeHandler(const std::error_code& ec,
        std::size_t bytes_transferred, CommunicationDataPtr comm_data,
        std::shared_ptr<std::vector<char>> data_ptr) const
    {
        if (ec)
        {
            log_->print(DebugType::Error) << "[CommunicationInterfaceLocal::writeHandler] "
                << "Failed to send data." << std::endl
                << "Error number: " << ec.value() << std::endl
                << "Error message: " << ec.message() << std::endl;

            return;
        }

        // delete corresponding buffer
        std::unique_lock<std::mutex> guard(comm_data->mutex_buffer_send_);
        for (unsigned int i = 0; i < comm_data->buffer_send_.size(); ++i)
        {
            if (comm_data->buffer_send_[i] == data_ptr)
                comm_data->buffer_send_.erase(comm_data->buffer_send_.begin() + i);
        }
    }

    void CommunicationInterfaceLocal::readHandler(const std::error_code& ec,
        const size_t& amountOfBytes, CommunicationDataPtr comm_data)
    {
        if (ec)
        {
            // handle the recognized disconnect
            if(comm_data->is_connected_) 
                handle_disconnect(comm_data);
            return;
        }

        // safe data
        const unsigned int pos_in_vector = static_cast<unsigned int>( comm_data->data_.size() );
        comm_data->data_.resize(comm_data->data_.size() + amountOfBytes, 0);
        std::copy(comm_data->buffer_read_.begin(), comm_data->buffer_read_.begin() + amountOfBytes, comm_data->data_.begin() + pos_in_vector);

        // only evaluate the data, if at least 4 bytes are read
        while(comm_data->data_.size() >= 4)
        {
            unsigned int size_of_packet = 0;
			unsigned int pos = 0;

            // read size of packet
            std::from_chars(&comm_data->data_.at(0), &comm_data->data_.at(Constants::SIZE_OF_HEADERS_), size_of_packet);

            // if enough data is read, evaluate it
            if (comm_data->data_.size() >= size_of_packet)
			{
                // create local copy of data
				auto data = std::make_shared<std::stringstream>();
				for (unsigned int k = Constants::SIZE_OF_HEADERS_; k < size_of_packet; ++k)
					(*data) << comm_data->data_[k];

				comm_data->data_.erase(comm_data->data_.begin(), comm_data->data_.begin() + size_of_packet);

                // evaluate data asynchronously
                asio::post([comm_data, this, data]() { evaluate_data(comm_data, data); });
            }
            else
                break;
        }

        // read data from socket
        if(comm_data->is_connected_) 
            async_read_some(comm_data);
    }

    void CommunicationInterfaceLocal::connect_to_neighbor(const CommunicationInfoPtr& info)
    {
        // check if it is actually a neighbor
        if( info->id_ == comm_info_local_.id_ ) 
            return;

        std::shared_lock<std::shared_mutex> soft_guard_comm_data_vec(mutex_comm_data_vec_);

        // check if connection is already established
        for( const auto& comm_data : comm_data_vec_ )
            if (comm_data->communication_info_->id_ == info->id_) 
                return;

        soft_guard_comm_data_vec.unlock();

        // it's always the agent with the higher id that is connecting to
        // the agent with lower id to prevent race conditions
        if( info->id_ > comm_info_local_.id_ ) 
            return;

        // create new connection
        std::unique_lock<std::shared_mutex> strong_guard_comm_data_vec(mutex_comm_data_vec_);
        comm_data_vec_.push_back(std::shared_ptr<CommunicationData>( new CommunicationData(asio_->ioService_)));
        strong_guard_comm_data_vec.unlock();

        // add communicationInfo to list
        comm_data_vec_.back()->communication_info_ = info;

        // connect to neighbor
        async_connect(comm_data_vec_.back(), info->ip_, info->port_);
    }

    const bool CommunicationInterfaceLocal::fromCommunication_deregistered_coupling(const CouplingInfo& coupling)
    {
		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_fromCommunication_deregistered_coupling>(coupling));
		const auto data = serialize(message);

        // send to agent
        const auto comm_data_agent = get_communicationData(coupling.agent_id_);
        if(comm_data_agent != nullptr )
            async_send(comm_data_agent, data);

        // send to neighbor
        const auto comm_data_neighbor = get_communicationData(coupling.neighbor_id_);
        if(comm_data_neighbor != nullptr )
            async_send(comm_data_neighbor, data);

        return true;
    }

    const CommunicationDataPtr CommunicationInterfaceLocal::get_communicationData(const int agent_id) const
    {
        std::shared_lock<std::shared_mutex> guard(mutex_comm_data_vec_);

        for(const auto& comm_data : comm_data_vec_)
        {
            if( comm_data->communication_info_->id_ == agent_id ) 
                return comm_data; 

            if( comm_data->agent_info_->id_ == agent_id )
                return comm_data;
        }
     
        return nullptr;
    }

    const CommunicationDataPtr CommunicationInterfaceLocal::get_communicationData(const std::string& agent_type) const
    {
        std::shared_lock<std::shared_mutex> guard(mutex_comm_data_vec_);

        for (const auto& comm_data : comm_data_vec_)
            if (comm_data->communication_info_->agent_type_ == agent_type)
                return comm_data;

        return nullptr;
    }

    void CommunicationInterfaceLocal::waitFor_flag_from_coordinator()
    {
        const auto comm_data = get_communicationData("coordinator");

        // check if agent is known
        if (comm_data == nullptr)
        { 
            log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::getFlagFromCoordinator] "
                << "Coordinator not found. First set coordinator." << std::endl;

            return;
        }

        // wait until flag is set from coordinator
        comm_data->flag_from_coordinator_ = false;
	    while (!comm_data->flag_from_coordinator_)
		    std::this_thread::sleep_for(std::chrono::seconds(general_waiting_time_s_));
    }

    void CommunicationInterfaceLocal::cap_stored_data(unsigned int data_points)
    {
        std::unique_lock<std::shared_mutex> guard_basic(mutex_basics_);
        agent_->get_solution()->maximum_number_of_data_points_ = data_points;
    }

    void CommunicationInterfaceLocal::set_passive()
    {
	    while (true)
		    std::this_thread::sleep_for(std::chrono::seconds(general_waiting_time_s_));
    }

    const SolutionPtr CommunicationInterfaceLocal::get_solution(unsigned int agent_id) const
    {
        if (agent_->get_id() != agent_id)
        {
            log_->print(DebugType::Error) << "[CommunicationInterfaceLocal::get_solution] "
                << "Unknown agent." << std::endl;
            return nullptr;
        }

        return agent_->get_solution();
    }

    const std::vector<SolutionPtr> CommunicationInterfaceLocal::get_solution(const std::string& agents) const
    {
        if (agents != "all")
        {
            log_->print(DebugType::Error) << "[CommunicationInterfaceLocal::get_solution] "
                << "Unknown set of agents." << std::endl;
            return std::vector<SolutionPtr>();
        }

        if (!coordinator_)
        {
            log_->print(DebugType::Error) << "[CommunicationInterfaceLocal::get_solution] "
                << "Coordinator is required to sample all solutions." << std::endl;

		    return std::vector<SolutionPtr>();
        }

        std::vector<CommunicationDataPtr> comm_data_to_wait_for;

        std::shared_lock<std::shared_mutex> guard(mutex_comm_data_vec_);
        std::unique_lock<std::mutex> guard_getSolutions(mutex_getSolutions_);
        numberOfNotifications_getSolutions_ = 0;

        // send request for solution
		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_get_solution>());
		const auto data = serialize(message);
        for (const auto& comm_data : comm_data_vec_)
        {
            if (comm_data->is_connected_)
            {
                ++numberOfNotifications_getSolutions_;
                comm_data_to_wait_for.push_back(comm_data);
                async_send(comm_data, data);
            }
        }

        guard.unlock();

        // wait for response
        conditionVariable_getSolutions_.wait_for(guard_getSolutions, std::chrono::seconds(general_waiting_time_s_),
            [this] {return numberOfNotifications_getSolutions_ == 0; });

        // check if agents did respond
        if (numberOfNotifications_getSolutions_ != 0)
            log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::getSolutions] "
            << "Some agents did not return solution in time." << std::endl;

        // collect solutions
        std::vector<SolutionPtr> solutions(comm_data_to_wait_for.size());
        for (unsigned int i = 0; i < comm_data_to_wait_for.size(); ++i)
            solutions[i] = comm_data_to_wait_for[i]->solution_;

        return solutions;
    }

    void CommunicationInterfaceLocal::reset_solution(unsigned int agent_id)
    {
        if (agent_->get_id() == agent_id)
            agent_->reset_solution();
        else
            log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::reset_solution] "
            << "Agent " << agent_id << " not found." << std::endl;
    }

    void CommunicationInterfaceLocal::reset_solution(const std::string& agents)
    {
        if (agents == "all")
            agent_->reset_solution();
        else
            log_->print(DebugType::Error) << "[CommunicationInterfaceLocal::reset_solution] "
            << "Set of agents not known." << std::endl;
    }

    void CommunicationInterfaceLocal::send_flag_to_agents(const std::string& agents) const
    {
        if (agents == "all")
        {
            std::shared_lock<std::shared_mutex> guard(mutex_comm_data_vec_);
			const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_send_flag_to_agents>());
			const auto data = serialize(message);

            for (const auto& comm_data : comm_data_vec_)
                async_send(comm_data, data);

            return;
        }

        log_->print(DebugType::Error) << "[CommunicationInterfaceLocal::send_flag_to_agents] "
            << "Unknown set of agents '" << agents << "'." << std::endl;
    }

    void CommunicationInterfaceLocal::send_flag_to_agents(const int agent_id) const
    {
        const auto comm_data = get_communicationData(agent_id);

        // check if agent is known
        if (comm_data == nullptr)
        {
            log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::send_flag_to_agents] "
                << "Agent with id " << agent_id << " not found." << std::endl;
            return;
        }

		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_send_flag_to_agents>());
		async_send(comm_data, serialize(message));
    }

    void CommunicationInterfaceLocal::send_flag_to_agents(const std::vector<int>& agent_ids) const
    {
        // prepare data
		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_send_flag_to_agents>());
		const auto data = serialize(message);

        // send flag to each agent
        for (unsigned int k = 0; k < agent_ids.size(); ++k)
        {
            const auto comm_data = get_communicationData(agent_ids[k]);

            // check if agent is known
            if (comm_data == nullptr)
                log_->print(DebugType::Warning) << "[CommunicationInterfaceLocal::send_flag_to_agents] "
                << "Agent with id " << agent_ids[k] << " not found." << std::endl;
            else
                async_send(comm_data, data);
        }      
    }

    void CommunicationInterfaceLocal::async_send(const CommunicationDataPtr& comm_data, const std::shared_ptr< std::vector<char> >& data) const
    {
        if (comm_data == nullptr)
            return;

        if (!comm_data->is_connected_)
            return;

        // as everything is ok, send the data
        std::unique_lock<std::mutex> guard(comm_data->mutex_buffer_send_);

        comm_data->buffer_send_.push_back(data);
        comm_data->socket_.async_write_some(asio::buffer(*comm_data->buffer_send_.back()),
            comm_data->strand_.wrap(std::bind(&CommunicationInterfaceLocal::writeHandler,
                this, std::placeholders::_1, std::placeholders::_2,
                comm_data, comm_data->buffer_send_.back())));
    }

    void CommunicationInterfaceLocal::async_connect(const CommunicationDataPtr& comm_data)
    {
        if (comm_data == nullptr)
            return;

        std::shared_lock<std::shared_mutex> guard_socket(comm_data->mutex_socket_);

        comm_data->socket_.async_connect(comm_data->socket_.remote_endpoint(),
            comm_data->strand_.wrap(std::bind(&CommunicationInterfaceLocal::connectHandler,
                this, std::placeholders::_1, comm_data)));
    }

    void CommunicationInterfaceLocal::async_connect(const CommunicationDataPtr& comm_data, const std::string& ip, const std::string& port)
    {
        if (comm_data == nullptr)
            return;

        if (comm_data->is_connected_)
            return;

        // do not connect if socket is already open
        std::shared_lock<std::shared_mutex> guard_mutex(comm_data->mutex_socket_);
        if (comm_data->socket_.is_open())
            return;

        // as everything is ok, connect 
        const asio::ip::tcp::endpoint endpoint(asio::ip::address::from_string(ip), std::stoi(port));

        comm_data->socket_.async_connect(endpoint,
            comm_data->strand_.wrap(std::bind(&CommunicationInterfaceLocal::connectHandler,
                this, std::placeholders::_1, comm_data ) ) );
    }

    void CommunicationInterfaceLocal::async_accept(const CommunicationDataPtr& comm_data)
    {
        if (comm_data == nullptr)
            return;

        std::unique_lock<std::shared_mutex> guard_socket(comm_data->mutex_socket_);
        if (comm_data->is_connected_ || comm_data->socket_.is_open())
            return;
    
        // as everything is ok, accept the connection
        asio_->acceptor_.async_accept(comm_data->socket_,
            comm_data->strand_.wrap(std::bind(&CommunicationInterfaceLocal::acceptHandler,
                this, std::placeholders::_1, comm_data ) ) );
    }

    void CommunicationInterfaceLocal::async_read_some(const CommunicationDataPtr& comm_data)
    {
        if (comm_data == nullptr)
            return;

        std::shared_lock<std::shared_mutex> guard_socket(comm_data->mutex_socket_);

        comm_data->socket_.async_read_some(asio::buffer(comm_data->buffer_read_),
            comm_data->strand_.wrap(std::bind(&CommunicationInterfaceLocal::readHandler,
                this, std::placeholders::_1, std::placeholders::_2, comm_data)));
    }

    void CommunicationInterfaceLocal::fromCommunication_send_communicationInfo(const CommunicationDataPtr& comm_data, const CommunicationInfoPtr& info)
    {
        log_->print(DebugType::Message) << "[CommunicationInterfaceLocal::fromCommunication_send_communicationInfo] "
            << "Accepted connection to agent with id " << info->id_ << "." << std::endl;

        comm_data->communication_info_ = info;

        std::unique_lock<std::shared_mutex> guard(mutex_comm_data_vec_);
        if(!DataConversion::is_element_in_vector(comm_data_vec_, comm_data))
            comm_data_vec_.push_back(comm_data);
    }

    void CommunicationInterfaceLocal::fromCommunication_send_flagToAgents(const CommunicationDataPtr& comm_data) const
    {
        comm_data->flag_from_coordinator_ = true;        
    }

    void CommunicationInterfaceLocal::fromCommunication_get_ping(const CommunicationDataPtr& comm_data) const
	{
		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_send_ping>());
		async_send(comm_data, serialize(message));
    }

    void CommunicationInterfaceLocal::fromCommunication_send_ping(const CommunicationDataPtr& comm_data) const
    {
        // reset timer
        comm_data->timer_check_connection_.cancel();
    }

    void CommunicationInterfaceLocal::fromCommunication_get_numberOfActiveCouplings(const CommunicationDataPtr& comm_data) const
    {
        std::shared_lock<std::shared_mutex> guard(mutex_infos_);
        int number = 0;

        // count number of couplings with an open connection
        for(const auto& info : coupling_infos_active_)
        {
            const auto comm_data_search = get_communicationData(info.neighbor_id_);

            if (comm_data_search != nullptr)
                if (comm_data_search->is_connected_)
                    ++number;
        }
		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_send_number_of_active_couplings>(number));
		async_send(comm_data, serialize(message));
    }

    void CommunicationInterfaceLocal::fromCommunication_send_numberOfActiveCouplings(const CommunicationDataPtr& comm_data, int number) const
    {
        std::unique_lock<std::mutex> guard(mutex_waitForConnection_);
        --numberOfNotifications_waitForConnection_;
        comm_data->number_of_connected_agents_ = number;
        conditionVariable_waitForConnection_.notify_one();
    }

    void CommunicationInterfaceLocal::fromCommunication_send_couplingState(const CommunicationDataPtr& comm_data, const CouplingState& state, int from)
    {
        std::unique_lock<std::shared_mutex> guard(mutex_basics_);
        agent_->fromCommunication_received_couplingState(state, from);
    }

    void CommunicationInterfaceLocal::fromCommunication_send_couplingState(const CommunicationDataPtr& comm_data, 
        const CouplingState& state1, const CouplingState& state2, int from)
    {
        std::unique_lock<std::shared_mutex> guard(mutex_basics_);
        agent_->fromCommunication_received_couplingState(state1, state2, from);
    }

    void CommunicationInterfaceLocal::fromCommunication_get_desiredAgentStateFromAgent(const CommunicationDataPtr& comm_data) const
    {
		// if mutex_basics_ is used here, there is a deadlock if neighbor approximation is used
		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_send_desired_agent_state>(agent_->get_desiredAgentState(), agent_->get_id()));
		async_send(comm_data, serialize(message));
    }

    void CommunicationInterfaceLocal::fromCommunication_send_simulatedState
    (
        const CommunicationDataPtr& comm_data, 
        const std::vector<typeRNum>& x_next, 
        typeRNum dt, 
        typeRNum t0,
        typeRNum cost
    )
    {
        std::unique_lock<std::shared_mutex> guard(mutex_basics_);
        agent_->set_updatedState(x_next, dt, t0, cost);
		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_acknowledge_executed_ADMM_step>());
		async_send(get_communicationData("coordinator"), serialize(message));
    }

    void CommunicationInterfaceLocal::fromCommunication_send_couplingModel(const CommunicationDataPtr& comm_data, const std::shared_ptr< std::map<int, grampcd::CouplingModelPtr> >& model) const
    {
        std::unique_lock<std::mutex> guard(comm_data->mutex_couplingModels_);
        comm_data->flag_couplingModels_ = true;
        comm_data->couplingModel_for_simulation_ = model;
        comm_data->cond_var_couplingModels_.notify_one();
    }

    void CommunicationInterfaceLocal::fromCommunication_get_couplingModelFromAgent(const CommunicationDataPtr& comm_data) const
    {
        std::shared_lock<std::shared_mutex> guard(mutex_basics_);
		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_send_coupling_models>(agent_->get_couplingModels(), agent_->get_id()));
		async_send(get_communicationData("coordinator"), serialize(message));
    }

    void CommunicationInterfaceLocal::fromCommunication_send_agentModel(const CommunicationDataPtr& comm_data, const grampcd::AgentModelPtr& model) const
    {
        std::unique_lock<std::mutex> guard(comm_data->mutex_agentModel_);
        comm_data->flag_agentModel_ = true;
        comm_data->agentModel_for_simulation_ = model;
        comm_data->cond_var_agentModel_.notify_one();
    }

    void CommunicationInterfaceLocal::fromCommunication_get_agentModelFromAgent(const CommunicationDataPtr& comm_data) const
    {
        std::shared_lock<std::shared_mutex> guard(mutex_basics_);
		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_send_agent_model>(agent_->get_agentModel(), agent_->get_id()));
		async_send(get_communicationData("coordinator"), serialize(message));
    }

    void CommunicationInterfaceLocal::fromCommunication_get_agentState_for_simulation(const CommunicationDataPtr& comm_data) const
    {
        // if mutex_basics_ is used here, there is a deadlock if neighbor approximation is used
		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_send_agent_state_for_simulation>(agent_->get_agentState(), agent_->get_id()));
		async_send(comm_data, serialize(message));
    }

    void CommunicationInterfaceLocal::fromCommunication_send_agentState_for_simulation(const CommunicationDataPtr& comm_data, const AgentStatePtr& state) const
    {
        std::unique_lock<std::mutex> guard(comm_data->mutex_agentState_for_simulation_);
        comm_data->flag_agentState_for_simulation_ = true;
        comm_data->agentState_for_simulation_ = state;
        comm_data->cond_var_agentState_for_simulation_.notify_all();
    }

    void CommunicationInterfaceLocal::fromCommunication_send_convergenceFlag(const CommunicationDataPtr& comm_data, bool converged, int from)
    {
        std::unique_lock<std::shared_mutex> guard_coordinator(mutex_coordinator_);
        coordinator_->fromCommunication_received_convergenceFlag(converged, from);

        std::unique_lock<std::mutex> guard_trigger(mutex_triggerStep_);
        --numberOfNotifications_triggerStep_;
        conditionVariable_triggerStep_.notify_one();
    }

    void CommunicationInterfaceLocal::fromCommunication_send_stoppedAlgFlag(const CommunicationDataPtr& comm_data, const bool flag, const int from)
    {
        // received flag that agent stopped admm
        std::unique_lock<std::shared_mutex> guard(mutex_basics_);
        agent_->fromCommunication_received_flagStoppedAdmm(flag, from);
    }

    void CommunicationInterfaceLocal::fromCommunication_send_stoppedAlgFlagCoordinator(const CommunicationDataPtr& comm_data, const bool flag, const int from)
    {
        // coordinator received flag that agent stopped admm
        std::unique_lock<std::shared_mutex> guard(mutex_coordinator_);
        coordinator_->fromCommunication_received_flagStoppedAlg(flag, from);
    }

    void CommunicationInterfaceLocal::fromCommunication_send_flagToStopAdmm(const CommunicationDataPtr& comm_data, const bool flag)
    {
        std::unique_lock<std::shared_mutex> guard(mutex_basics_);
        agent_->fromCommunication_received_flagToStopAdmm(flag);
    }

    void CommunicationInterfaceLocal::fromCommunication_send_multiplierPenaltyState(const CommunicationDataPtr& comm_data, 
        const MultiplierState& multiplier, const PenaltyState& penalty, const int from)
    {
        std::unique_lock<std::shared_mutex> guard(mutex_basics_);
        agent_->fromCommunication_received_multiplierState(multiplier, penalty, from);
    }

    void CommunicationInterfaceLocal::fromCommunication_send_desiredAgentState(const CommunicationDataPtr& comm_data, const grampcd::AgentStatePtr& state) const
    {
        std::unique_lock<std::mutex> guard(comm_data->mutex_desired_agentState_);
        comm_data->flag_desired_agentState_ = true;
        comm_data->desired_agentState_ = state;
        comm_data->cond_var_desired_agentState_.notify_all();
    }

    void CommunicationInterfaceLocal::fromCommunication_send_localCopies(const CommunicationDataPtr& comm_data, const AgentState& state, int from)
    {
        std::unique_lock<std::shared_mutex> guard(mutex_basics_);
        agent_->fromCommunication_received_localCopies(state, from);
    }

    void CommunicationInterfaceLocal::fromCommunication_send_agentState(const CommunicationDataPtr& comm_data, const AgentState& state, 
        const ConstraintState& constr_state, int from)
    {
        std::unique_lock<std::shared_mutex> guard(mutex_basics_);
        agent_->fromCommunication_received_agentState(state, constr_state, from);
    }

    void CommunicationInterfaceLocal::fromCommunication_send_numberOfNeighbors(const CommunicationDataPtr& comm_data, int number, int from)
    {
        std::unique_lock<std::shared_mutex> guard(mutex_basics_);
        agent_->fromCommunication_received_numberOfNeighbors(number, from);
    }

    void CommunicationInterfaceLocal::fromCommunication_deregister_coupling(const CommunicationDataPtr& comm_data, const CouplingInfoPtr& info)
    {
        std::unique_lock<std::shared_mutex> guard(mutex_coordinator_);
        coordinator_->deregister_coupling(info);
    }

    void CommunicationInterfaceLocal::fromCommunication_deregister_agent(const CommunicationDataPtr& comm_data, const AgentInfo& info)
    {
        std::unique_lock<std::shared_mutex> guard_coordinator(mutex_coordinator_);
        if (!coordinator_->deregister_agent(info))
        {
            log_->print(DebugType::Error) << "[CommunicationInterfaceLocal::fromCommunication_deregister_agent] "
                << "Deregister agent with id " << info.id_ << " failed." << std::endl;

            return;
        }
		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_successfully_deregistered_agent>(info));
		async_send(comm_data, serialize(message));
    }

    void CommunicationInterfaceLocal::fromCommunication_successfully_deregistered_agent(const CommunicationDataPtr& comm_data, const AgentInfo& info)
    {
        std::unique_lock<std::shared_mutex> guard(mutex_infos_);
        DataConversion::erase_element_from_vector(agent_infos_active_, info);
    }

    void CommunicationInterfaceLocal::fromCommunication_get_solution(const CommunicationDataPtr& comm_data) const
    {
        std::shared_lock<std::shared_mutex> guard(mutex_basics_);
		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_send_solution>(agent_->get_solution()));
		async_send(comm_data, serialize(message));
    }

    void CommunicationInterfaceLocal::fromCommunication_send_solution(const CommunicationDataPtr& comm_data, const SolutionPtr& solution) const
    {
        std::unique_lock<std::mutex> guard(mutex_getSolutions_);
        comm_data->solution_ = solution;
        --numberOfNotifications_getSolutions_;
        conditionVariable_getSolutions_.notify_one();
    }

    void CommunicationInterfaceLocal::fromCommunication_received_acknowledgement_executed_AlgStep(const CommunicationDataPtr& comm_data) const
    {
        std::unique_lock<std::mutex> guard(mutex_triggerStep_);
        --numberOfNotifications_triggerStep_;
        conditionVariable_triggerStep_.notify_one();
    }

    void CommunicationInterfaceLocal::fromCommunication_received_acknowledgement_received_optimizationInfo(const CommunicationDataPtr& comm_data) const
    {
        std::unique_lock<std::mutex> guard(mutex_config_optimizationInfo_);
	    --numberOfNotifications_config_optimizationInfo_;
	    comm_data->is_configured = true;
        conditionVariable_config_optimizationInfo_.notify_one();
    }

    void CommunicationInterfaceLocal::fromCommunication_triggerStep(AlgStep step)
    {
        std::unique_lock<std::shared_mutex> guard(mutex_basics_);
        agent_->fromCommunication_trigger_step(step);

        if (step != AlgStep::GEN_SEND_CONVERGENCE_FLAG)
        {
			const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_acknowledge_executed_ADMM_step>());
			async_send(get_communicationData("coordinator"), serialize(message));
        }
    }

    void CommunicationInterfaceLocal::fromCommunication_send_optimizationInfo(const CommunicationDataPtr& comm_data, const OptimizationInfo& optimization_info)
    {
        std::unique_lock<std::shared_mutex> guard(mutex_basics_);
        agent_->fromCommunication_configured_optimization(optimization_info);

		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_acknowledge_received_optimization_info>());
		async_send(get_communicationData("coordinator"), serialize(message));
    }

    void CommunicationInterfaceLocal::fromCommunication_successfully_registered_coupling(const CommunicationDataPtr& comm_data, 
        const CouplingInfo& coupling_info, const AgentInfo& agent_info, const CommunicationInfoPtr& communication_info)
    {
        std::unique_lock<std::shared_mutex> guard(mutex_infos_);

        // add coupling to list of active couplings
        coupling_infos_active_.push_back(coupling_info);

        // check if coupling is pending and delete it
	    DataConversion::erase_element_from_vector(coupling_infos_pending_, coupling_info);
	    guard.unlock();

        // register coupling in agent
        std::unique_lock<std::shared_mutex> guard_mutex_basics(mutex_basics_);

        agent_->fromCommunication_registered_coupling(coupling_info, agent_info);

        guard_mutex_basics.unlock();

        connect_to_neighbor(communication_info);
    }

    void CommunicationInterfaceLocal::fromCommunication_successfully_registered_agent(const CommunicationDataPtr& comm_data, const AgentInfo& info)
    {
        std::unique_lock<std::shared_mutex> guard(mutex_infos_);

        for (unsigned int i = 0; i < agent_infos_registering_.size(); ++i)
        {
            if (agent_infos_registering_[i].id_ == info.id_)
            {
                agent_infos_active_.push_back(agent_infos_registering_[i]);
                agent_infos_registering_.erase(agent_infos_registering_.begin() + i);
            }
        }

        // send communication info to coordinator
        std::shared_lock<std::shared_mutex> guard_socket(comm_data->mutex_socket_);

        comm_info_local_.ip_ = comm_data->socket_.local_endpoint().address().to_string();
		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_send_communication_info>(comm_info_local_));
		async_send(get_communicationData("coordinator"), serialize(message));
    }

    void CommunicationInterfaceLocal::fromCommunication_register_agent(const CommunicationDataPtr& comm_data, const AgentInfoPtr& info)
    {
        std::unique_lock<std::shared_mutex> guard(mutex_coordinator_);
        if (!coordinator_->register_agent(*info)) 
            return;
        guard.unlock();

        // send acknowledge
		const auto message = std::static_pointer_cast<Message>(std::make_shared<Message_successfully_registered_agent>(*info));
		async_send(comm_data, serialize(message));

        // add agent to internal list
        comm_data->agent_info_ = info;
        comm_data->communication_info_->id_ = info->id_;
    }

    void CommunicationInterfaceLocal::fromCommunication_register_coupling(const CommunicationDataPtr& comm_data, const CouplingInfo& coupling_info)
    {
        const auto comm_data_agent = get_communicationData(coupling_info.agent_id_);

        if (comm_data_agent == nullptr)
            return;

        // if the communication info of the agent has not been arrived, return
        if (comm_data_agent->communication_info_->port_ == "")
            return;

        const auto comm_data_neighbor = get_communicationData(coupling_info.neighbor_id_);
        if (comm_data_neighbor == nullptr)
            return;

        // if the communication info of the agent has not been arrived, return
        if (comm_data_neighbor->communication_info_->port_ == "")
            return;

        std::unique_lock<std::shared_mutex> guard(mutex_coordinator_);
        if (!coordinator_->register_coupling(coupling_info))
            return;
    
	    // send coupling to agent
		const auto message_agent = std::static_pointer_cast<Message>(std::make_shared<Message_successfully_registered_coupling>(coupling_info, *comm_data_neighbor->agent_info_, *comm_data_neighbor->communication_info_));
        async_send(comm_data_agent, serialize(message_agent));

		const auto message_neighbor = std::static_pointer_cast<Message>(std::make_shared<Message_successfully_registered_coupling>(coupling_info, *comm_data_agent->agent_info_, *comm_data_agent->communication_info_));
		async_send(comm_data_neighbor, serialize(message_neighbor));
    
    }

    void CommunicationInterfaceLocal::fromCommunication_deregistered_coupling(const CommunicationDataPtr& comm_data, const CouplingInfo& info)
    {
        std::unique_lock<std::shared_mutex> guard_basics(mutex_basics_);
        agent_->fromCommunication_deregistered_coupling(info);

        // check if coupling should be re-registered
        const bool reregister = (!DataConversion::is_element_in_vector(coupling_infos_blocked_, info))
            && (info.agent_id_ == agent_->get_id());

        std::unique_lock<std::shared_mutex> guard_infos(mutex_infos_);

        // if coupling should not be re-registered, delete it form the list pending
        if (!reregister)
            DataConversion::erase_element_from_vector(coupling_infos_pending_, info);
        else if(reregister && !DataConversion::is_element_in_vector(coupling_infos_pending_, info))
            coupling_infos_pending_.push_back(info);

        // delete it from the list deregistering
        DataConversion::erase_element_from_vector(coupling_infos_deregistering_, info);

        // delete it from the list active
        DataConversion::erase_element_from_vector(coupling_infos_active_, info);

        guard_infos.unlock();

        // Define neighbor id
        int neighbor_id;
        if (info.agent_id_ != agent_->get_id()) 
            neighbor_id = info.agent_id_;
        else 
            neighbor_id = info.neighbor_id_;

        // check if agent is still neighbor
        const bool isNeighbor = DataConversion::is_element_in_vector(agent_->get_neighbors(), neighbor_id);

        // check if there is still a open connection
        const auto comm_data_neighbor = get_communicationData(neighbor_id);

        if (comm_data_neighbor == nullptr)
            return;
    
        // if agent is not a neighbor anymore, close the connection
        if (!isNeighbor)
            handle_disconnect(comm_data_neighbor);
    }

    const unsigned int CommunicationInterfaceLocal::get_numberOfAgents() const
    {
        if (!coordinator_)
        {
            log_->print(DebugType::Error) << "CommunicationInterfaceCentral::get_numberOfAgents"
                << "This method requires the coordinator." << std::endl;
            return 1;
        }

        return coordinator_->get_numberOfAgents();
    }

    const LoggingPtr& CommunicationInterfaceLocal::get_log() const
    {
        return log_;
	}

	const std::shared_ptr< std::vector<char> > CommunicationInterfaceLocal::serialize(const MessagePtr message) const
	{
        try
		{
			std::stringstream sstream;
			cereal::BinaryOutputArchive oarchive(sstream);

			// serialize data
			oarchive(message);

			const auto data_string = sstream.str();

			// generate output
			const auto size_of_packet = static_cast<unsigned int>(data_string.size() + Constants::SIZE_OF_HEADERS_);
			const auto data = std::make_shared<std::vector<char>>(size_of_packet, 0);

			// insert header
			std::to_chars(&data->at(0), &data->at(Constants::SIZE_OF_HEADERS_), size_of_packet);

			// insert data
			for (int i = 0; i < data_string.size(); ++i)
				(*data)[Constants::SIZE_OF_HEADERS_ + i] = data_string[i];

			return data;
        }
        catch (cereal::Exception ec)
		{
			log_->print(DebugType::Base) << "[CommunicationInterfaceLocal::serialize]: "
                << "Serializing message with type " << static_cast<int>(message->get_message_type()) 
                << " failed." << std::endl;
            log_->print(DebugType::Base) << "Error message: " << ec.what() << std::endl;
        }

        return std::shared_ptr< std::vector<char> >();
	}

	void CommunicationInterfaceLocal::evaluate_data(CommunicationDataPtr comm_data, const std::shared_ptr<std::stringstream>& data)
	{
        try
		{
			cereal::BinaryInputArchive iarchive(*data);
			MessagePtr message;
			iarchive(message);

			message_handler_->handle_message(comm_data, message);
        }
        catch (cereal::Exception ec)
        {
			log_->print(DebugType::Base) << "[CommunicationInterfaceLocal::evaluate_data]: "
				<< "Deserialization failed" << std::endl;
			log_->print(DebugType::Base) << "Error message: " << ec.what() << std::endl;
        }
	}

   

}
