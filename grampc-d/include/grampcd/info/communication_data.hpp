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

#include "grampcd/optim/solution.hpp"

#include "grampcd/info/agent_info.hpp"
#include "grampcd/info/communication_info.hpp"

#include "grampcd/util/class_forwarding.hpp"

#include "asio.hpp"

#include <shared_mutex>

namespace grampcd
{

    /*@brief This class contains the data required for communication.*/
    class CommunicationData
    {
    public:
        CommunicationData(asio::io_service &ioService)
            : 
            strand_(ioService),
            socket_(ioService),
            timer_check_connection_(ioService),
            timer_ping_(ioService),
            timer_polling_(ioService),
            agent_info_(new AgentInfo),
            communication_info_(new CommunicationInfo),
            solution_(new Solution),
            desired_agentState_(new grampcd::AgentState),
            agentState_for_simulation_(new grampcd::AgentState),
            couplingModel_for_simulation_(new std::map<int, CouplingModelPtr>)
        {}

        /*Id of the corresponding agent.*/
        int id_ = -1;
        /*Agent info of the corresponding agent.*/
        AgentInfoPtr agent_info_;
        /*Flag that states wether flag from coordinator is received.*/
        bool flag_from_coordinator_ = false;

        /*Strand ensured serialization of function calls.*/
        asio::io_service::strand strand_;

        /*Mutex that protects the coscket*/
        std::shared_mutex mutex_socket_;
        /*TCP socket*/
	    asio::ip::tcp::socket socket_;

        /*Mutex that protects the buffer for sending data.*/
        std::mutex mutex_buffer_send_;
        /*Buffer for sending data.*/
	    std::vector< std::shared_ptr<std::vector<char>> > buffer_send_;

        /*Buffer to read data.*/
        std::vector<char> buffer_read_ = std::vector<char>(10000, 0);
        /*Read data.*/
        std::vector<char> data_;
        /*Timer for checking if connection is still alive.*/
        asio::basic_waitable_timer<std::chrono::system_clock> timer_check_connection_;
        /*Flag that states if connection is alive.*/
	    bool is_connected_ = false;
        /*Flag that states if corresponding agent is configured.*/
	    bool is_configured = false;
        /*Flag that states if corresponding agent is configuring.*/
        bool is_configuring = false;

	    /*Communication info of corresponding agent.*/
	    CommunicationInfoPtr communication_info_;
        /*Solution of corresponding agent.*/
	    SolutionPtr solution_;
        /*Desired agent state of corresponding agent.*/
	    grampcd::AgentStatePtr desired_agentState_;
        /*Agent state for simulation of corresponding agent.*/
	    grampcd::AgentStatePtr agentState_for_simulation_;
        /*Agent model of corresponding agent.*/
	    AgentModelPtr agentModel_for_simulation_;
        /*Coupling model for simulation of corresponding agent.*/
	    std::shared_ptr< std::map<int, CouplingModelPtr> > couplingModel_for_simulation_;
        /*Number of connected agents.*/
        int number_of_connected_agents_ = 0;
        /*Timer used for ping.*/
        asio::basic_waitable_timer<std::chrono::system_clock> timer_ping_;
        /*Timer used for polling.*/
        asio::basic_waitable_timer<std::chrono::system_clock> timer_polling_;

        /*Mutex that protects agent state for simulation.*/
        std::mutex mutex_agentState_for_simulation_;
        /*Condition variable used to wait for agent state for simulation.*/
        std::condition_variable cond_var_agentState_for_simulation_;
        /*Flag that states if agent state for simulation is received.*/
        bool flag_agentState_for_simulation_ = false;

        /*Mutex that protects desired agent state.*/
        std::mutex mutex_desired_agentState_;
        /*Condition variable used to wait for desired agent state.*/
        std::condition_variable cond_var_desired_agentState_;
        /*Flag that states whether desired agent state is received.*/
        bool flag_desired_agentState_ = false;

        /*Mutex that protects the agent model.*/
        std::mutex mutex_agentModel_;
        /*Condition variable used to wait for an agent model.*/
        std::condition_variable cond_var_agentModel_;
        /*Flag that states whether agent model is received.*/
        bool flag_agentModel_ = false;

        /*Mutex that protects the coupling models.*/
        std::mutex mutex_couplingModels_;
        /*Condition variable used to wait for coupling models.*/
        std::condition_variable cond_var_couplingModels_;
        /*Flag that states whether coupling models are received.*/
        bool flag_couplingModels_ = false;
    };

}