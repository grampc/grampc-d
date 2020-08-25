/* This file is part of GRAMPC-D - (https://github.com/DanielBurk/GRAMPC-D.git)
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

#ifndef COMMUNICATION_INTERFACE_LOCAL_HPP
#define COMMUNICATION_INTERFACE_LOCAL_HPP

#include "dmpc/comm/communication_interface.hpp"
#include "dmpc/coord/coordinator.hpp"
#include "dmpc/simulator/simulator.hpp"
#include "dmpc/info/communication_data.hpp"
#include "asio.hpp"
#include <shared_mutex>
#include <functional>

namespace dmpc
{

/*!
 * @brief Interface for communication between agents
 *
 * Implementation of a communication interface that has direct access to all agents.
 * This means that all agents run within the same process.
*/
class CommunicationInterfaceLocal : public CommunicationInterface
{
public:

    ~CommunicationInterfaceLocal();
    CommunicationInterfaceLocal(const LoggingPtr& log, const unsigned short port);
    CommunicationInterfaceLocal(const LoggingPtr& log, const CommunicationInfo& comm_info_coordinator);

	/*Returns the communication data of an agent.*/
	const CommunicationDataPtr get_communicationData(const int agent_id) const;
	/*Returns the communication data of the coordinator.*/
	const CommunicationDataPtr get_communicationData(const std::string& agent_type) const;

	/*Set the coordinator.*/
    void set_coordinator(const CoordinatorPtr& coordinator);
	/*Set an agent.*/
    void set_agent(const AgentPtr& agent);
	/*Set the simulator.*/
    void set_simulator(const SimulatorPtr& simulator);
	/*Set communication info for coordinator.*/
    void set_communicationInfo_for_coordinator(const CommunicationInfo& info);
	/*Returns a pointer to log*/
	const LoggingPtr& get_log() const;

	/*Register an agent.*/
	const bool register_agent(const AgentPtr& agent) override;
	/*De-register an agent.*/
	const bool deregister_agent(const AgentInfo& agent) override;
	/*Register coupling between agents.*/
	const bool register_coupling(const CouplingInfo& coupling) override;
	/*De-register coupling between agents.*/
	const bool deregister_coupling(const CouplingInfo& coupling) override;

	/*This function is called if a message arrived to de-register a coupling..*/
	const bool fromCommunication_deregistered_coupling(const CouplingInfo& coupling) override;

	/*Send number of neighbors to an agent.*/
	const bool send_numberOfNeighbors(const int number, const int from, const int to) override;
	const /*Send an agent state to an agent.*/
	bool send_agentState(const AgentState& state, const int from, const int to) override;
	/*Send desired agent state to an agent.*/
	const bool send_desiredAgentState(const AgentState& desired_state, const int from, const int to) override;
	/*Send coupling state to an agent.*/
	const bool send_couplingState(const CouplingState& state, const int from, const int to) override;
	/*Send two coupling states to an agent.*/
	const bool send_couplingState(const CouplingState& state, const CouplingState& state2, const int from, const int to) override;
	/*Send multiplier and penalty states to an agent.*/
	const bool send_multiplierState(const MultiplierState& state, PenaltyState penalty, const int from, const int to) override;
	/*Send convergence flag to the coordinator.*/
	const bool send_convergenceFlag(const bool converged, const int from) override;

	/*Configure optimization.*/
	const bool configure_optimization(const OptimizationInfo& info) override;
	void configureOptimization(const CommunicationDataPtr& comm_data);

	/*Trigger a step of the ADMM algorithm.*/
	const bool trigger_step(const ADMMStep& step) override;
	/*Trigger simulation.*/
	void trigger_simulation(const std::string& Integrator, const typeRNum dt) override;

	/*Return the agent state of an agent.*/
	const AgentStatePtr get_agentState_from_agent(const int agentId) const override;
	/*Return the desired agent state of an agent.*/
	const AgentStatePtr get_desiredAgentState_from_agent(const int agentId) const override;
	/*Returns a map with agent infos.*/
	const std::shared_ptr< std::map<unsigned int, AgentInfoPtr > > get_agentInfo_from_coordinator() const override;
	/*Return a map with sending neighbors.*/
	const std::shared_ptr< std::map< unsigned int, std::vector< CouplingInfoPtr > > > get_sendingNeighbors_from_coordinator() const override;
	/*Return a map with receiving neighbors.*/
	const std::shared_ptr< std::map< unsigned int, std::vector< CouplingInfoPtr > > > get_receivingNeighbors_from_coordinator() const override;
	/*Return the agent model of an agent.*/
	const AgentModelPtr get_agentModel(const int agentId) const override;
	/*Return all coupling models of an agent.*/
	const std::shared_ptr< std::map<int, CouplingModelPtr> > get_couplingModels_from_agent(const int agentId) const override;

	/*Set simulated state of an agent.*/
	void set_simulatedState_for_agent(const int agentId, const std::vector<typeRNum>& new_state, const typeRNum dt, const typeRNum t0) override;
	/*Returns the current solution of an agent.*/
	const SolutionPtr get_solution(const unsigned int agent_id) const override;
	/*Return the current solution of a set of agents.*/
	const std::vector< SolutionPtr > get_solution(const std::string& agents) const override;
	/*Resets the solution of an agent.*/
	void reset_solution(const unsigned int agent_id) override;
	/*Resets the solution of a set of agents.*/
	void reset_solution(const std::string& agents) override;

	/*Wait for connections.*/
	void waitFor_connection(const int agents, const int couplings) override;
	/*Wait for a flag of the coordinator.*/
	void waitFor_flag_from_coordinator() override;
	/*Send flag to a set of agents.*/
	void send_flag_to_agents(const std::string& agents) const override;
	/*Send flag to a set of agents.*/
	void send_flag_to_agents(const std::vector<int>& agent_ids) const override;
	/*Send flag to an agent.*/
	void send_flag_to_agents(const int agent_id) const override;
	/*Set to passive mode.*/
	void set_passive() override;

	/*Cap the stored data to a number of data points.*/
	void cap_stored_data(const unsigned int data_points) override;
	/*Returns the number of agents.*/
	const unsigned int get_numberOfAgents() const override;

	/*This function is called if number of neighbors is received.*/
	void fromCommunication_send_numberOfNeighbors(const CommunicationDataPtr& comm_data, const int number, const int from);
	/*This function is called convergence flag is received.*/
	void fromCommunication_send_convergenceFlag(const CommunicationDataPtr& comm_data, const bool converged, const int from);
	/*This function is called if requirement for ping is received.*/
	void fromCommunication_get_ping(const CommunicationDataPtr& comm_data) const;
	/*This function is called if ping is received.*/
	void fromCommunication_send_ping(const CommunicationDataPtr& comm_data) const;
	/*This function is called if requirement for number of active couplings is received.*/
	void fromCommunication_get_numberOfActiveCouplings(const CommunicationDataPtr& comm_data) const;
	/*This function is called if number of active couplings is received.*/
	void fromCommunication_send_numberOfActiveCouplings(const CommunicationDataPtr& comm_data, const int number) const;
	/*This function is called if flag is received.*/
    void fromCommunication_send_flagToAgents(const CommunicationDataPtr& comm_data) const;

	/*This function is called if requirement for agent model is received.*/
    void fromCommunication_get_agentModelFromAgent(const CommunicationDataPtr& comm_data) const;
	/*This function is called if agent model is received.*/
    void fromCommunication_send_agentModel(const CommunicationDataPtr& comm_data, const dmpc::AgentModelPtr& model) const;
	/*This function is called if requirement for coupling model is received.*/
    void fromCommunication_get_couplingModelFromAgent(const CommunicationDataPtr& comm_data) const;
	/*This function is called if coupling model is received.*/
    void fromCommunication_send_couplingModel(const CommunicationDataPtr& comm_data, const std::shared_ptr< std::map<int, dmpc::CouplingModelPtr> >& model) const;

	/*This function is called if agent state is received.*/
	void fromCommunication_send_agentState(const CommunicationDataPtr& comm_data, const dmpc::AgentStatePtr& state, const int from);
	/*This function is called if desired agent state is received.*/
	void fromCommunication_send_desiredAgentState(const CommunicationDataPtr& comm_data, const dmpc::AgentStatePtr& state) const;
	/*This function is called if multiplier state is received.*/
    void fromCommunication_send_multiplierPenaltyState( const CommunicationDataPtr& comm_data, const dmpc::MultiplierStatePtr& multiplier,
        const dmpc::PenaltyStatePtr& penalty, const int from);
	/*This function is called if requirement for agent state for simulation is received.*/
    void fromCommunication_get_agentState_for_simulation(const CommunicationDataPtr& comm_data) const;
	/*This function is called if agent state for simulation is received.*/
    void fromCommunication_send_agentState_for_simulation(const CommunicationDataPtr& comm_data, const dmpc::AgentStatePtr& state) const;
	/*This function is called if simulated state is received.*/
    void fromCommunication_send_simulatedState(const CommunicationDataPtr& comm_data, const std::shared_ptr< std::vector<typeRNum> >& x_next, const typeRNum dt, const typeRNum t0);
	/*This function is called if requirement for agent state is received.*/
	void fromCommunication_get_desiredAgentStateFromAgent(const CommunicationDataPtr& comm_data) const;
	/*This function is called if coupling state is received.*/
    void fromCommunication_send_couplingState(const CommunicationDataPtr& comm_data, const dmpc::CouplingStatePtr& state, const int from);
	/*This function is called if coupling state is received.*/
    void fromCommunication_send_couplingState(const CommunicationDataPtr& comm_data, const dmpc::CouplingStatePtr& state1, const dmpc::CouplingStatePtr& state2, const int from);

	/*This function is called agent should be registered.*/
	void fromCommunication_register_agent(const CommunicationDataPtr& comm_data, const AgentInfoPtr& info);
	/*This function is called if agent is successfully registered.*/
	void fromCommunication_successfully_registered_agent(const CommunicationDataPtr& comm_data, const AgentInfoPtr& info);
	/*This function is called if agent should be de-registered.*/
	void fromCommunication_deregister_agent(const CommunicationDataPtr& comm_data, const AgentInfoPtr& info);
	/*This function is called if agent is successfully de-registered.*/
	void fromCommunication_successfully_deregistered_agent(const CommunicationDataPtr& comm_data, const AgentInfoPtr& info);

	/*This function is called if optimization info is received.*/
    void fromCommunication_send_optimizationInfo(const CommunicationDataPtr& comm_data, const OptimizationInfoPtr& optimization_info);
	/*This function is called if communication info is received.*/
    void fromCommunication_send_communicationInfo(const CommunicationDataPtr& comm_data, const CommunicationInfoPtr& info);

	/*This function is called if acknowledgment for received optimization info is received.*/
	void fromCommunication_received_acknowledgement_received_optimizationInfo(const CommunicationDataPtr& comm_data) const;
	/*This function is called if acknowledgment for executed ADMM step is received.*/
    void fromCommunication_received_acknowledgement_executed_ADMMstep(const CommunicationDataPtr& comm_data) const;
	/*This function is called if ADMM step should be triggered.*/
    void fromCommunication_triggerStep(ADMMStep step);

	/*This function is called if coupling should be registered.*/
	void fromCommunication_register_coupling(const CommunicationDataPtr& comm_data, const CouplingInfoPtr& info);
	/*This function is called if coupling is successfully registered.*/
    void fromCommunication_successfully_registered_coupling(const CommunicationDataPtr& comm_data,
        const CouplingInfoPtr& coupling_info, const AgentInfoPtr& agent_info, const CommunicationInfoPtr& communication_info);
	/*This function is called if coupling should be de-registered.*/
	void fromCommunication_deregister_coupling(const CommunicationDataPtr& comm_data, const CouplingInfoPtr& info);
	/*This function is called if coupling is successfully de-registered.*/
	void fromCommunication_deregistered_coupling(const CommunicationDataPtr& comm_data, const CouplingInfoPtr& info);

	/*This function is called if requirement for solution is received.*/
    void fromCommunication_get_solution(const CommunicationDataPtr& comm_data) const;
	/*This function is called if solution is received.*/
	void fromCommunication_send_solution(const CommunicationDataPtr& comm_data, const SolutionPtr& solution) const;

private:
	/*Start an server.*/
	void start_server();
	/*Start a client.*/
	void start_client();

	/*Function handler for accepting connections.*/
    void acceptHandler(const std::error_code &ec, CommunicationDataPtr comm_data);
	/*Function handler for receiving messages.*/
    void readHandler(const std::error_code&ec, const size_t &amountOfBytes, CommunicationDataPtr comm_data);
	/*Function handler for sending messages.*/
	void writeHandler(const std::error_code& ec, std::size_t bytes_transferred, CommunicationDataPtr comm_data, std::shared_ptr<std::vector<char>> data_ptr) const;
	/*Function handler for connections.*/
    void connectHandler(const std::error_code &ec, CommunicationDataPtr comm_data);

	/*Start the polling.*/
    void start_polling(const CommunicationDataPtr& comm_data) const;
	/*Poll the registrations.*/
    void poll_registration(const std::error_code &ec, CommunicationDataPtr comm_data) const;
	/*Start the ping.*/
    void start_ping(const CommunicationDataPtr& comm_data);
	/*Poll the ping.*/
    void poll_ping(const std::error_code &ec, CommunicationDataPtr comm_data);
	/*Poll checking the connection.*/
    void poll_check_connection(const std::error_code &ec, CommunicationDataPtr comm_data);

	/*Handle a disconnect.*/
    void handle_disconnect(const CommunicationDataPtr& comm_data);
	/*Handle the disconnect of an agent.*/
    void handle_disconnect_as_agent(const CommunicationDataPtr& comm_data);
	/*Handle the disconnect of the coordinator.*/
    void handle_disconnect_as_coordinator(const CommunicationDataPtr& comm_data);

	/*Connect to a neighbor.*/
    void connect_to_neighbor(const CommunicationInfoPtr& comm_info);
	/*Close a socket.*/
    void close_socket(const CommunicationDataPtr& comm_data) const;
	/*Close and shutdown a socket.*/
    void close_shutdown_socket(const CommunicationDataPtr& comm_data) const;

	/*Send data asynchronously.*/
    void async_send(const CommunicationDataPtr& comm_data, const std::shared_ptr< std::vector<char> >& data) const;
	/*Connect asynchronously.*/
    void async_connect(const CommunicationDataPtr& comm_data);
	/*Connect asynchronously.*/
	void async_connect(const CommunicationDataPtr& comm_data, const std::string& ip, const std::string& port);
	/*Asynchronously accept a connection.*/
    void async_accept(const CommunicationDataPtr& comm_data);
	/*Read data asynchronously.*/
	void async_read_some(const CommunicationDataPtr& comm_data);

    mutable std::mutex mutex_stream_;
    LoggingPtr log_;
    std::ostringstream stream_;

    mutable std::shared_mutex mutex_basics_;
    AgentPtr agent_;
    SimulatorPtr simulator_;
    std::vector<AgentInfo> agents_;
    CommunicationInfo comm_info_local_;
    CommunicationInfo comm_info_coordinator_;
    const int number_of_threads_ = 4;

	mutable std::shared_mutex mutex_coordinator_;
	CoordinatorPtr coordinator_;

    // mutex to protect data
    mutable std::shared_mutex mutex_infos_;

    std::vector<AgentInfo> agent_infos_active_;
	std::vector<AgentInfo> agent_infos_registering_;
	std::vector<AgentInfo> agent_infos_deregistering_;
	std::vector<AgentInfo> agent_infos_blocked_;

    std::vector<CouplingInfo> coupling_infos_active_;
    std::vector<CouplingInfo> coupling_infos_pending_;
    std::vector<CouplingInfo> coupling_infos_deregistering_;
    std::vector<CouplingInfo> coupling_infos_blocked_;

	mutable std::shared_mutex mutex_comm_data_vec_;
	std::vector<CommunicationDataPtr> comm_data_vec_;
	std::vector<CommunicationDataPtr> comm_data_vec_to_delete_;

    asio::io_service ioService_;
    asio::ip::tcp::acceptor acceptor_;
    std::vector< std::thread > threads_for_communication_;
    asio::basic_waitable_timer<std::chrono::system_clock>  timer_waitForAck_;
    asio::basic_waitable_timer<std::chrono::system_clock>  timer_waitTrue_;
    const unsigned int general_waiting_time_s_ = 2;
    const unsigned int polling_period_s_ = 3;
    const unsigned int ping_period_s_ = 2;
    const unsigned int ping_waiting_time_ms_ = 1500;
    const bool disconnect_at_timetout_ = false;

    // variables to wait for: trigger step
    mutable std::mutex mutex_triggerStep_;
    mutable int numberOfNotifications_triggerStep_;
    mutable std::condition_variable conditionVariable_triggerStep_;

	// variables to wait for: get solutions
	mutable std::mutex mutex_getSolutions_;
    mutable int numberOfNotifications_getSolutions_;
    mutable std::condition_variable conditionVariable_getSolutions_;

	// variables to wait for: configure optimization info
    mutable std::mutex mutex_config_optimizationInfo_;
    mutable int numberOfNotifications_config_optimizationInfo_;
    mutable std::condition_variable conditionVariable_config_optimizationInfo_;

	// variables to wait for: wait for connection
    mutable std::mutex mutex_waitForConnection_;
    mutable int numberOfNotifications_waitForConnection_;
    mutable std::condition_variable conditionVariable_waitForConnection_;
};

}

#endif // COMMUNICATION_INTERFACE_LOCAL_HPP
