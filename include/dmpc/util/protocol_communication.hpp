
#ifndef PROTOCOL_COMMUNICATION_HPP
#define PROTOCOL_COMMUNICATION_HPP

#include "dmpc/optim/optim_util.hpp"
#include "dmpc/optim/solution.hpp"
#include "dmpc/info/communication_data.hpp"
#include "dmpc/comm/communication_interface_local.hpp"

#include "dmpc/info/agent_info.hpp"
#include "dmpc/info/coupling_info.hpp"
#include "dmpc/info/optimization_info.hpp"
#include "dmpc/info/communication_info.hpp"

#include "dmpc/state/agent_state.hpp"
#include "dmpc/state/coupling_state.hpp"
#include "dmpc/state/multiplier_state.hpp"
#include "dmpc/state/penalty_state.hpp"

namespace dmpc
{
	class ProtocolCommunication
	{
	public:
		// main function
		static void evaluateData(CommunicationDataPtr comm_data, CommunicationInterfaceLocal* communication_interface, const std::vector<char>& data);

		/*
		build protocol
		*/

		// common
		static const std::shared_ptr< std::vector<char> > buildProtocol_get_solution();
		static const std::shared_ptr< std::vector<char> > buildProtocol_send_solution(const SolutionPtr& solution);
		static const std::shared_ptr< std::vector<char> > buildProtocol_send_convergenceFlag(const bool converged, const int from);
		static const std::shared_ptr< std::vector<char> > buildProtocol_send_numberOfNeighbors(const int number, const int from);
		static const std::shared_ptr< std::vector<char> > buildProtocol_triggerStep(const ADMMStep& step);
		static const std::shared_ptr< std::vector<char> > buildProtocol_get_numberOfActiveCouplings();
		static const std::shared_ptr< std::vector<char> > buildProtocol_send_numberOfActiveCouplings(int number);
		static const std::shared_ptr< std::vector<char> > buildProtocol_get_ping();
		static const std::shared_ptr< std::vector<char> > buildProtocol_send_ping();
		static const std::shared_ptr< std::vector<char> > buildProtocol_send_flagToAgents();

		// info
		static const std::shared_ptr< std::vector<char> > buildProtocol_send_optimizationInfo(const OptimizationInfo& info);
		static const std::shared_ptr< std::vector<char> > buildProtocol_send_communicationInfo(const CommunicationInfo& info);

		// model
		static const std::shared_ptr< std::vector<char> > buildProtocol_get_agentModel_from_agent();
		static const std::shared_ptr< std::vector<char> > buildProtocol_get_couplingModel_from_agent();
		static const std::shared_ptr< std::vector<char> > buildProtocol_send_agentModel(const AgentModelPtr& model, const int from);
		static const std::shared_ptr< std::vector<char> > buildProtocol_send_couplingModel(const std::map<int, CouplingModelPtr>& map, const int from);

		// register and deregister
		static const std::shared_ptr< std::vector<char> > buildProtocol_register_agent(const AgentInfo& info);
		static const std::shared_ptr< std::vector<char> > buildProtocol_successfully_registered_agent(const AgentInfo& info);
		static const std::shared_ptr< std::vector<char> > buildProtocol_deregister_agent(const AgentInfo& info);
		static const std::shared_ptr< std::vector<char> > buildProtocol_successfully_deregistered_agent(const AgentInfo& info);
		static const std::shared_ptr< std::vector<char> > buildProtocol_deregister_coupling(const CouplingInfo& info);
		static const std::shared_ptr< std::vector<char> > buildProtocol_fromCommunication_deregistered_coupling(const CouplingInfo& info);
		static const std::shared_ptr< std::vector<char> > buildProtocol_successfully_registered_coupling(const CouplingInfo& coupling_info, const AgentInfo& agent_info, const CommunicationInfo& comm_info);
		static const std::shared_ptr< std::vector<char> > buildProtocol_register_coupling(const CouplingInfo& coupling_info);

		// states
		static const std::shared_ptr< std::vector<char> > buildProtocol_send_agentState(const AgentState& state, const int from);
		static const std::shared_ptr< std::vector<char> > buildProtocol_send_couplingState(const CouplingState& state, const int from);
		static const std::shared_ptr< std::vector<char> > buildProtocol_send_couplingState(const CouplingState& state1, const CouplingState& state2, const int from);
		static const std::shared_ptr< std::vector<char> > buildProtocol_send_multiplierPenaltyState(const MultiplierState& state, const PenaltyState& penalty, const int from);
		static const std::shared_ptr< std::vector<char> > buildProtocol_get_agentState_for_simulation();
		static const std::shared_ptr< std::vector<char> > buildProtocol_get_desiredAgentState_from_agent();
		static const std::shared_ptr< std::vector<char> > buildProtocol_send_agentState_for_simulation(const AgentState& state, const int from);
		static const std::shared_ptr< std::vector<char> > buildProtocol_send_desiredAgentState(const AgentState& state, const int from);
		static const std::shared_ptr< std::vector<char> > buildProtocol_send_simulatedState(const std::vector<typeRNum>& x_next, typeRNum dt, typeRNum t0);

		// acknowledgments
		static const std::shared_ptr< std::vector<char> > buildProtocol_acknowledge_received_optimizationInfo();
		static const std::shared_ptr< std::vector<char> > buildProtocol_acknowledge_executed_ADMMstep();

		/*
		build from protocol
		*/

		// common
		static const ADMMStep buildFromProtocol_ADMMStep(const std::vector<char>& data);
		static const SolutionPtr buildFromProtocol_solution(const std::vector<char>& data);
		static const int buildFromProtocol_numberOfNeighbors(const std::vector<char>& data);
		static const bool buildFromProtocol_convergenceFlag(const std::vector<char>& data);
		static const int buildFromProtocol_from(const std::vector<char>& data);
		static const int buildFromProtocol_numberOfActiveCouplings(const std::vector<char>& data);

		//infos
		static const std::shared_ptr< AgentInfo > buildFromProtocol_agentInfo(const std::vector<char>& data);
		static const std::shared_ptr< OptimizationInfo > buildFromProtocol_optimizationInfo(const std::vector<char>& data);
		static const std::shared_ptr< CouplingInfo > buildFromProtocol_couplingInfo(const std::vector<char>& data);
		static const std::shared_ptr< std::vector<CommunicationInfo> > buildFromProtocol_communicationInfoArray(const std::vector<char>& data);
		static const std::shared_ptr< CommunicationInfo > buildFromProtocol_communicationInfo(const std::vector<char>& data);
		static const std::shared_ptr< AgentInfo > buildFromProtocol_agentInfo_from_successfully_registered_coupling(const std::vector<char>& data);
		static const std::shared_ptr< CouplingInfo > buildFromProtocol_couplingInfo_from_successfully_registered_coupling(const std::vector<char>& data);
		static const std::shared_ptr< CommunicationInfo > buildFromProtocol_communicationInfo_from_successfully_registered_coupling(const std::vector<char>& data);

		// states
		static const AgentStatePtr buildFromProtocol_agentState(const std::vector<char>& data);
		static const CouplingStatePtr buildFromProtocol_couplingState(const std::vector<char>& data);
		static const CouplingStatePtr buildFromProtocol_couplingState_one(const std::vector<char>& data);
		static const CouplingStatePtr buildFromProtocol_couplingState_two(const std::vector<char>& data);
		static const MultiplierStatePtr buildFromProtocol_multiplierState(const std::vector<char>& data);
		static const PenaltyStatePtr buildFromProtocol_penaltyState(const std::vector<char>& data);
		static const typeRNum buildFromProtocol_t0_from_simulatedState(const std::vector<char>& data);
		static const typeRNum buildFromProtocol_dt_from_simulatedState(const std::vector<char>& data);
		static const std::shared_ptr< std::vector<typeRNum> > buildFromProtocol_xnext_from_simulatedState(const std::vector<char>& data);

		// models
		static const AgentModelPtr buildFromProtocol_agentModel(const std::vector<char>& data, const LoggingPtr& log);
		static const std::shared_ptr< std::map<int, CouplingModelPtr> > buildFromProtocol_couplingModel(const std::vector<char>& data, const LoggingPtr& log);
	};

	const char position_of_index_in_protocol_ = 4;
	const char first_element_with_data_ = 5;

	enum class index : char
	{
		register_agent = 1,
		register_coupling = 2,
		successfully_registered_coupling = 3,
		send_optimizationInfo = 6,
		send_numberOfNeighbors = 7,
		deregister_agent = 8,
		deregistered_coupling = 9,
		send_agentState = 10,
		send_desiredAgentState = 11,
		send_couplingState = 12,
		send_two_couplingStates = 13,
		send_multiplierPenaltyState = 14,
		get_agentState_for_simulation = 15,
		send_convergenceFlag = 16,
		deregister_coupling = 17,
		send_agentState_for_simulation = 19,
		triggerStep = 20,
		get_agentModelFromAgent = 22,
		send_agentModel = 23,
		get_couplingModelFromAgent = 24,
		send_couplingModel = 25,
		send_simulatedState = 26,
		get_desiredAgentStateFromAgent = 27,
		get_ping = 28,
		send_ping = 29,
		get_numberOfActiveCouplings = 30,
		send_numberOfActiveCouplings = 31,
		send_flagToAgents = 32,
		send_communicationInfo = 34,
		successfully_deregistered_agent = 35,
		successfully_registered_agent = 100,
		received_acknowledgement_received_optimizationInfo = 101,
		received_acknowledgement_executed_ADMMstep = 102,
		get_solution = 103,
		send_solution = 104
	};
}

#endif // PROTOCOL_COMMUNICATION_HPP
