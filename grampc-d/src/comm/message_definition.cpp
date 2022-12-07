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

#include "grampcd/comm/message.hpp"

#include "grampcd/info/optimization_info.hpp"
#include "grampcd/info/agent_info.hpp"
#include "grampcd/info/coupling_info.hpp"
#include "grampcd/info/communication_info.hpp"

#include "grampcd/optim/solution.hpp"
#include "grampcd/optim/optim_util.hpp"

#include "grampcd/model/agent_model.hpp"
#include "grampcd/model/coupling_model.hpp"
#include "grampcd/state/multiplier_state.hpp"
#include "grampcd/state/penalty_state.hpp"

#include "grampcd/state/agent_state.hpp"
#include "grampcd/state/coupling_state.hpp"
#include "grampcd/state/constraint_state.hpp"



namespace grampcd
{
	Message_send_communication_info::Message_send_communication_info::Message_send_communication_info() {}
	Message_send_communication_info::Message_send_communication_info(const CommunicationInfo communication_info)
		: communication_info_(std::make_shared<CommunicationInfo>(communication_info)) {}
	const Messagetype Message_send_communication_info::get_message_type() const { return Messagetype::SEND_COMMUNICATION_INFO; }

	Message_fromCommunication_deregistered_coupling::Message_fromCommunication_deregistered_coupling() {}
	Message_fromCommunication_deregistered_coupling::Message_fromCommunication_deregistered_coupling(const CouplingInfo coupling_info)
		: coupling_info_(std::make_shared<CouplingInfo>(coupling_info)) {}
	const Messagetype Message_fromCommunication_deregistered_coupling::get_message_type() const { return Messagetype::FROM_COMMUNICATION_DEREGISTERED_COUPLING; }

	const Messagetype Message_get_ping::get_message_type() const { return Messagetype::GET_PING; }

	Message_register_agent::Message_register_agent() {}
	Message_register_agent::Message_register_agent(const AgentInfo& agent_info)
		: agent_info_(std::make_shared<AgentInfo>(agent_info)) {}
	const Messagetype Message_register_agent::get_message_type() const { return Messagetype::REGISTER_AGENT; }

	Message_deregister_agent::Message_deregister_agent() {}
	Message_deregister_agent::Message_deregister_agent(const AgentInfo& agent_info)
		: agent_info_(std::make_shared<AgentInfo>(agent_info)) {}
	const Messagetype Message_deregister_agent::get_message_type() const { return Messagetype::DEREGISTER_AGENT; }

	Message_register_coupling::Message_register_coupling() {}
	Message_register_coupling::Message_register_coupling(const CouplingInfo& coupling_info)
		: coupling_info_(std::make_shared<CouplingInfo>(coupling_info)) {}
	const Messagetype Message_register_coupling::get_message_type() const { return Messagetype::REGISTER_COUPLING; }

	Message_deregister_coupling::Message_deregister_coupling() {}
	Message_deregister_coupling::Message_deregister_coupling(const CouplingInfo& coupling_info)
		: coupling_info_(std::make_shared<CouplingInfo>(coupling_info)) {}
	const Messagetype Message_deregister_coupling::get_message_type() const { return Messagetype::DEREGISTER_COUPLING; }

	Message_send_numberOfNeighbors::Message_send_numberOfNeighbors() {}
	Message_send_numberOfNeighbors::Message_send_numberOfNeighbors(const int number, const int from)
		: number_(number), from_(from) {}
	const Messagetype Message_send_numberOfNeighbors::get_message_type() const { return Messagetype::SEND_NUMBER_OF_NEIGHBORS; }

	Message_send_local_copies::Message_send_local_copies() {}
	Message_send_local_copies::Message_send_local_copies(const AgentState& agent_state, const int from)
		: agent_state_(std::make_shared<AgentState>(agent_state)), from_(from) {}
	const Messagetype Message_send_local_copies::get_message_type() const { return Messagetype::SEND_LOCAL_COPIES; }

	Message_send_desired_agent_state::Message_send_desired_agent_state() {}
	Message_send_desired_agent_state::Message_send_desired_agent_state(const AgentState& desired_agent_state, const int from)
		: desired_agent_state_(std::make_shared<AgentState>(desired_agent_state)),
		from_(from) {}
	const Messagetype Message_send_desired_agent_state::get_message_type() const { return Messagetype::SEND_DESIRED_AGENT_STATE; }

	Message_send_coupling_state::Message_send_coupling_state() {}
	Message_send_coupling_state::Message_send_coupling_state(const CouplingState& coupling_state, const int from)
		: coupling_state_(std::make_shared<CouplingState>(coupling_state)), from_(from) {}
	const Messagetype Message_send_coupling_state::get_message_type() const { return Messagetype::SEND_COUPLING_STATE; }

	Message_send_two_coupling_states::Message_send_two_coupling_states() {}
	Message_send_two_coupling_states::Message_send_two_coupling_states(const CouplingState& state1, const CouplingState& state2, const int from)
		: state1_(std::make_shared<CouplingState>(state1)),
		state2_(std::make_shared<CouplingState>(state2)),
		from_(from) {}
	const Messagetype Message_send_two_coupling_states::get_message_type() const { return Messagetype::SEND_TWO_COUPLING_STATES; }

	Message_send_multiplier_state::Message_send_multiplier_state() {}
	Message_send_multiplier_state::Message_send_multiplier_state(const MultiplierState& multiplier_state, const PenaltyState& penalty_state, const int from)
		: multiplier_state_(std::make_shared<MultiplierState>(multiplier_state)),
		penalty_state_(std::make_shared<PenaltyState>(penalty_state)),
		from_(from) {}
	const Messagetype Message_send_multiplier_state::get_message_type() const { return Messagetype::SEND_MULTIPLIER_STATE; }

	Message_send_agent_state::Message_send_agent_state() {}
	Message_send_agent_state::Message_send_agent_state(const AgentState& agent_state, const ConstraintState& constr_state, const int from)
		: agent_state_(std::make_shared<AgentState>(agent_state)),
		constr_state_(std::make_shared<ConstraintState>(constr_state)),
		from_(from) {}
	const Messagetype Message_send_agent_state::get_message_type() const { return Messagetype::SEND_AGENT_STATE; }

	Message_get_agentState_for_simulation::Message_get_agentState_for_simulation() {}
	const Messagetype Message_get_agentState_for_simulation::get_message_type() const { return Messagetype::GET_AGENT_STATE_FOR_SIMULATION; }

	const Messagetype Message_get_desired_agent_state_from_agent::get_message_type() const { return Messagetype::GET_DESIRED_AGENT_STATE_FROM_AGENT; }

	const Messagetype Message_get_agent_model_from_agent::get_message_type() const { return Messagetype::GET_AGENT_MODEL_FROM_AGENT; }

	const Messagetype Message_get_coupling_model_from_agent::get_message_type() const { return Messagetype::GET_COUPLING_MODEL_FROM_AGENT; }

	Message_send_convergenceFlag::Message_send_convergenceFlag() {}
	Message_send_convergenceFlag::Message_send_convergenceFlag(const bool converged, const int from)
		: converged_(converged), from_(from) {}
	const Messagetype Message_send_convergenceFlag::get_message_type() const { return Messagetype::SEND_CONVERGENCE_FLAG; }

	Message_send_optimizationInfo::Message_send_optimizationInfo() {}
	Message_send_optimizationInfo::Message_send_optimizationInfo(const OptimizationInfo& optimization_info)
		: optimization_info_(std::make_shared<OptimizationInfo>(optimization_info)) {}
	const Messagetype Message_send_optimizationInfo::get_message_type() const { return Messagetype::SEND_OPTIMIZATION_INFO; }

	const Messagetype Message_get_number_of_active_couplings::get_message_type() const { return Messagetype::GET_NUMBER_OF_ACTIVE_COUPLINGS; }

	Message_trigger_step::Message_trigger_step() {}
	Message_trigger_step::Message_trigger_step(const AlgStep& step)
		: step_(std::make_shared<AlgStep>(step)) {}
	const Messagetype Message_trigger_step::get_message_type() const { return Messagetype::TRIGGER_STEP; }

	Message_send_simulated_state::Message_send_simulated_state() {}
	Message_send_simulated_state::Message_send_simulated_state(const std::vector<typeRNum>& new_state, ctypeRNum dt, ctypeRNum t0, ctypeRNum cost)
		: new_state_(new_state), dt_(dt), t0_(t0), cost_(cost) {}
	const Messagetype Message_send_simulated_state::get_message_type() const { return Messagetype::SEND_SIMULATED_STATE; }

	const Messagetype Message_get_solution::get_message_type() const { return Messagetype::GET_SOLUTION; }

	const Messagetype Message_send_flag_to_agents::get_message_type() const { return Messagetype::SEND_FLAG_TO_AGENTS; }

	const Messagetype Message_send_ping::get_message_type() const { return Messagetype::SEND_PING; }

	Message_send_number_of_active_couplings::Message_send_number_of_active_couplings() {}
	Message_send_number_of_active_couplings::Message_send_number_of_active_couplings(const int number)
		: number_(number) {}
	const Messagetype Message_send_number_of_active_couplings::get_message_type() const { return Messagetype::SEND_NUMBER_OF_ACTIVE_COUPLINGS; }

	const Messagetype Message_acknowledge_executed_ADMM_step::get_message_type() const { return Messagetype::ACKNOWELDGE_EXECUTED_ADMM_STEP; }

	Message_send_coupling_models::Message_send_coupling_models() {}
	Message_send_coupling_models::Message_send_coupling_models(const std::shared_ptr<std::map<int, CouplingModelPtr>>& coupling_models, const int from)
		: coupling_models_(coupling_models), from_(from) {}
	const Messagetype Message_send_coupling_models::get_message_type() const { return Messagetype::SEND_COUPLING_MODELS; }

	Message_send_agent_model::Message_send_agent_model() {}
	Message_send_agent_model::Message_send_agent_model(const AgentModelPtr& agent_model, const int from)
		: agent_model_(agent_model), from_(from) {}
	const Messagetype Message_send_agent_model::get_message_type() const { return Messagetype::SEND_AGENT_MODEL; }

	Message_send_agent_state_for_simulation::Message_send_agent_state_for_simulation() {}
	Message_send_agent_state_for_simulation::Message_send_agent_state_for_simulation(const AgentState& agent_state, const int from)
		: agent_state_(std::make_shared<AgentState>(agent_state)), from_(from) {}
	const Messagetype Message_send_agent_state_for_simulation::get_message_type() const { return Messagetype::SEND_AGENT_STATE_FOR_SIMULATION; }

	Message_successfully_deregistered_agent::Message_successfully_deregistered_agent() {}
	Message_successfully_deregistered_agent::Message_successfully_deregistered_agent(const AgentInfo& agent_info)
		: agent_info_(std::make_shared<AgentInfo>(agent_info)) {}
	const Messagetype Message_successfully_deregistered_agent::get_message_type() const { return Messagetype::SUCCESSFULLY_DEREGISTERED_AGENT; }

	Message_send_solution::Message_send_solution() {}
	Message_send_solution::Message_send_solution(const SolutionPtr& solution)
		: solution_(solution) {}
	const Messagetype Message_send_solution::get_message_type() const { return Messagetype::SEND_SOLUTION; }

	const Messagetype Message_acknowledge_received_optimization_info::get_message_type() const { return Messagetype::ACKNOWLEDGE_RECEIVED_OPTIMIZATION_INFO; }

	Message_successfully_registered_agent::Message_successfully_registered_agent() {}
	Message_successfully_registered_agent::Message_successfully_registered_agent(const AgentInfo& agent_info)
		: agent_info_(std::make_shared<AgentInfo>(agent_info)) {}
	const Messagetype Message_successfully_registered_agent::get_message_type() const { return Messagetype::SUCCESSFULLY_REGISTERED_AGENT; }

	Message_successfully_registered_coupling::Message_successfully_registered_coupling() {}
	Message_successfully_registered_coupling::Message_successfully_registered_coupling(const CouplingInfo& coupling_info, const AgentInfo& agent_info, const CommunicationInfo& communication_info)
		: coupling_info_(std::make_shared<CouplingInfo>(coupling_info)),
		agent_info_(std::make_shared<AgentInfo>(agent_info)),
		communication_info_(std::make_shared<CommunicationInfo>(communication_info)) {}
	const Messagetype Message_successfully_registered_coupling::get_message_type() const { return Messagetype::SUCCESSFULLY_REGISTERED_COUPLING; }

	Message_send_flag_stopped_admm::Message_send_flag_stopped_admm() {}
	Message_send_flag_stopped_admm::Message_send_flag_stopped_admm(const bool flag_stoppedAdmm, const int from)
		:flag_stoppedAdmm_(flag_stoppedAdmm),
		 from_(from) {}
	const Messagetype Message_send_flag_stopped_admm::get_message_type() const { return Messagetype::SEND_FLAG_STOPPED_ADMM; }

	Message_send_flag_stopped_admm_coordinator::Message_send_flag_stopped_admm_coordinator() {}
	Message_send_flag_stopped_admm_coordinator::Message_send_flag_stopped_admm_coordinator(const bool flag_stoppedAdmm, const int from)
		: flag_stoppedAdmm_(flag_stoppedAdmm),
		  from_(from) {}
	const Messagetype Message_send_flag_stopped_admm_coordinator::get_message_type() const { return Messagetype::SEND_FLAG_STOPPED_ADMM_COORDINATOR; }

	Message_send_flag_to_stop_admm::Message_send_flag_to_stop_admm() {}
	Message_send_flag_to_stop_admm::Message_send_flag_to_stop_admm(const bool flag_toStopAdmm)
		: flag_toStopAdmm_(flag_toStopAdmm){}
	const Messagetype Message_send_flag_to_stop_admm::get_message_type() const { return Messagetype::SEND_FLAG_TO_STOP_ADMM; }

}