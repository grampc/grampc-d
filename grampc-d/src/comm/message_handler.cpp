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

#include "grampcd/comm/message_handler.hpp"
#include "grampcd/comm/message.hpp"
#include "grampcd/comm/communication_interface_local.hpp"

#include "grampcd/util/logging.hpp"

namespace grampcd
{
	MessageHandler::MessageHandler(CommunicationInterfaceLocal* communication_interface, const LoggingPtr& log)
		: communication_interface_(communication_interface), log_(log) {}

	void MessageHandler::handle_message(const CommunicationDataPtr& comm_data, const MessagePtr& message)
	{
		const auto message_type = message->get_message_type();

		switch (message_type)
		{
		case Messagetype::SEND_LOCAL_COPIES:
		{
			const auto my_message = std::static_pointer_cast<Message_send_local_copies>(message);
			communication_interface_->fromCommunication_send_localCopies(comm_data, *my_message->agent_state_, my_message->from_);

			break;
		}
		case Messagetype::SEND_COUPLING_STATE:
		{
			const auto my_message = std::static_pointer_cast<Message_send_coupling_state>(message);
			communication_interface_->fromCommunication_send_couplingState(comm_data, *my_message->coupling_state_, my_message->from_);

			break;
		}
		case Messagetype::SEND_TWO_COUPLING_STATES:
		{
			const auto my_message = std::static_pointer_cast<Message_send_two_coupling_states>(message);
			communication_interface_->fromCommunication_send_couplingState(comm_data, *my_message->state1_, 
				*my_message->state2_, my_message->from_);

			break;
		}
		case Messagetype::SEND_MULTIPLIER_STATE:
		{
			const auto my_message = std::static_pointer_cast<Message_send_multiplier_state>(message);
			communication_interface_->fromCommunication_send_multiplierPenaltyState(comm_data, *my_message->multiplier_state_, 
				*my_message->penalty_state_, my_message->from_);

			break;
		}
		case Messagetype::SEND_AGENT_STATE:
		{
			const auto my_message = std::static_pointer_cast<Message_send_agent_state>(message);
			communication_interface_->fromCommunication_send_agentState(comm_data, *my_message->agent_state_, *my_message->constr_state_, my_message->from_);

			break;
		}
		case Messagetype::SEND_DESIRED_AGENT_STATE:
		{
			const auto my_message = std::static_pointer_cast<Message_send_desired_agent_state>(message);
			communication_interface_->fromCommunication_send_desiredAgentState(comm_data, my_message->desired_agent_state_);

			break;
		}
		case Messagetype::FROM_COMMUNICATION_DEREGISTERED_COUPLING:
		{
			const auto my_message = std::static_pointer_cast<Message_fromCommunication_deregistered_coupling>(message);
			communication_interface_->fromCommunication_deregistered_coupling(comm_data, *my_message->coupling_info_);

			break;
		}
		case Messagetype::GET_PING:
		{
			communication_interface_->fromCommunication_get_ping(comm_data);

			break;
		}
		case Messagetype::REGISTER_AGENT:
		{
			const auto my_message = std::static_pointer_cast<Message_register_agent>(message);
			communication_interface_->fromCommunication_register_agent(comm_data, my_message->agent_info_);

			break;
		}
		case Messagetype::DEREGISTER_AGENT:
		{
			const auto my_message = std::static_pointer_cast<Message_deregister_agent>(message);
			communication_interface_->fromCommunication_deregister_agent(comm_data, *my_message->agent_info_);

			break;
		}
		case Messagetype::REGISTER_COUPLING:
		{
			const auto my_message = std::static_pointer_cast<Message_register_coupling>(message);
			communication_interface_->fromCommunication_register_coupling(comm_data, *my_message->coupling_info_);

			break;
		}
		case Messagetype::DEREGISTER_COUPLING:
		{
			const auto my_message = std::static_pointer_cast<Message_deregister_coupling>(message);
			communication_interface_->fromCommunication_deregister_coupling(comm_data, my_message->coupling_info_);

			break;
		}
		case Messagetype::SEND_NUMBER_OF_NEIGHBORS:
		{
			const auto my_message = std::static_pointer_cast<Message_send_numberOfNeighbors>(message);
			communication_interface_->fromCommunication_send_numberOfNeighbors(comm_data, my_message->number_, my_message->from_);

			break;
		}
		case Messagetype::GET_AGENT_STATE_FOR_SIMULATION:
		{
			const auto my_message = std::static_pointer_cast<Message_get_agentState_for_simulation>(message);
			communication_interface_->fromCommunication_get_agentState_for_simulation(comm_data);

			break;
		}
		case Messagetype::GET_DESIRED_AGENT_STATE_FROM_AGENT:
		{
			communication_interface_->fromCommunication_get_desiredAgentStateFromAgent(comm_data);

			break;
		}
		case Messagetype::GET_AGENT_MODEL_FROM_AGENT:
		{
			const auto my_message = std::static_pointer_cast<Message_get_agent_model_from_agent>(message);
			communication_interface_->fromCommunication_get_agentModelFromAgent(comm_data);

			break;
		}
		case Messagetype::GET_COUPLING_MODEL_FROM_AGENT:
		{
			const auto my_message = std::static_pointer_cast<Message_get_coupling_model_from_agent>(message);
			communication_interface_->fromCommunication_get_couplingModelFromAgent(comm_data);

			break;
		}
		case Messagetype::SEND_CONVERGENCE_FLAG:
		{
			const auto my_message = std::static_pointer_cast<Message_send_convergenceFlag>(message);
			communication_interface_->fromCommunication_send_convergenceFlag(comm_data, my_message->converged_, my_message->from_);

			break;
		}
		case Messagetype::SEND_OPTIMIZATION_INFO:
		{
			const auto my_message = std::static_pointer_cast<Message_send_optimizationInfo>(message);
			communication_interface_->fromCommunication_send_optimizationInfo(comm_data, *my_message->optimization_info_);

			break;
		}
		case Messagetype::GET_NUMBER_OF_ACTIVE_COUPLINGS:
		{
			communication_interface_->fromCommunication_get_numberOfActiveCouplings(comm_data);

			break;
		}
		case Messagetype::TRIGGER_STEP:
		{
			const auto my_message = std::static_pointer_cast<Message_trigger_step>(message);
			communication_interface_->fromCommunication_triggerStep(*my_message->step_);

			break;
		}
		case Messagetype::SEND_SIMULATED_STATE:
		{
			const auto my_message = std::static_pointer_cast<Message_send_simulated_state>(message);

			communication_interface_->fromCommunication_send_simulatedState(
				comm_data, my_message->new_state_, my_message->dt_,
				my_message->t0_, my_message->cost_);

			break;
		}
		case Messagetype::GET_SOLUTION:
		{
			communication_interface_->fromCommunication_get_solution(comm_data);
			break;
		}
		case Messagetype::SEND_FLAG_TO_AGENTS:
		{
			communication_interface_->fromCommunication_send_flagToAgents(comm_data);
			break;
		}
		case Messagetype::SEND_PING:
		{
			communication_interface_->fromCommunication_send_ping(comm_data);

			break;
		}
		case Messagetype::SEND_NUMBER_OF_ACTIVE_COUPLINGS:
		{
			const auto my_message = std::static_pointer_cast<Message_send_number_of_active_couplings>(message);
			communication_interface_->fromCommunication_send_numberOfActiveCouplings(comm_data, my_message->number_);

			break;
		}
		case Messagetype::ACKNOWELDGE_EXECUTED_ADMM_STEP:
		{
			communication_interface_->fromCommunication_received_acknowledgement_executed_AlgStep(comm_data);
			break;
		}
		case Messagetype::SEND_COUPLING_MODELS:
		{
			const auto my_message = std::static_pointer_cast<Message_send_coupling_models>(message);
			communication_interface_->fromCommunication_send_couplingModel(comm_data, my_message->coupling_models_);

			break;
		}
		case Messagetype::SEND_AGENT_MODEL:
		{
			const auto my_message = std::static_pointer_cast<Message_send_agent_model>(message);
			communication_interface_->fromCommunication_send_agentModel(comm_data, my_message->agent_model_);

			break;
		}
		case Messagetype::SEND_AGENT_STATE_FOR_SIMULATION:
		{
			const auto my_message = std::static_pointer_cast<Message_send_agent_state_for_simulation>(message);
			communication_interface_->fromCommunication_send_agentState_for_simulation(comm_data, my_message->agent_state_);

			break;
		}
		case Messagetype::SEND_SOLUTION:
		{
			const auto my_message = std::static_pointer_cast<Message_send_solution>(message);
			communication_interface_->fromCommunication_send_solution(comm_data, my_message->solution_);

			break;
		}
		case Messagetype::ACKNOWLEDGE_RECEIVED_OPTIMIZATION_INFO:
		{
			communication_interface_->fromCommunication_received_acknowledgement_received_optimizationInfo(comm_data);
			break;
		}
		case Messagetype::SUCCESSFULLY_REGISTERED_AGENT:
		{
			const auto my_message = std::static_pointer_cast<Message_successfully_registered_agent>(message);
			communication_interface_->fromCommunication_successfully_registered_agent(comm_data, *my_message->agent_info_);

			break;
		}
		case Messagetype::SUCCESSFULLY_REGISTERED_COUPLING:
		{
			const auto my_message = std::static_pointer_cast<Message_successfully_registered_coupling>(message);

			communication_interface_->fromCommunication_successfully_registered_coupling(comm_data,
				*my_message->coupling_info_, *my_message->agent_info_, my_message->communication_info_);

			break;
		}
		case Messagetype::SEND_COMMUNICATION_INFO:
		{
			const auto my_message = std::static_pointer_cast<Message_send_communication_info>(message);
			communication_interface_->fromCommunication_send_communicationInfo(comm_data, my_message->communication_info_);

			break;
		}
		case Messagetype::SUCCESSFULLY_DEREGISTERED_AGENT:
		{
			const auto my_message = std::static_pointer_cast<Message_successfully_deregistered_agent>(message);
			communication_interface_->fromCommunication_successfully_deregistered_agent(comm_data, *my_message->agent_info_);

			break;
		}

		case Messagetype::SEND_FLAG_STOPPED_ADMM:
		{
			const auto my_message = std::static_pointer_cast<Message_send_flag_stopped_admm>(message);
			communication_interface_->fromCommunication_send_stoppedAlgFlag(comm_data,my_message->flag_stoppedAdmm_, my_message->from_);

			break;
		}

		case Messagetype::SEND_FLAG_STOPPED_ADMM_COORDINATOR:
		{
			const auto my_message = std::static_pointer_cast<Message_send_flag_stopped_admm_coordinator>(message);
			communication_interface_->fromCommunication_send_stoppedAlgFlagCoordinator(comm_data, my_message->flag_stoppedAdmm_, my_message->from_);

			break;
		}

		case Messagetype::SEND_FLAG_TO_STOP_ADMM:
		{
			const auto my_message = std::static_pointer_cast<Message_send_flag_to_stop_admm>(message);
			communication_interface_->fromCommunication_send_flagToStopAdmm(comm_data, my_message->flag_toStopAdmm_);

			break;
		}


		default:
			log_->print(DebugType::Error) << "[MessageHandler::handle_message] Unknown message type." << std::endl;
			break;
		}
	}
}