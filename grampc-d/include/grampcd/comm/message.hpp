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

#include "grampcd/util/class_forwarding.hpp"

#include <cereal/types/polymorphic.hpp>

namespace grampcd
{
	enum class Messagetype
	{
		SEND_LOCAL_COPIES,
		SEND_COUPLING_STATE,
		SEND_TWO_COUPLING_STATES,
		SEND_MULTIPLIER_STATE,
		SEND_AGENT_STATE,
		SEND_DESIRED_AGENT_STATE,
		FROM_COMMUNICATION_DEREGISTERED_COUPLING,
		GET_PING,
		REGISTER_AGENT,
		DEREGISTER_AGENT,
		REGISTER_COUPLING,
		DEREGISTER_COUPLING,
		SEND_NUMBER_OF_NEIGHBORS,
		GET_AGENT_STATE_FOR_SIMULATION,
		GET_DESIRED_AGENT_STATE_FROM_AGENT,
		GET_AGENT_MODEL_FROM_AGENT,
		GET_COUPLING_MODEL_FROM_AGENT,
		SEND_CONVERGENCE_FLAG,
		SEND_OPTIMIZATION_INFO,
		GET_NUMBER_OF_ACTIVE_COUPLINGS,
		TRIGGER_STEP,
		SEND_SIMULATED_STATE,
		GET_SOLUTION,
		SEND_FLAG_TO_AGENTS,
		SEND_PING,
		SEND_NUMBER_OF_ACTIVE_COUPLINGS,
		ACKNOWELDGE_EXECUTED_ADMM_STEP,
		SEND_COUPLING_MODELS,
		SEND_AGENT_MODEL,
		SEND_AGENT_STATE_FOR_SIMULATION,
		SEND_SOLUTION,
		ACKNOWLEDGE_RECEIVED_OPTIMIZATION_INFO,
		SUCCESSFULLY_REGISTERED_AGENT,
		SUCCESSFULLY_REGISTERED_COUPLING,
		SEND_COMMUNICATION_INFO,
		SUCCESSFULLY_DEREGISTERED_AGENT,
		SEND_FLAG_STOPPED_ADMM,
		SEND_FLAG_STOPPED_ADMM_COORDINATOR,
		SEND_FLAG_TO_STOP_ADMM
	};

	struct Message
	{
		virtual const Messagetype get_message_type() const = 0;
	};

	struct Message_send_communication_info : public Message
	{
		Message_send_communication_info();
		Message_send_communication_info(const CommunicationInfo communication_info);
		const Messagetype get_message_type() const;

		CommunicationInfoPtr communication_info_;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(communication_info_);
		}
	};

	struct Message_fromCommunication_deregistered_coupling : public Message
	{
		Message_fromCommunication_deregistered_coupling();
		Message_fromCommunication_deregistered_coupling(const CouplingInfo coupling_info);
		const Messagetype get_message_type() const;

		CouplingInfoPtr coupling_info_;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(coupling_info_);
		}
	};

	struct Message_get_ping : public Message
	{
		const Messagetype get_message_type() const;

		template<class Archive>
		void serialize(Archive& ar)
		{}
	};

	struct Message_register_agent : public Message
	{
		Message_register_agent();
		Message_register_agent(const AgentInfo& agent_info);
		const Messagetype get_message_type() const;

		AgentInfoPtr agent_info_;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(agent_info_);
		}
	};

	struct Message_deregister_agent : public Message
	{
		Message_deregister_agent();
		Message_deregister_agent(const AgentInfo& agent_info);
		const Messagetype get_message_type() const;

		AgentInfoPtr agent_info_;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(agent_info_);
		}
	};

	struct Message_register_coupling : public Message
	{
		Message_register_coupling();
		Message_register_coupling(const CouplingInfo& coupling_info);
		const Messagetype get_message_type() const;

		CouplingInfoPtr coupling_info_;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(coupling_info_);
		}
	};

	struct Message_deregister_coupling : public Message
	{
		Message_deregister_coupling();
		Message_deregister_coupling(const CouplingInfo& coupling_info);
		const Messagetype get_message_type() const;

		CouplingInfoPtr coupling_info_;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(coupling_info_);
		}
	};

	struct Message_send_numberOfNeighbors : public Message
	{
		Message_send_numberOfNeighbors();
		Message_send_numberOfNeighbors(const int number, const int from);
		const Messagetype get_message_type() const;

		int number_ = 0, from_ = 0;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(number_, from_);
		}
	};

	struct Message_send_local_copies : public Message
	{
		Message_send_local_copies();
		Message_send_local_copies(const AgentState& agent_state, const int from);
		const Messagetype get_message_type() const;

		AgentStatePtr agent_state_;
		int from_ = 0;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(agent_state_, from_);
		}
	};

	struct Message_send_desired_agent_state : public Message
	{
		Message_send_desired_agent_state();
		Message_send_desired_agent_state(const AgentState& desired_agent_state, const int from);
		const Messagetype get_message_type() const;

		AgentStatePtr desired_agent_state_;
		int from_ = 0;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(desired_agent_state_, from_);
		}
	};

	struct Message_send_coupling_state : public Message
	{
		Message_send_coupling_state();
		Message_send_coupling_state(const CouplingState& coupling_state, const int from);
		const Messagetype get_message_type() const;

		CouplingStatePtr coupling_state_;
		int from_ = 0;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(coupling_state_, from_);
		}
	};

	struct Message_send_two_coupling_states : public Message
	{
		Message_send_two_coupling_states();
		Message_send_two_coupling_states(const CouplingState& state1, const CouplingState& state2, const int from);
		const Messagetype get_message_type() const;

		CouplingStatePtr state1_;
		CouplingStatePtr state2_;
		int from_ = 0;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(state1_, state2_, from_);
		}
	};

	struct Message_send_multiplier_state : public Message
	{
		Message_send_multiplier_state();
		Message_send_multiplier_state(const MultiplierState& multiplier_state, const PenaltyState& penalty_state, const int from);
		const Messagetype get_message_type() const;

		MultiplierStatePtr multiplier_state_;
		PenaltyStatePtr penalty_state_;
		int from_ = 0;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(multiplier_state_, penalty_state_, from_);
		}
	};

	struct Message_send_agent_state : public Message
	{
		Message_send_agent_state();
		Message_send_agent_state(const AgentState& agent_state, const ConstraintState& constr_state, const int from);
		const Messagetype get_message_type() const;

		AgentStatePtr agent_state_;
		ConstraintStatePtr constr_state_;
		int from_ = 0;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(agent_state_, constr_state_, from_);
		}
	};

	struct Message_get_agentState_for_simulation : public Message
	{
		Message_get_agentState_for_simulation();
		const Messagetype get_message_type() const;

		template<class Archive>
		void serialize(Archive& ar)
		{}
	};

	struct Message_get_desired_agent_state_from_agent : public Message
	{
		const Messagetype get_message_type() const;

		template<class Archive>
		void serialize(Archive& ar)
		{}
	};

	struct Message_get_agent_model_from_agent : public Message
	{
		const Messagetype get_message_type() const;

		template<class Archive>
		void serialize(Archive& ar)
		{}
	};

	struct Message_get_coupling_model_from_agent : public Message
	{
		const Messagetype get_message_type() const;

		template<class Archive>
		void serialize(Archive& ar)
		{}
	};

	struct Message_send_convergenceFlag : public Message
	{
		Message_send_convergenceFlag();
		Message_send_convergenceFlag(const bool converged, const int from);
		const Messagetype get_message_type() const;

		bool converged_ = false;
		int from_ = 0;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(converged_, from_);
		}
	};

	struct Message_send_optimizationInfo : public Message
	{
		Message_send_optimizationInfo();
		Message_send_optimizationInfo(const OptimizationInfo& optimization_info);
		const Messagetype get_message_type() const;

		OptimizationInfoPtr optimization_info_;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(optimization_info_);
		}
	};

	struct Message_get_number_of_active_couplings : public Message
	{
		const Messagetype get_message_type() const;

		template<class Archive>
		void serialize(Archive& ar)
		{}
	};

	struct Message_trigger_step : public Message
	{
		Message_trigger_step();
		Message_trigger_step(const AlgStep& step);
		const Messagetype get_message_type() const;

		std::shared_ptr<AlgStep> step_;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(step_);
		}
	};

	struct Message_send_simulated_state : public Message
	{
		Message_send_simulated_state();
		Message_send_simulated_state(const std::vector<typeRNum>& new_state, ctypeRNum dt, ctypeRNum t0, ctypeRNum cost);
		const Messagetype get_message_type() const;

		std::vector<typeRNum> new_state_;
		typeRNum dt_ = 0;
		typeRNum t0_ = 0;
		typeRNum cost_ = 0;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(new_state_, dt_, t0_, cost_);
		}
	};

	struct Message_get_solution : public Message
	{
		const Messagetype get_message_type() const;

		template<class Archive>
		void serialize(Archive& ar)
		{}
	};

	struct Message_send_flag_to_agents : public Message
	{
		const Messagetype get_message_type() const;

		template<class Archive>
		void serialize(Archive& ar)
		{}
	};

	struct Message_send_ping : public Message
	{
		const Messagetype get_message_type() const;

		template<class Archive>
		void serialize(Archive& ar)
		{}
	};

	struct Message_send_number_of_active_couplings : public Message
	{
		Message_send_number_of_active_couplings();
		Message_send_number_of_active_couplings(const int number);
		const Messagetype get_message_type() const;

		int number_ = 0;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(number_);
		}
	};

	struct Message_acknowledge_executed_ADMM_step : public Message
	{
		const Messagetype get_message_type() const;

		template<class Archive>
		void serialize(Archive& ar)
		{}
	};

	struct Message_send_coupling_models : public Message
	{
		Message_send_coupling_models();
		Message_send_coupling_models(const std::shared_ptr<std::map<int, CouplingModelPtr>>& coupling_models, const int from);
		const Messagetype get_message_type() const;

		std::shared_ptr<std::map<int, CouplingModelPtr>> coupling_models_;
		int from_ = 0;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(coupling_models_, from_);
		}
	};

	struct Message_send_agent_model : public Message
	{
		Message_send_agent_model();
		Message_send_agent_model(const AgentModelPtr& agent_model, const int from);
		const Messagetype get_message_type() const;

		AgentModelPtr agent_model_;
		int from_ = 0;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(agent_model_, from_);
		}
	};

	struct Message_send_agent_state_for_simulation : public Message
	{
		Message_send_agent_state_for_simulation();
		Message_send_agent_state_for_simulation(const AgentState& agent_state, const int from);
		const Messagetype get_message_type() const;

		AgentStatePtr agent_state_;
		int from_ = 0;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(agent_state_, from_);
		}
	};

	struct Message_successfully_deregistered_agent : public Message
	{
		Message_successfully_deregistered_agent();
		Message_successfully_deregistered_agent(const AgentInfo& agent_info);
		const Messagetype get_message_type() const;

		AgentInfoPtr agent_info_;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(agent_info_);
		}
	};

	struct Message_send_solution : public Message
	{
		Message_send_solution();
		Message_send_solution(const SolutionPtr& solution);
		const Messagetype get_message_type() const;

		SolutionPtr solution_;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(solution_);
		}
	};

	struct Message_acknowledge_received_optimization_info : public Message
	{
		const Messagetype get_message_type() const;

		template<class Archive>
		void serialize(Archive& ar) {}
	};

	struct Message_successfully_registered_agent : public Message
	{
		Message_successfully_registered_agent();
		Message_successfully_registered_agent(const AgentInfo& agent_info);
		const Messagetype get_message_type() const;

		AgentInfoPtr agent_info_;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(agent_info_);
		}
	};

	struct Message_successfully_registered_coupling : public Message
	{
		Message_successfully_registered_coupling();
		Message_successfully_registered_coupling(const CouplingInfo& coupling_info, const AgentInfo& agent_info, const CommunicationInfo& communication_info);
		const Messagetype get_message_type() const;

		CouplingInfoPtr coupling_info_;
		AgentInfoPtr agent_info_;
		CommunicationInfoPtr communication_info_;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(coupling_info_, agent_info_, communication_info_);
		}
	};

	
	struct Message_send_flag_stopped_admm : public Message
	{
		Message_send_flag_stopped_admm();
		Message_send_flag_stopped_admm(const bool flag_stoppedAdmm, const int from);
		const Messagetype get_message_type() const;
		
		 bool flag_stoppedAdmm_ = false;
		 int from_ = 0;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(flag_stoppedAdmm_,from_);
		}
	};

	struct Message_send_flag_stopped_admm_coordinator : public Message
	{
		Message_send_flag_stopped_admm_coordinator();
		Message_send_flag_stopped_admm_coordinator(const bool flag_stoppedAdmm, const int from);
		const Messagetype get_message_type() const;

		bool flag_stoppedAdmm_ = false;
		int from_ = 0;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(flag_stoppedAdmm_, from_);
		}
	};

	struct Message_send_flag_to_stop_admm : public Message
	{
		Message_send_flag_to_stop_admm();
		Message_send_flag_to_stop_admm(const bool flag_toStopAdmm);
		const Messagetype get_message_type() const;

		bool flag_toStopAdmm_ = false;

		template<class Archive>
		void serialize(Archive& ar)
		{
			ar(flag_toStopAdmm_);
		}
	};


}

CEREAL_FORCE_DYNAMIC_INIT(message)