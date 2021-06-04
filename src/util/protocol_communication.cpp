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

#include "dmpc/util/protocol_communication.hpp"
#include "dmpc/util/data_conversion.hpp"
#include "dmpc/util/logging.hpp"

#include "dmpc/model/agent_model.hpp"
#include "dmpc/model/coupling_model.hpp"

#include "dmpc/comm/communication_interface_local.hpp"

#include "dmpc/info/agent_info.hpp"
#include "dmpc/info/coupling_info.hpp"
#include "dmpc/info/optimization_info.hpp"

#include "dmpc/state/agent_state.hpp"
#include "dmpc/state/coupling_state.hpp"
#include "dmpc/state/multiplier_state.hpp"
#include "dmpc/state/penalty_state.hpp"

#include "dmpc/optim/solution.hpp"

#include "general_model_factory.hpp"

namespace dmpc
{
    void ProtocolCommunication::evaluateData(CommunicationDataPtr comm_data, CommunicationInterfaceLocal* communication_interface, const std::vector<char>& data)
    {
        const index dataType = static_cast<index>(data[position_of_index_in_protocol_]);
        switch (dataType)
        {
        case index::successfully_registered_agent:
            communication_interface->fromCommunication_successfully_registered_agent(comm_data, ProtocolCommunication::buildFromProtocol_agentInfo(data));
            break;
        case index::received_acknowledgement_received_optimizationInfo:
            communication_interface->fromCommunication_received_acknowledgement_received_optimizationInfo(comm_data);
            break;
        case index::received_acknowledgement_executed_ADMMstep:
            communication_interface->fromCommunication_received_acknowledgement_executed_ADMMstep(comm_data);
            break;
        case index::get_solution:
            communication_interface->fromCommunication_get_solution(comm_data);
            break;
        case index::send_solution:
            communication_interface->fromCommunication_send_solution(comm_data, ProtocolCommunication::buildFromProtocol_solution(data));
            break;
        case index::register_agent:
            communication_interface->fromCommunication_register_agent(comm_data, ProtocolCommunication::buildFromProtocol_agentInfo(data));
            break;
        case index::register_coupling:
            communication_interface->fromCommunication_register_coupling(comm_data, ProtocolCommunication::buildFromProtocol_couplingInfo(data));
            break;
        case index::successfully_registered_coupling:
            communication_interface->fromCommunication_successfully_registered_coupling(comm_data,
                ProtocolCommunication::buildFromProtocol_couplingInfo_from_successfully_registered_coupling(data),
                ProtocolCommunication::buildFromProtocol_agentInfo_from_successfully_registered_coupling(data),
                ProtocolCommunication::buildFromProtocol_communicationInfo_from_successfully_registered_coupling(data));
            break;
        case index::deregister_agent:
            communication_interface->fromCommunication_deregister_agent(comm_data, ProtocolCommunication::buildFromProtocol_agentInfo(data));
            break;
        case index::deregister_coupling:
            communication_interface->fromCommunication_deregister_coupling(comm_data, ProtocolCommunication::buildFromProtocol_couplingInfo(data));
            break;
        case index::send_optimizationInfo:
            communication_interface->fromCommunication_send_optimizationInfo(comm_data, ProtocolCommunication::buildFromProtocol_optimizationInfo(data));
            break;
        case index::send_numberOfNeighbors:
            communication_interface->fromCommunication_send_numberOfNeighbors(comm_data, ProtocolCommunication::buildFromProtocol_numberOfNeighbors(data),
                ProtocolCommunication::buildFromProtocol_from(data));
            break;
        case index::deregistered_coupling:
            communication_interface->fromCommunication_deregistered_coupling(comm_data, ProtocolCommunication::buildFromProtocol_couplingInfo(data));
            break;
        case index::send_agentState:
            communication_interface->fromCommunication_send_agentState(comm_data, ProtocolCommunication::buildFromProtocol_agentState(data),
                ProtocolCommunication::buildFromProtocol_from(data));
            break;
        case index::send_desiredAgentState:
            communication_interface->fromCommunication_send_desiredAgentState(comm_data, ProtocolCommunication::buildFromProtocol_agentState(data));
            break;
        case index::send_couplingState:
            communication_interface->fromCommunication_send_couplingState(
                comm_data,
                ProtocolCommunication::buildFromProtocol_couplingState(data),
                ProtocolCommunication::buildFromProtocol_from(data));
            break;
        case index::send_two_couplingStates:
            communication_interface->fromCommunication_send_couplingState(
                comm_data,
                ProtocolCommunication::buildFromProtocol_couplingState_one(data),
                ProtocolCommunication::buildFromProtocol_couplingState_two(data),
                ProtocolCommunication::buildFromProtocol_from(data));
            break;
        case index::send_multiplierPenaltyState:
            communication_interface->fromCommunication_send_multiplierPenaltyState(comm_data, ProtocolCommunication::buildFromProtocol_multiplierState(data),
                ProtocolCommunication::buildFromProtocol_penaltyState(data),
                ProtocolCommunication::buildFromProtocol_from(data));
            break;
        case index::get_agentState_for_simulation:
            communication_interface->fromCommunication_get_agentState_for_simulation(comm_data);
            break;
        case index::send_convergenceFlag:
            communication_interface->fromCommunication_send_convergenceFlag(comm_data, ProtocolCommunication::buildFromProtocol_convergenceFlag(data),
                ProtocolCommunication::buildFromProtocol_from(data));
            break;
        case index::send_agentState_for_simulation:
            communication_interface->fromCommunication_send_agentState_for_simulation(comm_data, ProtocolCommunication::buildFromProtocol_agentState(data));
            break;
        case index::triggerStep:
            communication_interface->fromCommunication_triggerStep(ProtocolCommunication::buildFromProtocol_ADMMStep(data));
            break;
        case index::get_agentModelFromAgent:
            communication_interface->fromCommunication_get_agentModelFromAgent(comm_data);
            break;
        case index::send_agentModel:
            communication_interface->fromCommunication_send_agentModel(comm_data, ProtocolCommunication::buildFromProtocol_agentModel(data, communication_interface->get_log()));
            break;
        case index::get_couplingModelFromAgent:
            communication_interface->fromCommunication_get_couplingModelFromAgent(comm_data);
            break;
        case index::send_couplingModel:
            communication_interface->fromCommunication_send_couplingModel(comm_data, ProtocolCommunication::buildFromProtocol_couplingModel(data, communication_interface->get_log()));
            break;
        case index::send_simulatedState:
            communication_interface->fromCommunication_send_simulatedState(comm_data, ProtocolCommunication::buildFromProtocol_xnext_from_simulatedState(data),
                ProtocolCommunication::buildFromProtocol_dt_from_simulatedState(data),
                ProtocolCommunication::buildFromProtocol_t0_from_simulatedState(data));
            break;
        case index::get_desiredAgentStateFromAgent:
            communication_interface->fromCommunication_get_desiredAgentStateFromAgent(comm_data);
            break;
        case index::get_ping:
            communication_interface->fromCommunication_get_ping(comm_data);
            break;
        case index::send_ping:
            communication_interface->fromCommunication_send_ping(comm_data);
            break;
        case index::get_numberOfActiveCouplings:
            communication_interface->fromCommunication_get_numberOfActiveCouplings(comm_data);
            break;
        case index::send_numberOfActiveCouplings:
            communication_interface->fromCommunication_send_numberOfActiveCouplings(comm_data, buildFromProtocol_numberOfActiveCouplings(data));
            break;
        case index::send_flagToAgents:
            communication_interface->fromCommunication_send_flagToAgents(comm_data);
            break;
        case index::send_communicationInfo:
            communication_interface->fromCommunication_send_communicationInfo(comm_data, buildFromProtocol_communicationInfo(data));
            break;
        case index::successfully_deregistered_agent:
            communication_interface->fromCommunication_successfully_deregistered_agent(comm_data, buildFromProtocol_agentInfo(data));
            break;
        default:
            const LoggingPtr& log = communication_interface->get_log();
            log->print(DebugType::Error) << "[ProtocolCommunication::evaluateData] Unknown index in protocol." << std::endl;
            break;
        }
    }

    /****************************************************
     * functions to build protocol
    ****************************************************/

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_send_flagToAgents()
    {
        const char index = static_cast<char>(index::send_flagToAgents);
        const auto size_of_header = static_cast<char>(first_element_with_data_);
        const unsigned int size_of_data = size_of_header;
        std::shared_ptr< std::vector<char> > data(new std::vector<char>(size_of_data, 0));
        unsigned int pos = 0;

        // size of data
        DataConversion::insert_into_charArray(data, pos, size_of_data);

        // index
        DataConversion::insert_into_charArray(data, pos, index);

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_get_ping()
    {
		const char index = static_cast<char>(index::get_ping);
		const auto size_of_header = static_cast<char>(first_element_with_data_);
        const unsigned int size_of_data = size_of_header;
        std::shared_ptr< std::vector<char> > data(new std::vector<char>(size_of_data, 0));
        unsigned int pos = 0;

        // size of data
        DataConversion::insert_into_charArray(data, pos, size_of_data);

        // index
        DataConversion::insert_into_charArray(data, pos, index);

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_send_ping()
    {
		const char index = static_cast<char>(index::send_ping);
		const auto size_of_header = static_cast<char>(first_element_with_data_);
        const unsigned int size_of_data = size_of_header;
        std::shared_ptr< std::vector<char> > data(new std::vector<char>(size_of_data, 0));
        unsigned int pos = 0;

        // size of data
        DataConversion::insert_into_charArray(data, pos, size_of_data);

        // index
        DataConversion::insert_into_charArray(data, pos, index);

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_get_numberOfActiveCouplings()
    {
		const char index = static_cast<char>(index::get_numberOfActiveCouplings);
		const auto size_of_header = static_cast<char>(first_element_with_data_);
        const unsigned int size_of_data = size_of_header;

        std::shared_ptr< std::vector<char> > data(new std::vector<char>(size_of_data, 0));
        unsigned int pos = 0;

        // size of data
        DataConversion::insert_into_charArray(data, pos, size_of_data);

        // index
        DataConversion::insert_into_charArray(data, pos, index);

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_send_numberOfActiveCouplings(int number)
    {
		const char index = static_cast<char>(index::send_numberOfActiveCouplings);
		const auto size_of_header = static_cast<char>(first_element_with_data_);
        const unsigned int size_of_data = size_of_header + sizeof(int);

        std::shared_ptr< std::vector<char> > data(new std::vector<char>(size_of_data, 0));
        unsigned int pos = 0;

        // size of data
        DataConversion::insert_into_charArray(data, pos, size_of_data);

        // index
        DataConversion::insert_into_charArray(data, pos, index);

        DataConversion::insert_into_charArray(data, pos, number);

        return data;
    }

    const std::shared_ptr<std::vector<char>> ProtocolCommunication::buildProtocol_get_solution()
    {
		const char index = static_cast<char>(index::get_solution);
		const auto size_of_header = static_cast<char>(first_element_with_data_);
        const unsigned int size_of_data = size_of_header;
        std::shared_ptr<std::vector<char>> data(new std::vector<char>(size_of_data, 0));

        // include sizeOfData
        unsigned int pos = 0;
        DataConversion::insert_into_charArray(data, pos, size_of_data);

        // index
        DataConversion::insert_into_charArray(data, pos, index);

        return data;
    }

    const std::shared_ptr<std::vector<char>> ProtocolCommunication::buildProtocol_send_solution(const SolutionPtr& solution)
    {
        AgentState& agent_state = solution->agentState_;
        AgentState& predicted_agentState = solution->predicted_agentState_;
        std::vector<typeRNum>& cost = solution->cost_;

		const auto size_of_header = static_cast<char>(first_element_with_data_);
        const unsigned int sizeOfData = size_of_header // header
            + static_cast<unsigned int>(
                // agent state
                sizeof(agent_state.i_)
                + sizeof(agent_state.t0_)
                + sizeof(typeRNum) * agent_state.t_.size() + sizeof(int)
                + sizeof(typeRNum) * agent_state.u_.size() + sizeof(int)
                + sizeof(typeRNum) * agent_state.v_.size() + sizeof(int)
                + sizeof(typeRNum) * agent_state.x_.size() + sizeof(int)
                // predicted agent state
                + sizeof(predicted_agentState.i_)
                + sizeof(predicted_agentState.t0_)
                + sizeof(typeRNum) * predicted_agentState.t_.size() + sizeof(int)
                + sizeof(typeRNum) * predicted_agentState.u_.size() + sizeof(int)
                + sizeof(typeRNum) * predicted_agentState.v_.size() + sizeof(int)
                + sizeof(typeRNum) * predicted_agentState.x_.size() + sizeof(int)
                // cost 
                + sizeof(typeRNum) * cost.size() + sizeof(int));

        std::shared_ptr<std::vector<char>> data(new std::vector<char>(sizeOfData, 0));
        unsigned int pos = 0;
        const char index = static_cast<char>(index::send_solution);

        // include sizeOfData
        DataConversion::insert_into_charArray(data, pos, sizeOfData);

        // index
        DataConversion::insert_into_charArray(data, pos, index);

        // include agent state
        DataConversion::insert_into_charArray(data, pos, agent_state.i_);
        DataConversion::insert_into_charArray(data, pos, agent_state.t0_);
        DataConversion::insert_into_charArray(data, pos, agent_state.t_);
        DataConversion::insert_into_charArray(data, pos, agent_state.u_);
        DataConversion::insert_into_charArray(data, pos, agent_state.v_);
        DataConversion::insert_into_charArray(data, pos, agent_state.x_);

        // include predicted agent state
        DataConversion::insert_into_charArray(data, pos, predicted_agentState.i_);
        DataConversion::insert_into_charArray(data, pos, predicted_agentState.t0_);
        DataConversion::insert_into_charArray(data, pos, predicted_agentState.t_);
        DataConversion::insert_into_charArray(data, pos, predicted_agentState.u_);
        DataConversion::insert_into_charArray(data, pos, predicted_agentState.v_);
        DataConversion::insert_into_charArray(data, pos, predicted_agentState.x_);

        // cost
        DataConversion::insert_into_charArray(data, pos, cost);

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_successfully_deregistered_agent(const AgentInfo& info)
    {
        auto data = buildProtocol_register_agent(info);

		// change index
        (*data)[position_of_index_in_protocol_] = static_cast<char>(index::successfully_deregistered_agent);

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_successfully_registered_agent(const AgentInfo& info)
    {
        auto data = buildProtocol_register_agent(info);

		// change index
        (*data)[position_of_index_in_protocol_] = static_cast<char>(index::successfully_registered_agent);

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_get_desiredAgentState_from_agent()
    {
        const char index = static_cast<char>(index::get_desiredAgentStateFromAgent);
		unsigned int pos = 0;
		const auto size_of_header = static_cast<char>(first_element_with_data_);
        const unsigned int size_of_data = size_of_header;

        std::shared_ptr< std::vector<char> > data(new std::vector<char>(size_of_data, 0));

        // size of data
        DataConversion::insert_into_charArray(data, pos, size_of_data);

        // index
        DataConversion::insert_into_charArray(data, pos, index);

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_send_simulatedState(const std::vector<typeRNum>& x_next, typeRNum dt, typeRNum t0)
    {
        const char index = static_cast<char>(index::send_simulatedState);
		unsigned int pos = 0;
		const auto size_of_header = static_cast<char>(first_element_with_data_);
        const unsigned int size_of_data = size_of_header
            + static_cast<unsigned int>(
                // x_next
                sizeof(typeRNum) * x_next.size() + sizeof(int)
                // dt
                + sizeof(typeRNum)
                // t0
                + sizeof(typeRNum));

        std::shared_ptr< std::vector<char> > data(new std::vector<char>(size_of_data, 0));

        // size of data
        DataConversion::insert_into_charArray(data, pos, size_of_data);

        // index
        DataConversion::insert_into_charArray(data, pos, index);

        //xnext
        DataConversion::insert_into_charArray(data, pos, x_next);
        //dt
        DataConversion::insert_into_charArray(data, pos, dt);
        // t0
        DataConversion::insert_into_charArray(data, pos, t0);

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_send_agentModel(const AgentModelPtr& model, const int from)
    {
        const char index = static_cast<char>(index::send_agentModel);
		unsigned int pos = 0;
		const auto size_of_header = static_cast<char>(first_element_with_data_);
        const unsigned int size_of_data = size_of_header
            + static_cast<unsigned int>(sizeof(from)
                + model->get_modelName().size() + sizeof(int)
                + sizeof(typeRNum) * model->get_modelParameters().size() + sizeof(int)
                + sizeof(typeRNum) * model->get_costParameters().size() + sizeof(int));

        std::shared_ptr< std::vector<char> > data(new std::vector<char>(size_of_data, 0));

        // size of data
        DataConversion::insert_into_charArray(data, pos, size_of_data);

        // index
        DataConversion::insert_into_charArray(data, pos, index);

        // from
        DataConversion::insert_into_charArray(data, pos, from);

        // agent model
        DataConversion::insert_into_charArray(data, pos, model->get_modelName());
        DataConversion::insert_into_charArray(data, pos, model->get_modelParameters());
        DataConversion::insert_into_charArray(data, pos, model->get_costParameters());

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_send_couplingModel(const std::map<int, CouplingModelPtr>& map, const int from)
    {
        const char index = static_cast<char>(index::send_couplingModel);
        unsigned int pos = 0;
		// calculate size of data
		const auto size_of_header = static_cast<char>(first_element_with_data_);
        unsigned int size_of_data = size_of_header
            + sizeof(from)
            + sizeof(map.size());
        for (auto iterator = map.begin(); iterator != map.end(); ++iterator)
        {
            size_of_data += sizeof(iterator->first);
            size_of_data += static_cast<unsigned int>(iterator->second->get_modelName().size()) + sizeof(int);
            size_of_data += sizeof(typeRNum) * static_cast<unsigned int>(iterator->second->get_modelParameters().size()) + sizeof(int);
        }

        std::shared_ptr< std::vector<char> > data(new std::vector<char>(size_of_data, 0));

        // size of data
        DataConversion::insert_into_charArray(data, pos, size_of_data);

        // index
        DataConversion::insert_into_charArray(data, pos, index);

        //from
        DataConversion::insert_into_charArray(data, pos, from);

        // number of coupling models
        const unsigned int size_of_map = static_cast<unsigned int>(map.size());
        DataConversion::insert_into_charArray(data, pos, size_of_map);

        // maps
        for (auto iterator = map.begin(); iterator != map.end(); ++iterator)
        {
            DataConversion::insert_into_charArray(data, pos, iterator->first);
            DataConversion::insert_into_charArray(data, pos, iterator->second->get_modelName());
            DataConversion::insert_into_charArray(data, pos, iterator->second->get_modelParameters());
        }

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_get_agentModel_from_agent()
    {
		const char index = static_cast<char>(index::get_agentModelFromAgent);
		const auto size_of_header = static_cast<char>(first_element_with_data_);
        const unsigned int sizeOfData = size_of_header;
        unsigned int pos = 0;

        std::shared_ptr< std::vector<char> > data(new std::vector<char>(size_of_header, 0));

        // size of data
        DataConversion::insert_into_charArray(data, pos, sizeOfData);

        // index
        DataConversion::insert_into_charArray(data, pos, index);

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_get_couplingModel_from_agent()
    {
		const char index = static_cast<char>(index::get_couplingModelFromAgent);
		const auto size_of_header = static_cast<char>(first_element_with_data_);
        const unsigned int sizeOfData = size_of_header;
        unsigned int pos = 0;

        std::shared_ptr< std::vector<char> > data(new std::vector<char>(sizeOfData, 0));

        // size of data
        DataConversion::insert_into_charArray(data, pos, sizeOfData);

        // index
        DataConversion::insert_into_charArray(data, pos, index);

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_send_agentState_for_simulation(const AgentState& state, const int from)
    {
        auto data = buildProtocol_send_agentState(state, from);

		// change index
        (*data)[position_of_index_in_protocol_] = static_cast<char>(index::send_agentState_for_simulation);

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_get_agentState_for_simulation()
    {
        const char index = static_cast<char>(index::get_agentState_for_simulation);
		unsigned int pos = 0;
		const auto size_of_header = static_cast<char>(first_element_with_data_);
        const unsigned int size_of_data = size_of_header;

        std::shared_ptr< std::vector<char> > data(new std::vector<char>(size_of_data, 0));

        // size of data
        DataConversion::insert_into_charArray(data, pos, size_of_data);

        // index
        DataConversion::insert_into_charArray(data, pos, index);

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_acknowledge_received_optimizationInfo()
    {
        const char index = static_cast<char>(index::received_acknowledgement_received_optimizationInfo);
		unsigned int pos = 0;
		const auto size_of_header = static_cast<char>(first_element_with_data_);
        const unsigned int size_of_data = size_of_header;

        std::shared_ptr< std::vector<char> > data(new std::vector<char>(size_of_data, 0));

        // size of data
        DataConversion::insert_into_charArray(data, pos, size_of_data);

        // index
        DataConversion::insert_into_charArray(data, pos, index);

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_acknowledge_executed_ADMMstep()
    {
        const char index = static_cast<char>(index::received_acknowledgement_executed_ADMMstep);
		unsigned int pos = 0;
		const auto size_of_header = static_cast<char>(first_element_with_data_);
        const unsigned int size_of_data = size_of_header;

        std::shared_ptr< std::vector<char> > data(new std::vector<char>(size_of_data, 0));

        // size of data
        DataConversion::insert_into_charArray(data, pos, size_of_data);

        // index
        DataConversion::insert_into_charArray(data, pos, index);

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_send_communicationInfo(const CommunicationInfo& info)
    {
        const char index = static_cast<char>(index::send_communicationInfo);
		unsigned int pos = 0;
		const auto size_of_header = static_cast<char>(first_element_with_data_);
        const unsigned int sizeOfData = size_of_header
            + static_cast<unsigned int>(sizeof(info.id_)
                + info.ip_.size() + sizeof(int)
                + info.port_.size() + sizeof(int)
                + info.agent_type_.size() + sizeof(int));

        std::shared_ptr< std::vector<char> > data(new std::vector<char>(sizeOfData, 0));

        // size of data
        DataConversion::insert_into_charArray(data, pos, sizeOfData);

        // index
        DataConversion::insert_into_charArray(data, pos, index);

        // communication info
        DataConversion::insert_into_charArray(data, pos, info.id_);
        DataConversion::insert_into_charArray(data, pos, info.ip_);
        DataConversion::insert_into_charArray(data, pos, info.port_);
        DataConversion::insert_into_charArray(data, pos, info.agent_type_);

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_triggerStep(const ADMMStep& step)
	{
		const auto size_of_header = static_cast<char>(first_element_with_data_);
        const unsigned int sizeOfData = size_of_header + sizeof(int);
        std::shared_ptr< std::vector<char> > data(new std::vector<char>(sizeOfData, 0));
        unsigned int pos = 0;

        // include sizeOfData
        DataConversion::insert_into_charArray(data, pos, sizeOfData);

        (*data)[position_of_index_in_protocol_] = static_cast<char>(index::triggerStep);
        ++pos;

        // include step
        DataConversion::insert_into_charArray(data, pos, DataConversion::ADMMStep_to_int(step));

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_fromCommunication_deregistered_coupling(const CouplingInfo& info)
    {
        auto data = buildProtocol_deregister_coupling(info);

		// change index
        (*data)[position_of_index_in_protocol_] = static_cast<char>(index::deregistered_coupling);

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_deregister_coupling(const CouplingInfo& info)
    {
        const char index = static_cast<char>(index::deregister_coupling);
		unsigned int pos = 0;
		const auto size_of_header = static_cast<char>(first_element_with_data_);
        const unsigned int size_of_data = size_of_header
            + static_cast<unsigned int>(sizeof(info.agent_id_)
                + sizeof(info.neighbor_id_)
                + info.model_name_.size() + sizeof(int)
                + sizeof(typeRNum) * info.model_parameters_.size() + sizeof(int));

        std::shared_ptr< std::vector<char> > data(new std::vector<char>(size_of_data, 0));

        // size of data
        DataConversion::insert_into_charArray(data, pos, size_of_data);

        // index
        DataConversion::insert_into_charArray(data, pos, index);

        // coupling info
        DataConversion::insert_into_charArray(data, pos, info.agent_id_);
        DataConversion::insert_into_charArray(data, pos, info.neighbor_id_);
        DataConversion::insert_into_charArray(data, pos, info.model_name_);
        DataConversion::insert_into_charArray(data, pos, info.model_parameters_);

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_register_coupling(const CouplingInfo& coupling_info)
    {
        const char index = static_cast<char>(index::register_coupling);
		unsigned int pos = 0;
		const auto size_of_header = static_cast<char>(first_element_with_data_);
        const unsigned int size_of_data = size_of_header
            + static_cast<unsigned int>(sizeof(coupling_info.agent_id_)
                + sizeof(coupling_info.neighbor_id_)
                + coupling_info.model_name_.size() + sizeof(int)
                + sizeof(typeRNum) * coupling_info.model_parameters_.size() + sizeof(int));

        std::shared_ptr< std::vector<char> > data(new std::vector<char>(size_of_data, 0));

        // size of data
        DataConversion::insert_into_charArray(data, pos, size_of_data);

        // index
        DataConversion::insert_into_charArray(data, pos, index);

        // coupling info
        DataConversion::insert_into_charArray(data, pos, coupling_info.agent_id_);
        DataConversion::insert_into_charArray(data, pos, coupling_info.neighbor_id_);
        DataConversion::insert_into_charArray(data, pos, coupling_info.model_name_);
        DataConversion::insert_into_charArray(data, pos, coupling_info.model_parameters_);

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_deregister_agent(const AgentInfo& info)
    {
        // start with protocol for registering an agent
        auto data = buildProtocol_register_agent(info);

		// change index
        (*data)[position_of_index_in_protocol_] = static_cast<char>(index::deregister_agent);

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_register_agent(const AgentInfo& info)
    {
        const char index = static_cast<char>(index::register_agent);
        unsigned int pos = 0;

		const auto size_of_header = static_cast<char>(first_element_with_data_);
        const int size_of_data = size_of_header
            + static_cast<unsigned int>(sizeof(typeRNum) * info.cost_parameters_.size() + sizeof(int)
                + sizeof(info.id_)
                + info.model_name_.size() + sizeof(int)
                + sizeof(typeRNum) * info.model_parameters_.size() + sizeof(int));

        std::shared_ptr< std::vector<char> > data(new std::vector<char>(size_of_data, 0));

        // size of data
        DataConversion::insert_into_charArray(data, pos, size_of_data);

        // index
        DataConversion::insert_into_charArray(data, pos, index);

        // agent info
        DataConversion::insert_into_charArray(data, pos, info.cost_parameters_);
        DataConversion::insert_into_charArray(data, pos, info.id_);
        DataConversion::insert_into_charArray(data, pos, info.model_name_);
        DataConversion::insert_into_charArray(data, pos, info.model_parameters_);

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_successfully_registered_coupling(const CouplingInfo& coupling_info, const AgentInfo& agent_info, const CommunicationInfo& comm_info)
    {
        const char index = static_cast<char>(index::successfully_registered_coupling);
		unsigned int pos = 0;
		const auto size_of_header = static_cast<char>(first_element_with_data_);
        const unsigned int size_of_data = size_of_header
            // coupling info
            + static_cast<unsigned int>(sizeof(coupling_info.agent_id_)
                + sizeof(coupling_info.neighbor_id_)
                + coupling_info.model_name_.size() + sizeof(int)
                + sizeof(typeRNum) * coupling_info.model_parameters_.size() + sizeof(int)
                // agent info
                + sizeof(agent_info.id_)
                + agent_info.model_name_.size() + sizeof(int)
                + sizeof(typeRNum) * agent_info.model_parameters_.size() + sizeof(int)
                + sizeof(typeRNum) * agent_info.cost_parameters_.size() + sizeof(int)
                // communication info
                + sizeof(comm_info.id_)
                + comm_info.ip_.size() + sizeof(int)
                + comm_info.port_.size() + sizeof(int)
                + comm_info.agent_type_.size() + sizeof(int));

        std::shared_ptr< std::vector<char> > data(new std::vector<char>(size_of_data, 0));

        // size of data
        DataConversion::insert_into_charArray(data, pos, size_of_data);

        // index
        DataConversion::insert_into_charArray(data, pos, index);

        // coupling info
        DataConversion::insert_into_charArray(data, pos, coupling_info.agent_id_);
        DataConversion::insert_into_charArray(data, pos, coupling_info.neighbor_id_);
        DataConversion::insert_into_charArray(data, pos, coupling_info.model_name_);
        DataConversion::insert_into_charArray(data, pos, coupling_info.model_parameters_);

        // agent info
        DataConversion::insert_into_charArray(data, pos, agent_info.id_);
        DataConversion::insert_into_charArray(data, pos, agent_info.model_name_);
        DataConversion::insert_into_charArray(data, pos, agent_info.model_parameters_);
        DataConversion::insert_into_charArray(data, pos, agent_info.cost_parameters_);

        // communication info
        DataConversion::insert_into_charArray(data, pos, comm_info.id_);
        DataConversion::insert_into_charArray(data, pos, comm_info.ip_);
        DataConversion::insert_into_charArray(data, pos, comm_info.port_);
        DataConversion::insert_into_charArray(data, pos, comm_info.agent_type_);

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_send_agentState(const AgentState& state, const int from)
    {
        const char index = static_cast<char>(index::send_agentState);
		unsigned int pos = 0;
		const auto size_of_header = static_cast<char>(first_element_with_data_);
        const unsigned int size_of_data = size_of_header
            + static_cast<unsigned int>(sizeof(from)
                + sizeof(state.i_)
                + sizeof(state.t0_)
                + sizeof(typeRNum) * state.t_.size() + sizeof(int)
                + sizeof(typeRNum) * state.x_.size() + sizeof(int)
                + sizeof(typeRNum) * state.u_.size() + sizeof(int)
                + sizeof(typeRNum) * state.v_.size() + sizeof(int));

        std::shared_ptr< std::vector<char> > data(new std::vector<char>(size_of_data, 0));

        // size of data
        DataConversion::insert_into_charArray(data, pos, size_of_data);

        // index
        DataConversion::insert_into_charArray(data, pos, index);

        // from
        DataConversion::insert_into_charArray(data, pos, from);

        // agent state
        DataConversion::insert_into_charArray(data, pos, state.i_);
        DataConversion::insert_into_charArray(data, pos, state.t0_);
        DataConversion::insert_into_charArray(data, pos, state.t_);
        DataConversion::insert_into_charArray(data, pos, state.x_);
        DataConversion::insert_into_charArray(data, pos, state.u_);
        DataConversion::insert_into_charArray(data, pos, state.v_);

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_send_desiredAgentState(const AgentState& state, const int from)
    {
        auto data = buildProtocol_send_agentState(state, from);

		// change index
        (*data)[position_of_index_in_protocol_] = static_cast<char>(index::send_desiredAgentState);

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_send_couplingState(const CouplingState& state, const int from)
    {
        const char index = static_cast<char>(index::send_couplingState);
        unsigned int pos = 0;
        const auto size_of_header = static_cast<char>(first_element_with_data_);
        const unsigned int size_of_data = size_of_header
            + static_cast<unsigned int>(sizeof(from)
                + sizeof(state.i_)
                + sizeof(state.t0_)
                + sizeof(typeRNum) * state.t_.size() + sizeof(int)
                + sizeof(typeRNum) * state.z_x_.size() + sizeof(int)
                + sizeof(typeRNum) * state.z_u_.size() + sizeof(int)
                + sizeof(typeRNum) * state.z_v_.size() + sizeof(int));

        std::shared_ptr< std::vector<char> > data(new std::vector<char>(size_of_data, 0));

        // size of data
        DataConversion::insert_into_charArray(data, pos, size_of_data);

        // index
        DataConversion::insert_into_charArray(data, pos, index);

        // from
        DataConversion::insert_into_charArray(data, pos, from);

        // coupling state
        DataConversion::insert_into_charArray(data, pos, state.i_);
        DataConversion::insert_into_charArray(data, pos, state.t0_);
        DataConversion::insert_into_charArray(data, pos, state.t_);
        DataConversion::insert_into_charArray(data, pos, state.z_x_);
        DataConversion::insert_into_charArray(data, pos, state.z_u_);
        DataConversion::insert_into_charArray(data, pos, state.z_v_);

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_send_couplingState(
        const CouplingState& state1,
        const CouplingState& state2,
        const int from)
    {
        const char index = static_cast<char>(index::send_two_couplingStates);
		unsigned int pos = 0;
		const auto size_of_header = static_cast<char>(first_element_with_data_);
        const unsigned int size_of_data = size_of_header
            + static_cast<unsigned int>(sizeof(from)
                // state1
                + sizeof(state1.i_)
                + sizeof(state1.t0_)
                + sizeof(typeRNum) * state1.t_.size() + sizeof(int)
                + sizeof(typeRNum) * state1.z_x_.size() + sizeof(int)
                + sizeof(typeRNum) * state1.z_u_.size() + sizeof(int)
                + sizeof(typeRNum) * state1.z_v_.size() + sizeof(int)
                // state2
                + sizeof(state2.i_)
                + sizeof(state2.t0_)
                + sizeof(typeRNum) * state2.t_.size() + sizeof(int)
                + sizeof(typeRNum) * state2.z_x_.size() + sizeof(int)
                + sizeof(typeRNum) * state2.z_u_.size() + sizeof(int)
                + sizeof(typeRNum) * state2.z_v_.size() + sizeof(int));

        std::shared_ptr< std::vector<char> > data(new std::vector<char>(size_of_data, 0));

        // size of data
        DataConversion::insert_into_charArray(data, pos, size_of_data);

        // index
        DataConversion::insert_into_charArray(data, pos, index);

        // from
        DataConversion::insert_into_charArray(data, pos, from);

        // coupling state 1
        DataConversion::insert_into_charArray(data, pos, state1.i_);
        DataConversion::insert_into_charArray(data, pos, state1.t0_);
        DataConversion::insert_into_charArray(data, pos, state1.t_);
        DataConversion::insert_into_charArray(data, pos, state1.z_x_);
        DataConversion::insert_into_charArray(data, pos, state1.z_u_);
        DataConversion::insert_into_charArray(data, pos, state1.z_v_);

        // coupling state 2
        DataConversion::insert_into_charArray(data, pos, state2.i_);
        DataConversion::insert_into_charArray(data, pos, state2.t0_);
        DataConversion::insert_into_charArray(data, pos, state2.t_);
        DataConversion::insert_into_charArray(data, pos, state2.z_x_);
        DataConversion::insert_into_charArray(data, pos, state2.z_u_);
        DataConversion::insert_into_charArray(data, pos, state2.z_v_);

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_send_multiplierPenaltyState(const MultiplierState& multiplier, const PenaltyState& penalty, const int from)
    {
        const char index = static_cast<char>(index::send_multiplierPenaltyState);
		unsigned int pos = 0;
		const auto size_of_header = static_cast<char>(first_element_with_data_);
        const unsigned int size_of_data = size_of_header
            // from
            + static_cast<unsigned int>(static_cast<unsigned int>(sizeof(from)
                // multiplier state
                + sizeof(multiplier.i_)
                + sizeof(multiplier.t0_)
                + sizeof(typeRNum) * multiplier.t_.size() + sizeof(int)
                + sizeof(typeRNum) * multiplier.mu_x_.size() + sizeof(int)
                + sizeof(typeRNum) * multiplier.mu_u_.size() + sizeof(int)
                + sizeof(typeRNum) * multiplier.mu_v_.size() + sizeof(int)
                // penalty state
                + sizeof(penalty.i_)
                + sizeof(penalty.t0_)
                + sizeof(typeRNum) * penalty.t_.size() + sizeof(int)
                + sizeof(typeRNum) * penalty.rho_x_.size() + sizeof(int)
                + sizeof(typeRNum) * penalty.rho_u_.size() + sizeof(int)
                + sizeof(typeRNum) * penalty.rho_v_.size() + sizeof(int)));

        std::shared_ptr< std::vector<char> > data(new std::vector<char>(size_of_data, 0));

        // size of data
        DataConversion::insert_into_charArray(data, pos, size_of_data);

        // index
        DataConversion::insert_into_charArray(data, pos, index);

        // from
        DataConversion::insert_into_charArray(data, pos, from);

        // multiplier state
        DataConversion::insert_into_charArray(data, pos, multiplier.i_);
        DataConversion::insert_into_charArray(data, pos, multiplier.t0_);
        DataConversion::insert_into_charArray(data, pos, multiplier.t_);
        DataConversion::insert_into_charArray(data, pos, multiplier.mu_x_);
        DataConversion::insert_into_charArray(data, pos, multiplier.mu_u_);
        DataConversion::insert_into_charArray(data, pos, multiplier.mu_v_);

        // penalty state
        DataConversion::insert_into_charArray(data, pos, penalty.i_);
        DataConversion::insert_into_charArray(data, pos, penalty.t0_);
        DataConversion::insert_into_charArray(data, pos, penalty.t_);
        DataConversion::insert_into_charArray(data, pos, penalty.rho_x_);
        DataConversion::insert_into_charArray(data, pos, penalty.rho_u_);
        DataConversion::insert_into_charArray(data, pos, penalty.rho_v_);

        return data;
    }

    const std::shared_ptr<std::vector<char>> ProtocolCommunication::buildProtocol_send_convergenceFlag(const bool converged, const int from)
	{
		const auto size_of_header = static_cast<char>(first_element_with_data_);
        const unsigned int sizeOfData = size_of_header + 1 + sizeof(int);
        const char index = static_cast<char>(index::send_convergenceFlag);

        //safe space for size of data
        std::shared_ptr< std::vector<char> > data(new std::vector<char>(sizeOfData, 0));
        unsigned int pos = 0;

        // include sizeOfData
        DataConversion::insert_into_charArray(data, pos, sizeOfData);

        // index 
        DataConversion::insert_into_charArray(data, pos, index);

        // include from
        DataConversion::insert_into_charArray(data, pos, from);

        // include convergence flag
        DataConversion::insert_into_charArray(data, pos, converged);

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_send_numberOfNeighbors(const int number, const int from)
	{
		const auto size_of_header = static_cast<char>(first_element_with_data_);
        const unsigned int sizeOfData = size_of_header + 2 * sizeof(int);
        const char index = static_cast<char>(index::send_numberOfNeighbors);

        std::shared_ptr< std::vector<char> > data(new std::vector<char>(sizeOfData, 0));
        unsigned int pos = 0;

        // include sizeOfdata
        DataConversion::insert_into_charArray(data, pos, sizeOfData);

        // index 
        DataConversion::insert_into_charArray(data, pos, index);

        // include from
        DataConversion::insert_into_charArray(data, pos, from);

        // send number of neighbor
        DataConversion::insert_into_charArray(data, pos, number);

        return data;
    }

    const std::shared_ptr< std::vector<char> > ProtocolCommunication::buildProtocol_send_optimizationInfo(const OptimizationInfo& info)
	{
		const auto size_of_header = static_cast<char>(first_element_with_data_);
        const unsigned int size_of_data = size_of_header // header
            // COMMON
            + static_cast<unsigned int>(sizeof(info.COMMON_Thor_)
                + sizeof(info.COMMON_dt_)
                + sizeof(info.COMMON_Nhor_)
                + sizeof(info.COMMON_ShiftControl_)
                + info.COMMON_Integrator_.size() + sizeof(int)
                // GRAMPC
                + sizeof(info.GRAMPC_MaxGradIter_)
                + sizeof(info.GRAMPC_MaxMultIter_)
                + sizeof(info.GRAMPC_PenaltyMin_)
                + sizeof(info.GRAMPC_PenaltyMax_)
                + sizeof(info.GRAMPC_AugLagUpdateGradientRelTol_)
                + info.GRAMPC_Integrator_.size() + sizeof(int)
                + info.GRAMPC_LineSearchType_.size() + sizeof(int)
                + sizeof(info.GRAMPC_PenaltyIncreaseFactor_)
                + sizeof(info.GRAMPC_PenaltyDecreaseFactor_)
                + sizeof(info.GRAMPC_LineSearchMax_)
                + sizeof(info.GRAMPC_LineSearchMin_)
                + info.GRAMPC_ConvergenceCheck_.size() + sizeof(int)
                + sizeof(info.GRAMPC_ConvergenceGradientRelTol_)
                + sizeof(typeRNum) * info.GRAMPC_ConstraintsAbsTol_.size() + sizeof(int)
                + sizeof(info.GRAMPC_PenaltyIncreaseThreshold_)
                + sizeof(info.GRAMPC_LineSearchInit_)
                //ADMM
                + sizeof(info.ADMM_maxIterations_)
                + sizeof(info.ADMM_innerIterations_)
                + sizeof(info.ADMM_ConvergenceTolerance_)
                + sizeof(info.ADMM_PenaltyIncreaseFactor_)
                + sizeof(info.ADMM_PenaltyDecreaseFactor_)
                + sizeof(info.ADMM_PenaltyMin_)
                + sizeof(info.ADMM_PenaltyMax_)
                + sizeof(info.ADMM_PenaltyInit_)
                + sizeof(info.ADMM_AdaptPenaltyParameter_)
                // APPROX
                + sizeof(info.APPROX_ApproximateCost_)
                + sizeof(info.APPROX_ApproximateConstraints_)
                + sizeof(info.APPROX_ApproximateDynamics_));

        std::shared_ptr< std::vector<char> >data(new std::vector<char>(size_of_data, 0));
        unsigned int pos = 0;
        const char index = static_cast<char>(index::send_optimizationInfo);

        // size of data
        DataConversion::insert_into_charArray(data, pos, size_of_data);

        // index 
        DataConversion::insert_into_charArray(data, pos, index);

        /**************************
         * COMMON
        **************************/

        DataConversion::insert_into_charArray(data, pos, info.COMMON_Thor_);
        DataConversion::insert_into_charArray(data, pos, info.COMMON_dt_);
        DataConversion::insert_into_charArray(data, pos, info.COMMON_Nhor_);
        DataConversion::insert_into_charArray(data, pos, info.COMMON_ShiftControl_);
        DataConversion::insert_into_charArray(data, pos, info.COMMON_Integrator_);

        /**************************
         * GRAMPC
        **************************/

        DataConversion::insert_into_charArray(data, pos, info.GRAMPC_MaxGradIter_);
        DataConversion::insert_into_charArray(data, pos, info.GRAMPC_MaxMultIter_);
        DataConversion::insert_into_charArray(data, pos, info.GRAMPC_PenaltyMin_);
        DataConversion::insert_into_charArray(data, pos, info.GRAMPC_PenaltyMax_);
        DataConversion::insert_into_charArray(data, pos, info.GRAMPC_AugLagUpdateGradientRelTol_);
        DataConversion::insert_into_charArray(data, pos, info.GRAMPC_Integrator_);
        DataConversion::insert_into_charArray(data, pos, info.GRAMPC_LineSearchType_);
        DataConversion::insert_into_charArray(data, pos, info.GRAMPC_PenaltyIncreaseFactor_);
        DataConversion::insert_into_charArray(data, pos, info.GRAMPC_PenaltyDecreaseFactor_);
        DataConversion::insert_into_charArray(data, pos, info.GRAMPC_LineSearchMax_);
        DataConversion::insert_into_charArray(data, pos, info.GRAMPC_LineSearchMin_);
        DataConversion::insert_into_charArray(data, pos, info.GRAMPC_ConvergenceCheck_);
        DataConversion::insert_into_charArray(data, pos, info.GRAMPC_ConvergenceGradientRelTol_);
        DataConversion::insert_into_charArray(data, pos, info.GRAMPC_ConstraintsAbsTol_);
        DataConversion::insert_into_charArray(data, pos, info.GRAMPC_PenaltyIncreaseThreshold_);
        DataConversion::insert_into_charArray(data, pos, info.GRAMPC_LineSearchInit_);

        /**************************
         * ADMM
        **************************/

        DataConversion::insert_into_charArray(data, pos, info.ADMM_maxIterations_);
        DataConversion::insert_into_charArray(data, pos, info.ADMM_innerIterations_);
        DataConversion::insert_into_charArray(data, pos, info.ADMM_ConvergenceTolerance_);
        DataConversion::insert_into_charArray(data, pos, info.ADMM_PenaltyIncreaseFactor_);
        DataConversion::insert_into_charArray(data, pos, info.ADMM_PenaltyDecreaseFactor_);
        DataConversion::insert_into_charArray(data, pos, info.ADMM_PenaltyMin_);
        DataConversion::insert_into_charArray(data, pos, info.ADMM_PenaltyMax_);
        DataConversion::insert_into_charArray(data, pos, info.ADMM_PenaltyInit_);
        DataConversion::insert_into_charArray(data, pos, info.ADMM_AdaptPenaltyParameter_);

        /**************************
         * APPROX
        **************************/

        DataConversion::insert_into_charArray(data, pos, info.APPROX_ApproximateCost_);
        DataConversion::insert_into_charArray(data, pos, info.APPROX_ApproximateConstraints_);
        DataConversion::insert_into_charArray(data, pos, info.APPROX_ApproximateDynamics_);

        return data;
    }

    /****************************************************
     * functions to rebuild protocol
	****************************************************/

    const SolutionPtr ProtocolCommunication::buildFromProtocol_solution(const std::vector<char>& data)
	{
		unsigned int pos = static_cast<char>(first_element_with_data_);

		SolutionPtr solution(new Solution);

		DataConversion::read_from_charArray(data, pos, solution->agentState_.i_);
		DataConversion::read_from_charArray(data, pos, solution->agentState_.t0_);
		DataConversion::read_from_charArray(data, pos, solution->agentState_.t_);
		DataConversion::read_from_charArray(data, pos, solution->agentState_.u_);
		DataConversion::read_from_charArray(data, pos, solution->agentState_.v_);
		DataConversion::read_from_charArray(data, pos, solution->agentState_.x_);

		DataConversion::read_from_charArray(data, pos, solution->predicted_agentState_.i_);
		DataConversion::read_from_charArray(data, pos, solution->predicted_agentState_.t0_);
		DataConversion::read_from_charArray(data, pos, solution->predicted_agentState_.t_);
		DataConversion::read_from_charArray(data, pos, solution->predicted_agentState_.u_);
		DataConversion::read_from_charArray(data, pos, solution->predicted_agentState_.v_);
		DataConversion::read_from_charArray(data, pos, solution->predicted_agentState_.x_);

		DataConversion::read_from_charArray(data, pos, solution->cost_);

		return solution;
	}

    const int ProtocolCommunication::buildFromProtocol_numberOfActiveCouplings(const std::vector<char>& data)
    {
        unsigned int pos = static_cast<char>(first_element_with_data_);
        int number = 0;
        DataConversion::read_from_charArray(data, pos, number);
        return number;
    }

    const typeRNum ProtocolCommunication::buildFromProtocol_t0_from_simulatedState(const std::vector<char>& data)
    {
        unsigned int pos = static_cast<char>(first_element_with_data_);

        // skip xnext
        DataConversion::skip_in_charArray_vecTypeRNum(data, pos);

        // skip dt
        DataConversion::skip_in_charArray_typeRNum(data, pos);

        // read t0
        typeRNum t0 = 0;
        DataConversion::read_from_charArray(data, pos, t0);

        return t0;
    }

    const typeRNum ProtocolCommunication::buildFromProtocol_dt_from_simulatedState(const std::vector<char>& data)
    {
        unsigned int pos = static_cast<char>(first_element_with_data_);

        // skip xnext
        DataConversion::skip_in_charArray_vecTypeRNum(data, pos);

        // read dt
        typeRNum dt = 0;
        DataConversion::read_from_charArray(data, pos, dt);
        return dt;
    }

    const std::shared_ptr< std::vector<typeRNum> > ProtocolCommunication::buildFromProtocol_xnext_from_simulatedState(const std::vector<char>& data)
    {
        unsigned int pos = static_cast<char>(first_element_with_data_);

        std::vector<typeRNum> vec(0, 0);
        DataConversion::read_from_charArray(data, pos, vec);

        // read xnext
        return std::make_shared < std::vector<typeRNum>>(vec);
    }

    const AgentModelPtr ProtocolCommunication::buildFromProtocol_agentModel(const std::vector<char>& data, const LoggingPtr& log)
    {
        AgentInfo info;
        unsigned int pos = static_cast<char>(first_element_with_data_);
        // skip from
        DataConversion::skip_in_charArray_int(data, pos);

        DataConversion::read_from_charArray(data, pos, info.model_name_);
        DataConversion::read_from_charArray(data, pos, info.model_parameters_);
        DataConversion::read_from_charArray(data, pos, info.cost_parameters_);

        ModelFactoryPtr factory(new GeneralModelFactory(log));

        return factory->create_agentModel(info);
    }

    const std::shared_ptr< std::map<int, CouplingModelPtr> > ProtocolCommunication::buildFromProtocol_couplingModel(const std::vector<char>& data, const LoggingPtr& log)
    {
        std::shared_ptr< std::map<int, CouplingModelPtr> > map(new std::map<int, CouplingModelPtr>);
        CouplingInfo info;
        ModelFactoryPtr factory(new GeneralModelFactory(log));

        unsigned int pos = static_cast<char>(first_element_with_data_);
        // skip from
        DataConversion::skip_in_charArray_int(data, pos);

        // read number of entries
        int entries = 0;
        DataConversion::read_from_charArray(data, pos, entries);
        int id = 0;

        for (int i = 0; i < entries; ++i)
        {
            DataConversion::read_from_charArray(data, pos, id);
            DataConversion::read_from_charArray(data, pos, info.model_name_);
            DataConversion::read_from_charArray(data, pos, info.model_parameters_);

            CouplingModelPtr coupling_model = factory->create_couplingModel(info);

            // insert map entry
            map->insert(std::make_pair(id, coupling_model));
        }

        return map;
    }

    const int ProtocolCommunication::buildFromProtocol_numberOfNeighbors(const std::vector<char>& data)
    {
        unsigned int pos = static_cast<char>(first_element_with_data_);
        // skip from
        DataConversion::skip_in_charArray_int(data, pos);

        int number = 0;
        DataConversion::read_from_charArray(data, pos, number);

        return number;
    }

    const bool ProtocolCommunication::buildFromProtocol_convergenceFlag(const std::vector<char>& data)
    {
        unsigned int pos = static_cast<char>(first_element_with_data_);
        // skip from
        DataConversion::skip_in_charArray_int(data, pos);

        bool flag = false;
        DataConversion::read_from_charArray(data, pos, flag);

        return flag;
    }
    const int ProtocolCommunication::buildFromProtocol_from(const std::vector<char>& data)
    {
        unsigned int pos = static_cast<char>(first_element_with_data_);
        int from = 0;
        DataConversion::read_from_charArray(data, pos, from);
        return from;
    }

    const std::shared_ptr< CommunicationInfo > ProtocolCommunication::buildFromProtocol_communicationInfo(const std::vector<char>& data)
    {
        std::shared_ptr< CommunicationInfo > info(new CommunicationInfo);
        unsigned int pos = static_cast<char>(first_element_with_data_);

        DataConversion::read_from_charArray(data, pos, info->id_);
        DataConversion::read_from_charArray(data, pos, info->ip_);
        DataConversion::read_from_charArray(data, pos, info->port_);
        DataConversion::read_from_charArray(data, pos, info->agent_type_);

        return info;
    }

    const std::shared_ptr< std::vector<CommunicationInfo> > ProtocolCommunication::buildFromProtocol_communicationInfoArray(const std::vector<char>& data)
    {
        unsigned int pos = static_cast<char>(first_element_with_data_);
        unsigned int sizeOfData = 0;
        DataConversion::read_from_charArray(data, pos, sizeOfData);

        std::shared_ptr< std::vector<CommunicationInfo> > infos(new std::vector<CommunicationInfo>);

        while (pos < sizeOfData)
        {
            infos->push_back(CommunicationInfo());

            DataConversion::read_from_charArray(data, pos, infos->back().id_);
            DataConversion::read_from_charArray(data, pos, infos->back().ip_);
            DataConversion::read_from_charArray(data, pos, infos->back().port_);
            DataConversion::read_from_charArray(data, pos, infos->back().agent_type_);
        }

        return infos;
    }

    const ADMMStep ProtocolCommunication::buildFromProtocol_ADMMStep(const std::vector<char>& data)
    {
        unsigned int pos = static_cast<char>(first_element_with_data_);
        int step = 0;
        DataConversion::read_from_charArray(data, pos, step);

        return DataConversion::Int_to_ADMMStep(step);
    }

    const std::shared_ptr< AgentInfo > ProtocolCommunication::buildFromProtocol_agentInfo(const std::vector<char>& data)
    {
        std::shared_ptr< AgentInfo > info(new AgentInfo);
        unsigned int pos = static_cast<char>(first_element_with_data_);

        // read cost parameters
        DataConversion::read_from_charArray(data, pos, info->cost_parameters_);
        DataConversion::read_from_charArray(data, pos, info->id_);
        DataConversion::read_from_charArray(data, pos, info->model_name_);
        DataConversion::read_from_charArray(data, pos, info->model_parameters_);

        return info;
    }

    const std::shared_ptr< CouplingInfo > ProtocolCommunication::buildFromProtocol_couplingInfo(const std::vector<char>& data)
    {
        std::shared_ptr< CouplingInfo > info(new CouplingInfo);
        unsigned int pos = static_cast<char>(first_element_with_data_);

        DataConversion::read_from_charArray(data, pos, info->agent_id_);
        DataConversion::read_from_charArray(data, pos, info->neighbor_id_);
        DataConversion::read_from_charArray(data, pos, info->model_name_);
        DataConversion::read_from_charArray(data, pos, info->model_parameters_);

        return info;
    }

    const std::shared_ptr< CommunicationInfo > ProtocolCommunication::buildFromProtocol_communicationInfo_from_successfully_registered_coupling(const std::vector<char>& data)
    {
        std::shared_ptr< CommunicationInfo > info(new CommunicationInfo);
        unsigned int pos = static_cast<char>(first_element_with_data_);

        // skip coupling info
        DataConversion::skip_in_charArray_int(data, pos);
        DataConversion::skip_in_charArray_int(data, pos);
        DataConversion::skip_in_charArray_string(data, pos);
        DataConversion::skip_in_charArray_vecTypeRNum(data, pos);

        // skip agent info
        DataConversion::skip_in_charArray_int(data, pos);
        DataConversion::skip_in_charArray_string(data, pos);
        DataConversion::skip_in_charArray_vecTypeRNum(data, pos);
        DataConversion::skip_in_charArray_vecTypeRNum(data, pos);

        // read communication info
        DataConversion::read_from_charArray(data, pos, info->id_);
        DataConversion::read_from_charArray(data, pos, info->ip_);
        DataConversion::read_from_charArray(data, pos, info->port_);
        DataConversion::read_from_charArray(data, pos, info->agent_type_);

        return info;
    }

    const std::shared_ptr< AgentInfo > ProtocolCommunication::buildFromProtocol_agentInfo_from_successfully_registered_coupling(const std::vector<char>& data)
    {
        std::shared_ptr< AgentInfo > info(new AgentInfo);
        unsigned int pos = static_cast<char>(first_element_with_data_);

        // skip coupling info
        DataConversion::skip_in_charArray_int(data, pos);
        DataConversion::skip_in_charArray_int(data, pos);
        DataConversion::skip_in_charArray_string(data, pos);
        DataConversion::skip_in_charArray_vecTypeRNum(data, pos);

        // read agent info
        DataConversion::read_from_charArray(data, pos, info->id_);
        DataConversion::read_from_charArray(data, pos, info->model_name_);
        DataConversion::read_from_charArray(data, pos, info->model_parameters_);
        DataConversion::read_from_charArray(data, pos, info->cost_parameters_);

        return info;
    }

    const std::shared_ptr< CouplingInfo > ProtocolCommunication::buildFromProtocol_couplingInfo_from_successfully_registered_coupling(const std::vector<char>& data)
    {
        std::shared_ptr< CouplingInfo > info(new CouplingInfo);
        unsigned int pos = static_cast<char>(first_element_with_data_);

        DataConversion::read_from_charArray(data, pos, info->agent_id_);
        DataConversion::read_from_charArray(data, pos, info->neighbor_id_);
        DataConversion::read_from_charArray(data, pos, info->model_name_);
        DataConversion::read_from_charArray(data, pos, info->model_parameters_);

        return info;

    }

    const std::shared_ptr< OptimizationInfo > ProtocolCommunication::buildFromProtocol_optimizationInfo(const std::vector<char>& data)
    {
        std::shared_ptr< OptimizationInfo > info(new OptimizationInfo);
        unsigned int pos = static_cast<char>(first_element_with_data_);

        /**************************
         * COMMON
        **************************/
        DataConversion::read_from_charArray(data, pos, info->COMMON_Thor_);
        DataConversion::read_from_charArray(data, pos, info->COMMON_dt_);
        DataConversion::read_from_charArray(data, pos, info->COMMON_Nhor_);
        DataConversion::read_from_charArray(data, pos, info->COMMON_ShiftControl_);
        DataConversion::read_from_charArray(data, pos, info->COMMON_Integrator_);

        /**************************
         * GRAMPC
        **************************/
        DataConversion::read_from_charArray(data, pos, info->GRAMPC_MaxGradIter_);
        DataConversion::read_from_charArray(data, pos, info->GRAMPC_MaxMultIter_);
        DataConversion::read_from_charArray(data, pos, info->GRAMPC_PenaltyMin_);
        DataConversion::read_from_charArray(data, pos, info->GRAMPC_PenaltyMax_);
        DataConversion::read_from_charArray(data, pos, info->GRAMPC_AugLagUpdateGradientRelTol_);
        DataConversion::read_from_charArray(data, pos, info->GRAMPC_Integrator_);
        DataConversion::read_from_charArray(data, pos, info->GRAMPC_LineSearchType_);
        DataConversion::read_from_charArray(data, pos, info->GRAMPC_PenaltyIncreaseFactor_);
        DataConversion::read_from_charArray(data, pos, info->GRAMPC_PenaltyDecreaseFactor_);
        DataConversion::read_from_charArray(data, pos, info->GRAMPC_LineSearchMax_);
        DataConversion::read_from_charArray(data, pos, info->GRAMPC_LineSearchMin_);
        DataConversion::read_from_charArray(data, pos, info->GRAMPC_ConvergenceCheck_);
        DataConversion::read_from_charArray(data, pos, info->GRAMPC_ConvergenceGradientRelTol_);
        DataConversion::read_from_charArray(data, pos, info->GRAMPC_ConstraintsAbsTol_);
        DataConversion::read_from_charArray(data, pos, info->GRAMPC_PenaltyIncreaseThreshold_);
        DataConversion::read_from_charArray(data, pos, info->GRAMPC_LineSearchInit_);

        /**************************
         * ADMM
        **************************/
        DataConversion::read_from_charArray(data, pos, info->ADMM_maxIterations_);
        DataConversion::read_from_charArray(data, pos, info->ADMM_innerIterations_);
        DataConversion::read_from_charArray(data, pos, info->ADMM_ConvergenceTolerance_);
        DataConversion::read_from_charArray(data, pos, info->ADMM_PenaltyIncreaseFactor_);
        DataConversion::read_from_charArray(data, pos, info->ADMM_PenaltyDecreaseFactor_);
        DataConversion::read_from_charArray(data, pos, info->ADMM_PenaltyMin_);
        DataConversion::read_from_charArray(data, pos, info->ADMM_PenaltyMax_);
        DataConversion::read_from_charArray(data, pos, info->ADMM_PenaltyInit_);
        DataConversion::read_from_charArray(data, pos, info->ADMM_AdaptPenaltyParameter_);

        /**************************
         * APPROX
        **************************/

        DataConversion::read_from_charArray(data, pos, info->APPROX_ApproximateCost_);
        DataConversion::read_from_charArray(data, pos, info->APPROX_ApproximateConstraints_);
        DataConversion::read_from_charArray(data, pos, info->APPROX_ApproximateDynamics_);

        return info;
    }

    const AgentStatePtr ProtocolCommunication::buildFromProtocol_agentState(const std::vector<char>& data)
    {
        AgentStatePtr state(new AgentState);
        unsigned int pos = static_cast<char>(first_element_with_data_) + sizeof(int);

        DataConversion::read_from_charArray(data, pos, state->i_);
        DataConversion::read_from_charArray(data, pos, state->t0_);
        DataConversion::read_from_charArray(data, pos, state->t_);
        DataConversion::read_from_charArray(data, pos, state->x_);
        DataConversion::read_from_charArray(data, pos, state->u_);
        DataConversion::read_from_charArray(data, pos, state->v_);

        return state;
    }

    const CouplingStatePtr ProtocolCommunication::buildFromProtocol_couplingState(const std::vector<char>& data)
    {
        CouplingStatePtr state(new CouplingState);
        unsigned int pos = static_cast<char>(first_element_with_data_);

        // skip from
        DataConversion::skip_in_charArray_int(data, pos);

        DataConversion::read_from_charArray(data, pos, state->i_);
        DataConversion::read_from_charArray(data, pos, state->t0_);
        DataConversion::read_from_charArray(data, pos, state->t_);
        DataConversion::read_from_charArray(data, pos, state->z_x_);
        DataConversion::read_from_charArray(data, pos, state->z_u_);
        DataConversion::read_from_charArray(data, pos, state->z_v_);

        return state;
    }

    const CouplingStatePtr ProtocolCommunication::buildFromProtocol_couplingState_one(const std::vector<char>& data)
    {
        return buildFromProtocol_couplingState(data);
    }

    const CouplingStatePtr ProtocolCommunication::buildFromProtocol_couplingState_two(const std::vector<char>& data)
    {
        CouplingStatePtr state(new CouplingState);
        unsigned int pos = static_cast<char>(first_element_with_data_);

        // skip from
        DataConversion::skip_in_charArray_int(data, pos);

        // skip coupling state one
        DataConversion::skip_in_charArray_int(data, pos);
        DataConversion::skip_in_charArray_typeRNum(data, pos);
        DataConversion::skip_in_charArray_vecTypeRNum(data, pos);
        DataConversion::skip_in_charArray_vecTypeRNum(data, pos);
        DataConversion::skip_in_charArray_vecTypeRNum(data, pos);
        DataConversion::skip_in_charArray_vecTypeRNum(data, pos);

        // coupling state two
        DataConversion::read_from_charArray(data, pos, state->i_);
        DataConversion::read_from_charArray(data, pos, state->t0_);
        DataConversion::read_from_charArray(data, pos, state->t_);
        DataConversion::read_from_charArray(data, pos, state->z_x_);
        DataConversion::read_from_charArray(data, pos, state->z_u_);
        DataConversion::read_from_charArray(data, pos, state->z_v_);

        return state;
    }

    const MultiplierStatePtr ProtocolCommunication::buildFromProtocol_multiplierState(const std::vector<char>& data)
    {
        MultiplierStatePtr state(new MultiplierState);
        unsigned int pos = static_cast<char>(first_element_with_data_);

        // skip from
        DataConversion::skip_in_charArray_int(data, pos);

        DataConversion::read_from_charArray(data, pos, state->i_);
        DataConversion::read_from_charArray(data, pos, state->t0_);
        DataConversion::read_from_charArray(data, pos, state->t_);
        DataConversion::read_from_charArray(data, pos, state->mu_x_);
        DataConversion::read_from_charArray(data, pos, state->mu_u_);
        DataConversion::read_from_charArray(data, pos, state->mu_v_);

        return state;
    }

    const PenaltyStatePtr ProtocolCommunication::buildFromProtocol_penaltyState(const std::vector<char>& data)
    {
        PenaltyStatePtr state(new PenaltyState);
        unsigned int pos = static_cast<char>(first_element_with_data_);

        // skip from
        DataConversion::skip_in_charArray_int(data, pos);

        // skip multiplier state
        DataConversion::skip_in_charArray_int(data, pos);
        DataConversion::skip_in_charArray_typeRNum(data, pos);
        DataConversion::skip_in_charArray_vecTypeRNum(data, pos);
        DataConversion::skip_in_charArray_vecTypeRNum(data, pos);
        DataConversion::skip_in_charArray_vecTypeRNum(data, pos);
        DataConversion::skip_in_charArray_vecTypeRNum(data, pos);

        // read penalty state
        DataConversion::read_from_charArray(data, pos, state->i_);
        DataConversion::read_from_charArray(data, pos, state->t0_);
        DataConversion::read_from_charArray(data, pos, state->t_);
        DataConversion::read_from_charArray(data, pos, state->rho_x_);
        DataConversion::read_from_charArray(data, pos, state->rho_u_);
        DataConversion::read_from_charArray(data, pos, state->rho_v_);

        return state;
    }

}