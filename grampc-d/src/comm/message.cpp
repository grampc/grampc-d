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

#include <cereal/types/map.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

CEREAL_REGISTER_TYPE(grampcd::Message_send_communication_info);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_send_communication_info)

CEREAL_REGISTER_TYPE(grampcd::Message_fromCommunication_deregistered_coupling);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_fromCommunication_deregistered_coupling)

CEREAL_REGISTER_TYPE(grampcd::Message_get_ping);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_get_ping)

CEREAL_REGISTER_TYPE(grampcd::Message_register_agent);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_register_agent)

CEREAL_REGISTER_TYPE(grampcd::Message_deregister_agent);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_deregister_agent)

CEREAL_REGISTER_TYPE(grampcd::Message_register_coupling);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_register_coupling)

CEREAL_REGISTER_TYPE(grampcd::Message_deregister_coupling);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_deregister_coupling)

CEREAL_REGISTER_TYPE(grampcd::Message_send_numberOfNeighbors);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_send_numberOfNeighbors)

CEREAL_REGISTER_TYPE(grampcd::Message_send_local_copies);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_send_local_copies)

CEREAL_REGISTER_TYPE(grampcd::Message_send_desired_agent_state);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_send_desired_agent_state)

CEREAL_REGISTER_TYPE(grampcd::Message_get_agentState_for_simulation);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_get_agentState_for_simulation)

CEREAL_REGISTER_TYPE(grampcd::Message_get_desired_agent_state_from_agent);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_get_desired_agent_state_from_agent)

CEREAL_REGISTER_TYPE(grampcd::Message_get_agent_model_from_agent);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_get_agent_model_from_agent)

CEREAL_REGISTER_TYPE(grampcd::Message_get_coupling_model_from_agent);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_get_coupling_model_from_agent)

CEREAL_REGISTER_TYPE(grampcd::Message_send_coupling_state);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_send_coupling_state)

CEREAL_REGISTER_TYPE(grampcd::Message_send_multiplier_state);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_send_multiplier_state)

CEREAL_REGISTER_TYPE(grampcd::Message_send_agent_state);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_send_agent_state)

CEREAL_REGISTER_TYPE(grampcd::Message_send_optimizationInfo);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_send_optimizationInfo)

CEREAL_REGISTER_TYPE(grampcd::Message_get_number_of_active_couplings);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_get_number_of_active_couplings)

CEREAL_REGISTER_TYPE(grampcd::Message_trigger_step);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_trigger_step)

CEREAL_REGISTER_TYPE(grampcd::Message_send_simulated_state);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_send_simulated_state)

CEREAL_REGISTER_TYPE(grampcd::Message_get_solution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_get_solution)

CEREAL_REGISTER_TYPE(grampcd::Message_send_flag_to_agents);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_send_flag_to_agents)

CEREAL_REGISTER_TYPE(grampcd::Message_acknowledge_executed_ADMM_step);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_acknowledge_executed_ADMM_step)

CEREAL_REGISTER_TYPE(grampcd::Message_send_coupling_models);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_send_coupling_models)

CEREAL_REGISTER_TYPE(grampcd::Message_send_agent_model);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_send_agent_model)

CEREAL_REGISTER_TYPE(grampcd::Message_send_agent_state_for_simulation);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_send_agent_state_for_simulation)

CEREAL_REGISTER_TYPE(grampcd::Message_successfully_deregistered_agent);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_successfully_deregistered_agent)

CEREAL_REGISTER_TYPE(grampcd::Message_send_solution);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_send_solution)

CEREAL_REGISTER_TYPE(grampcd::Message_successfully_registered_agent);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_successfully_registered_agent)

CEREAL_REGISTER_TYPE(grampcd::Message_successfully_registered_coupling);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_successfully_registered_coupling)

CEREAL_REGISTER_TYPE(grampcd::Message_send_ping);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_send_ping)

CEREAL_REGISTER_TYPE(grampcd::Message_send_two_coupling_states);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_send_two_coupling_states)

CEREAL_REGISTER_TYPE(grampcd::Message_send_convergenceFlag);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_send_convergenceFlag)

CEREAL_REGISTER_TYPE(grampcd::Message_send_number_of_active_couplings);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_send_number_of_active_couplings)

CEREAL_REGISTER_TYPE(grampcd::Message_acknowledge_received_optimization_info);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_acknowledge_received_optimization_info)

CEREAL_REGISTER_TYPE(grampcd::Message_send_flag_stopped_admm);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_send_flag_stopped_admm)

CEREAL_REGISTER_TYPE(grampcd::Message_send_flag_stopped_admm_coordinator);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_send_flag_stopped_admm_coordinator)

CEREAL_REGISTER_TYPE(grampcd::Message_send_flag_to_stop_admm);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::Message, grampcd::Message_send_flag_to_stop_admm)


CEREAL_REGISTER_DYNAMIC_INIT(message)