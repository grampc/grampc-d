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

#include "general_model_factory.hpp"

#include "dmpc/info/agent_info.hpp"
#include "dmpc/info/coupling_info.hpp"

#include "dmpc/model/agent_model.hpp"
#include "dmpc/model/coupling_model.hpp"

#include "water_tank/include/water_tank_agent_model.hpp"
#include "water_tank/include/water_tank_coupling_model.hpp"

#include "van_der_pol_oscillator/include/vdp_agent_model.hpp"
#include "van_der_pol_oscillator/include/vdp_linear_coupling_model.hpp"
#include "van_der_pol_oscillator/include/vdp_nonlinear_coupling_model.hpp"

#include "scalable_spring_mass_system/include/ssms_agent_model.hpp"
#include "scalable_spring_mass_system/include/ssms_coupling_model.hpp"

#include "scalable_spring_mass_system_2D/include/ssms2d_agent_model.hpp"
#include "scalable_spring_mass_system_2D/include/ssms2d_coupling_model.hpp"

#include "scalable_spring_mass_system_3D/include/ssms3d_agent_model.hpp"
#include "scalable_spring_mass_system_3D/include/ssms3d_coupling_model.hpp"

#include "oscillators/include/oscillators_agent_model.hpp"
#include "oscillators/include/oscillators_coupling_model.hpp"

#include "smart_grid/include/smartGrid_agentModel.hpp"
#include "smart_grid/include/smartGrid_couplingModel.hpp"

GeneralModelFactory::GeneralModelFactory(const dmpc::LoggingPtr& log) :
	log_(log)
{
	map_agentModels_["vdp_agentModel"] = VDPAgentModel::create;
	map_couplingModels_["vdp_linear_couplingModel"] = VDPLinearCouplingModel::create;
	map_couplingModels_["vdp_nonlinear_couplingModel"] = VDPNonlinearCouplingModel::create;

	map_agentModels_["water_tank_agentModel"] = WaterTankAgentModel::create;
	map_couplingModels_["water_tank_couplingModel"] = WaterTankCouplingModel::create;
	map_couplingModels_["vdp_linear_couplingModel"] = VDPLinearCouplingModel::create;

	map_agentModels_["ssms_agentModel"] = SSMSAgentModel::create;
	map_couplingModels_["ssms_couplingModel"] = SSMSCouplingModel::create;

	map_agentModels_["ssms2d_agentModel"] = SSMS2DAgentModel::create;
	map_couplingModels_["ssms2d_couplingModel"] = SSMS2DCouplingModel::create;

	map_agentModels_["ssms3d_agentModel"] = SSMS3DAgentModel::create;
	map_couplingModels_["ssms3d_couplingModel"] = SSMS3DCouplingModel::create;

	map_agentModels_["oscillators_agentModel"] = OscillatorsAgentModel::create;
	map_couplingModels_["oscillators_couplingModel"] = OscillatorsCouplingModel::create;

	map_agentModels_["smartGrid_agentModel"] = SmartGridAgentModel::create;
	map_couplingModels_["smartGrid_couplingModel"] = SmartGridCouplingModel::create;
}

dmpc::AgentModelPtr GeneralModelFactory::create_agentModel(const dmpc::AgentInfo& info) const
{
	auto iterator = map_agentModels_.find(info.model_name_);

	if (iterator == map_agentModels_.end())
		log_->print(dmpc::DebugType::Error) << "Invalid model name '" << info.model_name_ << "'" << std::endl;

	return (*(iterator->second))(info.model_parameters_, info.cost_parameters_, info.model_name_, log_);
}

dmpc::CouplingModelPtr GeneralModelFactory::create_couplingModel(const dmpc::CouplingInfo& info) const
{
	auto iterator = map_couplingModels_.find(info.model_name_);

	if (iterator == map_couplingModels_.end())
		log_->print(dmpc::DebugType::Error) << "Invalid model name '" << info.model_name_ << "'" << std::endl;

	return (*(iterator->second))(info.model_parameters_, info.model_name_);
}
