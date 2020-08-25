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

#ifndef GENERAL_MODEL_FACTORY_HPP
#define GENERAL_MODEL_FACTORY_HPP

#include "dmpc/model/model_factory.hpp"

class GeneralModelFactory : public dmpc::ModelFactory
{
public:
	using create_agentModel_functionPtr = dmpc::AgentModelPtr(*)(const std::vector<typeRNum>&, const std::vector<typeRNum>&, const std::string&, const LoggingPtr&);
	using create_couplingModel_functionPtr = dmpc::CouplingModelPtr(*)(const std::vector<typeRNum>&, const std::string&);

	GeneralModelFactory(const LoggingPtr& log);

	virtual dmpc::AgentModelPtr create_agentModel(const dmpc::AgentInfo& info) const override;

	virtual dmpc::CouplingModelPtr create_couplingModel(const dmpc::CouplingInfo& info) const override;

private:
	LoggingPtr log_;
	std::map<std::string, create_agentModel_functionPtr> map_agentModels_;
	std::map<std::string, create_couplingModel_functionPtr> map_couplingModels_;
};

#endif // GENERAL_MODEL_FACTORY_HPP
