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

#include "grampcd/model/model_factory.hpp"

#include "grampcd/util/class_forwarding.hpp"
#include "grampcd/util/logging.hpp"

class GeneralModelFactory : public grampcd::ModelFactory
{
public:
	using create_agentModel_functionPtr = grampcd::AgentModelPtr(*)(const std::vector<typeRNum>&, const std::vector<typeRNum>&, const std::string&, const grampcd::LoggingPtr&);
	using create_couplingModel_functionPtr = grampcd::CouplingModelPtr(*)(const std::vector<typeRNum>&, const std::vector<typeRNum>&, const std::string&);

	GeneralModelFactory(const grampcd::LoggingPtr& log);

	virtual grampcd::AgentModelPtr create_agentModel(const grampcd::AgentInfo& info) const override;

	virtual grampcd::CouplingModelPtr create_couplingModel(const grampcd::CouplingInfo& info) const override;

private:
	const grampcd::LoggingPtr log_;
	std::map<std::string, create_agentModel_functionPtr> map_agentModels_;
	std::map<std::string, create_couplingModel_functionPtr> map_couplingModels_;
};