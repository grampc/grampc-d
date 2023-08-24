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

#include "../include/water_tank_agent_model.hpp"

WaterTankAgentModel::WaterTankAgentModel(
	const std::vector<typeRNum>& model_parameters,
	const std::vector<typeRNum>& cost_parameters,
	const std::string& name,
	const grampcd::LoggingPtr& log)
	: AgentModel(1, 1, 0, 1, { 0.0 }, { 0.2 },
		model_parameters,
		cost_parameters,
		name,
		log)
{
	Ai_ = model_parameters[0]; // tank area
	ci_ = model_parameters[1]; // inflow
	di_ = model_parameters[2]; // outflow

	P_.push_back(cost_parameters[0]); // terminal state weight
	Q_.push_back(cost_parameters[1]); // integral state weight
	R_.push_back(cost_parameters[2]); // integral control weight
}

grampcd::AgentModelPtr WaterTankAgentModel::create(
	const std::vector<typeRNum>& model_parameters,
	const std::vector<typeRNum>& cost_parameters,
	const std::string& name,
	const grampcd::LoggingPtr& log)
{
	return std::shared_ptr<AgentModel>(new WaterTankAgentModel(model_parameters, cost_parameters, name, log));
}

void WaterTankAgentModel::ffct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u)
{
	out[0] += (ci_ * u[0] - di_) / Ai_;
}

void WaterTankAgentModel::dfdx_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec)
{
	out[0] += 0.0;
}

void WaterTankAgentModel::dfdu_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec)
{
	out[0] += vec[0] * ci_ / Ai_;
}

void WaterTankAgentModel::lfct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes)
{
	for (unsigned int i = 0; i < get_Nxi(); ++i)
	{
		out[0] += Q_[i] * (x[i] - xdes[i]) * (x[i] - xdes[i]);
	}
	for (unsigned int i = 0; i < get_Nui(); ++i)
	{
		out[0] += R_[i] * u[i] * u[i];
	}
}

void WaterTankAgentModel::dldx(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes)
{
	for (unsigned int i = 0; i < get_Nxi(); ++i)
	{
		out[i] += Q_[i] * 2.0 * (x[i] - xdes[i]);
	}
}

void WaterTankAgentModel::dldu(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes)
{
	for (unsigned int i = 0; i < get_Nui(); ++i)
	{
		out[i] += R_[i] * 2.0 * u[i];
	}
}

void WaterTankAgentModel::Vfct(typeRNum* out, ctypeRNum T, ctypeRNum* x, ctypeRNum* xdes)
{
	for (unsigned int i = 0; i < get_Nxi(); ++i)
	{
		out[0] += P_[i] * (x[i] - xdes[i]) * (x[i] - xdes[i]);
	}
}

void WaterTankAgentModel::dVdx(typeRNum* out, ctypeRNum T, ctypeRNum* x, ctypeRNum* xdes)
{
	for (unsigned int i = 0; i < get_Nxi(); ++i)
	{
		out[i] += P_[i] * 2.0 * (x[i] - xdes[i]);
	}
}

void WaterTankAgentModel::hfct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u)
{
	out[0] += x[0] - 3;
}

void WaterTankAgentModel::dhdx_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec)
{
	 out[0] += vec[0];
}

void WaterTankAgentModel::dhdu_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec)
{
	out[0] += 0.0;
}


