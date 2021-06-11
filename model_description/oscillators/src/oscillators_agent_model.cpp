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

#include "../include/oscillators_agent_model.hpp"

OscillatorsAgentModel::OscillatorsAgentModel
(
	const std::vector<typeRNum>& model_parameters,
	const std::vector<typeRNum>& cost_parameters,
	const std::string& name,
	const grampcd::LoggingPtr& log
)
	: AgentModel(2, 1, 0, 1, { -100 }, { 100 },
		model_parameters,
		cost_parameters,
		name,
		log)
{
	mi_ = model_parameters[0]; // tank area
	ci_ = model_parameters[1]; // inflow
	di_ = model_parameters[2]; // outflow
	pi_ = model_parameters[3];

	P_.push_back(cost_parameters[0]); // terminal state weight
	Q_.push_back(cost_parameters[1]); // integral state weight
	R_.push_back(cost_parameters[2]); // integral control weight
}

grampcd::AgentModelPtr OscillatorsAgentModel::create
(
	const std::vector<typeRNum>& model_parameters,
	const std::vector<typeRNum>& cost_parameters,
	const std::string& name,
	const grampcd::LoggingPtr& log
)
{
	return std::shared_ptr<AgentModel>(new OscillatorsAgentModel(model_parameters, cost_parameters, name, log));
}

void OscillatorsAgentModel::ffct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u)
{
	out[0] += x[1];
	out[1] += -ci_ / mi_ * x[0] - di_ / mi_ * x[1] + pi_ / mi_ * u[0];
}

void OscillatorsAgentModel::dfdx_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec)
{
	out[0] += -ci_ / mi_ * vec[1];
	out[1] += vec[0] - di_ / mi_ * vec[1];
}

void OscillatorsAgentModel::dfdu_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec)
{
	out[0] += pi_ / mi_ * vec[1];
}

void OscillatorsAgentModel::lfct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes)
{
	for (unsigned int i = 0; i < get_Nxi(); ++i)
	{
		out[0] += Q_[i] * x[i] * x[i];
	}
	for (unsigned int i = 0; i < get_Nui(); ++i)
	{
		out[0] += R_[i] * u[i] * u[i];
	}
}

void OscillatorsAgentModel::dldx(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes)
{
	for (unsigned int i = 0; i < get_Nxi(); ++i)
	{
		out[i] += Q_[i] * 2.0 * x[i];
	}
}

void OscillatorsAgentModel::dldu(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes)
{
	for (unsigned int i = 0; i < get_Nui(); ++i)
	{
		out[i] += R_[i] * 2.0 * u[i];
	}
}

void OscillatorsAgentModel::Vfct(typeRNum* out, ctypeRNum T, ctypeRNum* x, ctypeRNum* xdes)
{
	for (unsigned int i = 0; i < get_Nxi(); ++i)
	{
		out[0] += P_[i] * x[i] * x[i];
	}
}

void OscillatorsAgentModel::dVdx(typeRNum* out, ctypeRNum T, ctypeRNum* x, ctypeRNum* xdes)
{
	for (unsigned int i = 0; i < get_Nxi(); ++i)
	{
		out[i] += P_[i] * 2.0 * x[i];
	}
}

void OscillatorsAgentModel::hfct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u)
{
	out[0] += 0;//x[1] - 0.25;
}

void OscillatorsAgentModel::dhdx_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec)
{
	out[0] += 0;
	out[1] += 0;//vec[0];
}

void OscillatorsAgentModel::dhdu_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec)
{
	out[0] += 0;
	out[1] += 0;
}
