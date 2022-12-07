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

#include "../include/ssms_agent_model.hpp"

SSMSAgentModel::SSMSAgentModel(
	const std::vector<typeRNum>& model_parameters,
	const std::vector<typeRNum>& cost_parameters,
	const std::string& name,
	const grampcd::LoggingPtr& log)
	: AgentModel(2, 1, 0, 0, { -100 }, { 100 },
		model_parameters,
		cost_parameters,
		name,
		log)
{
	m_ = model_parameters[0]; //m_i
	c_ = model_parameters[1]; //m_i
	d_ = 1;
	P_.push_back(cost_parameters[0]);
	P_.push_back(cost_parameters[1]);
	Q_.push_back(cost_parameters[2]);
	Q_.push_back(cost_parameters[3]);
	R_.push_back(cost_parameters[4]);
}

grampcd::AgentModelPtr SSMSAgentModel::create(
	const std::vector<typeRNum>& model_parameters,
	const std::vector<typeRNum>& cost_parameters,
	const std::string& name,
	const grampcd::LoggingPtr& log)
{
	return std::shared_ptr<AgentModel>(new SSMSAgentModel(model_parameters, cost_parameters, name, log));
}

void SSMSAgentModel::ffct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u)
{
	out[0] += x[1];
	out[1] += -d_ / m_ * x[1] + c_ / m_ * u[0];
}

void SSMSAgentModel::dfdx_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec)
{
	out[0] += 0;
	out[1] += vec[0] - d_ / m_ * vec[1];
}

void SSMSAgentModel::dfdu_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec)
{
	out[0] += c_ / m_ * vec[1];
}

void SSMSAgentModel::lfct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes)
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

void SSMSAgentModel::dldx(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes)
{
	for (unsigned int i = 0; i < get_Nxi(); ++i)
	{
		out[i] += Q_[i] * 2.0 * (x[i] - xdes[i]);
	}
}

void SSMSAgentModel::dldu(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes)
{
	for (unsigned int i = 0; i < get_Nui(); ++i)
	{
		out[i] += R_[i] * 2.0 * u[i];
	}
}

void SSMSAgentModel::Vfct(typeRNum* out, ctypeRNum T, ctypeRNum* x, ctypeRNum* xdes)
{
	for (unsigned int i = 0; i < get_Nxi(); ++i)
	{
		out[0] += P_[i] * (x[i] - xdes[i]) * (x[i] - xdes[i]);
	}
}

void SSMSAgentModel::dVdx(typeRNum* out, ctypeRNum T, ctypeRNum* x, ctypeRNum* xdes)
{
	for (unsigned int i = 0; i < get_Nxi(); ++i)
	{
		out[i] += P_[i] * 2.0 * (x[i] - xdes[i]);
	}
}