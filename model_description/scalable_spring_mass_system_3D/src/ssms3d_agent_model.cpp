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

#include "../include/ssms3d_agent_model.hpp"

SSMS3DAgentModel::SSMS3DAgentModel(
	const std::vector<typeRNum>& model_parameters,
	const std::vector<typeRNum>& cost_parameters,
	const std::string& name,
	const grampcd::LoggingPtr& log)
	: AgentModel(6, 3, 0, 0, { -100, -100, -100 }, { 100, 100, 100 },
		model_parameters,
		cost_parameters,
		name,
		log)
{
	g_ = 9.81;
	p1_ = model_parameters[0]; //m_i
	p2_ = model_parameters[1]; //m_i
	p3_ = model_parameters[2];
	P_.push_back(cost_parameters[0]);
	P_.push_back(cost_parameters[1]);
	P_.push_back(cost_parameters[2]);
	P_.push_back(cost_parameters[3]);
	P_.push_back(cost_parameters[4]);
	P_.push_back(cost_parameters[5]);
	Q_.push_back(cost_parameters[6]);
	Q_.push_back(cost_parameters[7]);
	Q_.push_back(cost_parameters[8]);
	Q_.push_back(cost_parameters[9]);
	Q_.push_back(cost_parameters[10]);
	Q_.push_back(cost_parameters[11]);
	R_.push_back(cost_parameters[12]);
	R_.push_back(cost_parameters[13]);
	R_.push_back(cost_parameters[14]);
}

grampcd::AgentModelPtr SSMS3DAgentModel::create(
	const std::vector<typeRNum>& model_parameters,
	const std::vector<typeRNum>& cost_parameters,
	const std::string& name,
	const grampcd::LoggingPtr& log)
{
	return std::shared_ptr<AgentModel>(new SSMS3DAgentModel(model_parameters, cost_parameters, name, log));
}

void SSMS3DAgentModel::ffct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u)
{
	out[0] += x[1];
	out[1] += p3_ / p1_ * u[0];
	out[2] += x[3];
	out[3] += p3_ / p1_ * u[1];
	out[4] += x[5];
	out[5] += p3_ / p1_ * u[2] - p2_ * g_;
}

void SSMS3DAgentModel::dfdx_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec)
{
	out[0] += 0;
	out[1] += vec[0];
	out[2] += 0;
	out[3] += vec[2];
	out[4] += 0;
	out[5] += vec[4];
}

void SSMS3DAgentModel::dfdu_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec)
{
	out[0] += p3_ / p1_ * vec[1];
	out[1] += p3_ / p1_ * vec[3];
	out[2] += p3_ / p1_ * vec[5];
}

void SSMS3DAgentModel::lfct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes)
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

void SSMS3DAgentModel::dldx(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes)
{
	for (unsigned int i = 0; i < get_Nxi(); ++i)
	{
		out[i] += Q_[i] * 2.0 * (x[i] - xdes[i]);
	}
}

void SSMS3DAgentModel::dldu(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes)
{
	for (unsigned int i = 0; i < get_Nui(); ++i)
	{
		out[i] += R_[i] * 2.0 * u[i];
	}
}

void SSMS3DAgentModel::Vfct(typeRNum* out, ctypeRNum T, ctypeRNum* x, ctypeRNum* xdes)
{
	for (unsigned int i = 0; i < get_Nxi(); ++i)
	{
		out[0] += P_[i] * (x[i] - xdes[i]) * (x[i] - xdes[i]);
	}
}

void SSMS3DAgentModel::dVdx(typeRNum* out, ctypeRNum T, ctypeRNum* x, ctypeRNum* xdes)
{
	for (unsigned int i = 0; i < get_Nxi(); ++i)
	{
		out[i] += P_[i] * 2.0 * (x[i] - xdes[i]);
	}
}