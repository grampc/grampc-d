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

#include "../include/vdp_agent_model.hpp"

VDPAgentModel::VDPAgentModel(
	const std::vector<typeRNum>& model_parameters,
	const std::vector<typeRNum>& cost_parameters,
	const std::string& name,
	const grampcd::LoggingPtr& log)
	: AgentModel(2, 1, 0, 0, { -1.0 }, { 1.0 },
		model_parameters,
		cost_parameters,
		name,
		log)
{
	p1_ = model_parameters[0];
	p2_ = model_parameters[1];
	p3_ = model_parameters[2];
	P_.push_back(cost_parameters[0]);
	P_.push_back(cost_parameters[1]);
	Q_.push_back(cost_parameters[2]);
	Q_.push_back(cost_parameters[3]);
	R_.push_back(cost_parameters[4]);
}

grampcd::AgentModelPtr VDPAgentModel::create(const std::vector<typeRNum>& model_parameters,
	const std::vector<typeRNum>& cost_parameters,
	const std::string& name,
	const grampcd::LoggingPtr& log)
{
	return std::shared_ptr<AgentModel>(new VDPAgentModel(model_parameters, cost_parameters, name, log));
}

void VDPAgentModel::ffct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u)
{
	out[0] += x[1];
	out[1] += p1_ * (1.0 - p2_ * x[0] * x[0]) * x[1] - p3_ * x[0] + u[0];
}

void VDPAgentModel::dfdx_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec)
{
	out[0] += vec[1] * (-p1_ * p2_ * 2.0 * x[0] * x[1] - p3_);
	out[1] += vec[0] + vec[1] * p1_ * (1.0 - p2_ * x[0] * x[0]);
}

void VDPAgentModel::dfdu_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec)
{
	out[0] += vec[1];
}

void VDPAgentModel::lfct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes)
{
	for (unsigned int i = 0; i < get_Nxi(); ++i)
		out[0] += Q_[i] * (x[i] - xdes[i]) * (x[i] - xdes[i]);

	for (unsigned int i = 0; i < get_Nui(); ++i)
		out[0] += R_[i] * u[i] * u[i];
}

void VDPAgentModel::dldx(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes)
{
	for (unsigned int i = 0; i < get_Nxi(); ++i)
		out[i] += Q_[i] * 2.0 * (x[i] - xdes[i]);
}

void VDPAgentModel::dldu(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes)
{
	for (unsigned int i = 0; i < get_Nui(); ++i)
		out[i] += R_[i] * 2.0 * u[i];
}

void VDPAgentModel::Vfct(typeRNum* out, ctypeRNum T, ctypeRNum* x, ctypeRNum* xdes)
{
	for (unsigned int i = 0; i < get_Nxi(); ++i)
		out[0] += P_[i] * (x[i] - xdes[i]) * (x[i] - xdes[i]);
}

void VDPAgentModel::dVdx(typeRNum* out, ctypeRNum T, ctypeRNum* x, ctypeRNum* xdes)
{
	for (unsigned int i = 0; i < get_Nxi(); ++i)
		out[i] += 2 * P_[i] * (x[i] - xdes[i]);
}