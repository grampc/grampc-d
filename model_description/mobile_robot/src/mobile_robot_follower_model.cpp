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

#include "../include/mobile_robot_follower_model.hpp"

MobileRobotFollowerModel::MobileRobotFollowerModel(
	const std::vector<typeRNum>& model_parameters,
	const std::vector<typeRNum>& cost_parameters,
	const std::string& name,
	const grampcd::LoggingPtr& log)
	: AgentModel(3, 2, 0, 0, {-2, -2}, {2, 2},
		model_parameters,
		cost_parameters,
		name,
		log)
{
}

grampcd::AgentModelPtr MobileRobotFollowerModel::create(
	const std::vector<typeRNum>& model_parameters,
	const std::vector<typeRNum>& cost_parameters,
	const std::string& name,
	const grampcd::LoggingPtr& log)
{
	return std::shared_ptr<AgentModel>(new MobileRobotFollowerModel(model_parameters, cost_parameters, name, log));
}

void MobileRobotFollowerModel::ffct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u)
{
	out[0] += u[0] * cos(x[2]);
	out[1] += u[0] * sin(x[2]);
	out[2] += u[1];
}

void MobileRobotFollowerModel::dfdx_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec)
{
	out[0] += 0;
	out[1] += 0;
	out[2] += u[0] * (vec[1] * cos(x[2]) - vec[0] * sin(x[2]));
}

void MobileRobotFollowerModel::dfdu_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec)
{
    out[0] += vec[0] * cos(x[2]) + vec[1] * sin(x[2]);
    out[1] += vec[2];
}

void MobileRobotFollowerModel::lfct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes)
{
}

void MobileRobotFollowerModel::dldx(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes)
{
}

void MobileRobotFollowerModel::dldu(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes)
{
}

void MobileRobotFollowerModel::Vfct(typeRNum* out, ctypeRNum T, ctypeRNum* x, ctypeRNum* xdes)
{
}

void MobileRobotFollowerModel::dVdx(typeRNum* out, ctypeRNum T, ctypeRNum* x, ctypeRNum* xdes)
{
}

void MobileRobotFollowerModel::hfct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u)
{
}

void MobileRobotFollowerModel::dhdx_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec)
{
}

void MobileRobotFollowerModel::dhdu_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec)
{
}
