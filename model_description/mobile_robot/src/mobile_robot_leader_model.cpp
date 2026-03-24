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

#include "../include/mobile_robot_leader_model.hpp"

// Utility macros for integer powers
#define POW2(x) (x)*(x)
#define POW3(x) (x)*(x)*(x)
#define POW4(x) (x)*(x)*(x)*(x)

MobileRobotLeaderModel::MobileRobotLeaderModel(
	const std::vector<typeRNum>& model_parameters,
	const std::vector<typeRNum>& cost_parameters,
	const std::string& name,
	const grampcd::LoggingPtr& log)
	: AgentModel(4, 3, 0, 0, {-2, -2, -1}, {2, 2, 1},
		model_parameters,
		cost_parameters,
		name,
		log)
{
}

grampcd::AgentModelPtr MobileRobotLeaderModel::create(
	const std::vector<typeRNum>& model_parameters,
	const std::vector<typeRNum>& cost_parameters,
	const std::string& name,
	const grampcd::LoggingPtr& log)
{
	return std::shared_ptr<AgentModel>(new MobileRobotLeaderModel(model_parameters, cost_parameters, name, log));
}

void MobileRobotLeaderModel::ffct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u)
{
	out[0] += u[0] * cos(x[2]);
	out[1] += u[0] * sin(x[2]);
	out[2] += u[1];
    out[3] += u[2];
}

void MobileRobotLeaderModel::dfdx_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec)
{
	out[0] += 0;
	out[1] += 0;
	out[2] += u[0] * (vec[1] * cos(x[2]) - vec[0] * sin(x[2]));
    out[3] += 0;
}

void MobileRobotLeaderModel::dfdu_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec)
{
    out[0] += vec[0] * cos(x[2]) + vec[1] * sin(x[2]);
    out[1] += vec[2];
    out[2] += vec[3];
}

void MobileRobotLeaderModel::lfct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes)
{	
    // desired path for position and orientation
	typeRNum p[3] = {cos(x[3]), sin(x[3]), M_PI_2 + x[3]};
    // longitudinal and lateral errors
    typeRNum elon =  cos(p[2]) * (x[0] - p[0]) + sin(p[2]) * (x[1] - p[1]);
    typeRNum elat = -sin(p[2]) * (x[0] - p[0]) + cos(p[2]) * (x[1] - p[1]);
    // non-quadratic cost function
    out[0] += cost_parameters_[0] * POW4(elon)
            + cost_parameters_[1] * POW2(elat)
            + cost_parameters_[2] * POW4(sin(x[2] - p[2]) / 2)
            + cost_parameters_[3] * POW4(u[0] - u[2])
            + cost_parameters_[4] * POW4(u[1] - u[2])
            + cost_parameters_[5] * POW4(u[2] - 1.0);
}

void MobileRobotLeaderModel::dldx(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes)
{
    // desired path for position and orientation
	typeRNum p[3] = {cos(x[3]), sin(x[3]), M_PI_2 + x[3]};
    // longitudinal and lateral errors
    typeRNum elon =  cos(p[2]) * (x[0] - p[0]) + sin(p[2]) * (x[1] - p[1]);
    typeRNum elat = -sin(p[2]) * (x[0] - p[0]) + cos(p[2]) * (x[1] - p[1]);
    // gradient of non-quadratic cost function
    out[0] += cost_parameters_[0] * 4 * POW3(elon) * cos(p[2]) + cost_parameters_[1] * 2 * elat * -sin(p[2]);
    out[1] += cost_parameters_[0] * 4 * POW3(elon) * sin(p[2]) + cost_parameters_[1] * 2 * elat *  cos(p[2]);
    out[2] += cost_parameters_[2] * 4 * POW3(sin((x[2] - p[2]) / 2)) * (cos((x[2] - p[2]) / 2) / 2);
    out[3] += cost_parameters_[0] * 4 * POW3(elon) * (-sin(p[2]) * x[0] + cos(p[2]) * x[1])
            + cost_parameters_[1] * 2 *      elat  * (-cos(p[2]) * x[0] - sin(p[2]) * x[1])
            + cost_parameters_[2] * 4 * POW3(sin((x[2] - p[2]) / 2)) * -(cos((x[2] - p[2]) / 2) / 2);
}

void MobileRobotLeaderModel::dldu(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes)
{
    // gradient of non-quadratic cost function
    out[0] += cost_parameters_[3] * 4 * POW3(u[0] - u[2]);
    out[1] += cost_parameters_[4] * 4 * POW3(u[1] - u[2]);
    out[2] += cost_parameters_[3] * 4 * POW3(u[0] - u[2]) * (-1)
            + cost_parameters_[4] * 4 * POW3(u[1] - u[2]) * (-1)
            + cost_parameters_[5] * 4 * POW3(u[2] - 1.0);
}

void MobileRobotLeaderModel::Vfct(typeRNum* out, ctypeRNum T, ctypeRNum* x, ctypeRNum* xdes)
{
}

void MobileRobotLeaderModel::dVdx(typeRNum* out, ctypeRNum T, ctypeRNum* x, ctypeRNum* xdes)
{
}

void MobileRobotLeaderModel::hfct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u)
{
}

void MobileRobotLeaderModel::dhdx_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec)
{
}

void MobileRobotLeaderModel::dhdu_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec)
{
}
