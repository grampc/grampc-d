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

#include "../include/smartGrid_agentModel.hpp"

SmartGridAgentModel::SmartGridAgentModel(
    const std::vector<typeRNum>& model_parameters,
    const std::vector<typeRNum>& cost_parameters,
	const std::string& name,
	const LoggingPtr& log)
    : AgentModel(2, 1, 0, 0, {-1e6}, {1e6},
                 model_parameters,
                 cost_parameters,
                 name,
                 log)
{
    I_ = model_parameters[0]; // inertia
    Omega_ = model_parameters[1]; // frequency
    kappa_ = model_parameters[2]; // friction term
    P0_ = model_parameters[3];
    p_ = model_parameters[4]; // control

    P_.push_back(cost_parameters[0]); // terminal state weight
    P_.push_back(cost_parameters[1]); // terminal state weight
    Q_.push_back(cost_parameters[2]); // integral state weight
    Q_.push_back(cost_parameters[3]); // integral state weight
    R_.push_back(cost_parameters[4]); // integral control weight
}

dmpc::AgentModelPtr SmartGridAgentModel::create(
	const std::vector<typeRNum>& model_parameters,
	const std::vector<typeRNum>& cost_parameters,
	const std::string& name,
	const LoggingPtr& log)
{
	return dmpc::AgentModelPtr(new SmartGridAgentModel(model_parameters, cost_parameters, name, log));
}

void SmartGridAgentModel::ffct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u)
{
    out[0] += x[1];
    out[1] += 1 / (I_ * Omega_) * ( p_ * u[0] + P0_ - kappa_ * Omega_*Omega_) - 2*kappa_ / I_ * x[1];
}

void SmartGridAgentModel::dfdx_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec)
{
    out[0] += 0.0;
    out[1] += vec[0] - 2*kappa_ / I_ * vec[1];
}

void SmartGridAgentModel::dfdu_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec)
{
    out[0] += p_ / ( I_ * Omega_ ) * vec[1];
}

void SmartGridAgentModel::lfct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes)
{
    for(unsigned int i = 0; i < get_Nxi(); ++i)
        out[0] += Q_[i] * (x[i] - xdes[i]) * (x[i] - xdes[i]);

    for(unsigned int i = 0; i < get_Nui(); ++i)
        out[0] += R_[i] * u[i] * u[i];
}

void SmartGridAgentModel::dldx(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes)
{
    for(unsigned int i = 0; i < get_Nxi(); ++i)
        out[i] += Q_[i] * 2.0 * (x[i] - xdes[i]);
}

void SmartGridAgentModel::dldu(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes)
{
    for(unsigned int i = 0; i < get_Nui(); ++i)
        out[i] += R_[i] * 2.0 * u[i];
}

void SmartGridAgentModel::Vfct(typeRNum* out, ctypeRNum T, ctypeRNum* x, ctypeRNum* xdes)
{
    for(unsigned int i = 0; i < get_Nxi(); ++i)
        out[0] += P_[i] * (x[i] - xdes[i]) * (x[i] - xdes[i]);
}

void SmartGridAgentModel::dVdx(typeRNum* out, ctypeRNum T, ctypeRNum* x, ctypeRNum* xdes)
{
    for (unsigned int i = 0; i < get_Nxi(); ++i)
        out[i] += P_[i] * 2.0 * (x[i] - xdes[i]);
}



