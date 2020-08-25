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

#include "../include/ssms2d_agent_model.hpp"

SSMS2DAgentModel::SSMS2DAgentModel(
    const std::vector<typeRNum>& model_parameters,                             
    const std::vector<typeRNum>& cost_parameters,
	const std::string& name,
	const LoggingPtr& log)
    : AgentModel(4, 2, 0, 0, { -10, -10 }, { 10, 10 },
        model_parameters,
        cost_parameters,
        name,
        log)
{
    p1_ = model_parameters[0]; //m_i
    p2_ = model_parameters[1];
    P_.push_back(cost_parameters[0]);
    P_.push_back(cost_parameters[1]);
    P_.push_back(cost_parameters[2]);
    P_.push_back(cost_parameters[3]);
    Q_.push_back(cost_parameters[4]);
    Q_.push_back(cost_parameters[5]);
    Q_.push_back(cost_parameters[6]);
    Q_.push_back(cost_parameters[7]);
    R_.push_back(cost_parameters[8]);
    R_.push_back(cost_parameters[9]);
}

dmpc::AgentModelPtr SSMS2DAgentModel::create(
	const std::vector<typeRNum>& model_parameters,
	const std::vector<typeRNum>& cost_parameters,
	const std::string& name,
	const LoggingPtr& log)
{
	return std::shared_ptr<AgentModel>(new SSMS2DAgentModel(model_parameters, cost_parameters, name, log));
}

void SSMS2DAgentModel::ffct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u)
{
    out[0] += x[1];
    out[1] += p2_ / p1_ * u[0];
    out[2] += x[3];
    out[3] += p2_ / p1_ * u[1];
}

void SSMS2DAgentModel::dfdx_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec)
{
     out[0] += 0;
    out[1] += vec[0];
    out[2] += 0;
    out[3] += vec[2];
}

void SSMS2DAgentModel::dfdu_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec)
{
    out[0] += p2_ / p1_ * vec[1];
    out[1] += p2_ / p1_ * vec[3];
}

void SSMS2DAgentModel::lfct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes)
{
    for(unsigned int i = 0; i < get_Nxi(); ++i)
    {
        out[0] += Q_[i] * ( x[i] - xdes[i] ) * ( x[i] - xdes[i] );
    }
    for(unsigned int i = 0; i < get_Nui(); ++i)
    {
        out[0] += R_[i] * u[i] * u[i];
    }
}

void SSMS2DAgentModel::dldx(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes)
{
    for(unsigned int i = 0; i < get_Nxi(); ++i)
    {
        out[i] += Q_[i] * 2.0 * ( x[i] - xdes[i] );
    }
}

void SSMS2DAgentModel::dldu(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes)
{
    for(unsigned int i = 0; i < get_Nui(); ++i)
    {
        out[i] += R_[i] * 2.0 * u[i];
    }
}

void SSMS2DAgentModel::Vfct(typeRNum* out, ctypeRNum T, ctypeRNum* x, ctypeRNum* xdes)
{
    for(unsigned int i = 0; i < get_Nxi(); ++i)
    {
        out[0] += P_[i] * ( x[i] - xdes[i] ) * ( x[i] - xdes[i] );
    }
}

void SSMS2DAgentModel::dVdx(typeRNum* out, ctypeRNum T, ctypeRNum* x, ctypeRNum* xdes)
{
    for(unsigned int i = 0; i < get_Nxi(); ++i)
    {
        out[i] += P_[i] * 2.0 * ( x[i] - xdes[i] );
    }
}
