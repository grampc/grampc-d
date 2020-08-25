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

#include "../include/ssms_coupling_model.hpp"
#include "math.h"

SSMSCouplingModel::SSMSCouplingModel(const std::vector<typeRNum>& model_parameters, const std::string& name)
    : CouplingModel(2, 1, 2, 1, 0, 0,
        model_parameters,
        name)
{
    p1_ = model_parameters[0]; // m_i
    p2_ = model_parameters[1]; // c
}

dmpc::CouplingModelPtr SSMSCouplingModel::create(const std::vector<typeRNum>& model_parameters, const std::string& name)
{
	return std::shared_ptr<CouplingModel>(new SSMSCouplingModel(model_parameters, name));
}

void SSMSCouplingModel::ffct(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{
    const typeRNum dx = fabs( xj[0] - xi[0] );
    if (dx > 1e-6)
    {
        out[0] += 0.0;
        out[1] += p2_ / p1_ * (1.0 - 1.0 / dx) * (xj[0] - xi[0]);
    }
    else
    {
        out[0] += 0.0;
        out[1] += 0;
    }
}

void SSMSCouplingModel::dfdxi_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
    out[0] += -p2_ / p1_ * vec[1]; 
    out[1] += 0;
}

void SSMSCouplingModel::dfdui_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
    out[0] += 0.0;
}

void SSMSCouplingModel::dfdxj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
    out[0] += p2_ / p1_ * vec[1];
    out[1] += 0;
}

void SSMSCouplingModel::dfduj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
    out[0] += 0.0;
}
