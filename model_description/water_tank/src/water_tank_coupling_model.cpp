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

#include "../include/water_tank_coupling_model.hpp"
#include <cmath>

const int signum(typeRNum a)
{
    return (a >= 0) ? 1 : -1;
}

WaterTankCouplingModel::WaterTankCouplingModel(const std::vector<typeRNum>& model_parameters, const std::string& name)
    : CouplingModel(1, 1, 1, 1, 0, 0,
                    model_parameters,
                    name)
{
    Ai_ = model_parameters[0]; // tank area
    aij_ = model_parameters[1]; // pipe area
    g_ = 9.81; // gravitation
    eps_ = 0.01; // region where Toricelli is replaced by linear approximation

    // parameters for polynomial fitting
    const typeRNum q = aij_ / Ai_ * signum(eps_) * std::sqrt(2 * g_ * std::abs(eps_));
    const typeRNum dqdeps = -(std::sqrt(2) * aij_ * g_ * signum(eps_)) / (2 * Ai_ * std::sqrt(g_ * std::abs(eps_)));
    poly_param1_ = (3 * q) / (2 * eps_) - dqdeps / 2;
    poly_param2_ = dqdeps / (2 * eps_*eps_) - q / (2 * eps_*eps_*eps_);
}

dmpc::CouplingModelPtr WaterTankCouplingModel::create(const std::vector<typeRNum>& model_parameters, const std::string& name)
{
    return std::shared_ptr<CouplingModel>(new WaterTankCouplingModel(model_parameters, name));
}

void WaterTankCouplingModel::ffct(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{
    const typeRNum dx = xj[0] - xi[0];
    if(std::abs(dx) < eps_)
        out[0] += poly_param2_ * dx * dx * dx + poly_param1_ * dx;
    else
        out[0] += (aij_ / Ai_) * signum(dx) * std::sqrt(2.0 * g_ * std::abs(dx));
}

void WaterTankCouplingModel::dfdxi_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
    const typeRNum dx = xj[0] - xi[0];
    if(std::abs(dx) < eps_)

        out[0] += vec[0] * (-poly_param1_ - 3 * poly_param2_ * dx * dx);
    else
        out[0] += vec[0] * (aij_ / Ai_) * g_ / std::sqrt(2.0 * g_ * std::abs(dx)) * (-1.0);
}

void WaterTankCouplingModel::dfdui_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
    out[0] += 0.0;
}

void WaterTankCouplingModel::dfdxj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
    const typeRNum dx = xj[0] - xi[0];
    if(std::abs(dx) < eps_)
        out[0] += vec[0] * (poly_param1_ + 3 * poly_param2_ * dx*dx);
    else
        out[0] += vec[0] * (aij_ / Ai_) * g_ / std::sqrt(2.0 * g_ * std::abs(dx));
}

void WaterTankCouplingModel::dfduj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
    out[0] += 0.0;
}
