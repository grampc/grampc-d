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

#include "../include/ssms3d_coupling_model.hpp"

SSMS3DCouplingModel::SSMS3DCouplingModel(const std::vector<typeRNum>& model_parameters, const std::string& name)
    : CouplingModel(6, 3, 6, 3, 0, 0,
        model_parameters,
        name)
{
    p1_ = model_parameters[0]; // m_i
    p2_ = model_parameters[1]; // c

    d0_ = 1;
}

dmpc::CouplingModelPtr SSMS3DCouplingModel::create(const std::vector<typeRNum>& model_parameters, const std::string& name)
{
	return std::shared_ptr<CouplingModel>(new SSMS3DCouplingModel(model_parameters, name));
}

void SSMS3DCouplingModel::ffct(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{
    typeRNum d = sqrt( ( xi[0] - xj[0] )*( xi[0] - xj[0] ) + ( xi[2] - xj[2] )*( xi[2] - xj[2] ) + ( xi[4] - xj[4] )*( xi[4] - xj[4] ) );
    typeRNum dx = ( xj[0] - xi[0] );
    typeRNum dy = ( xj[2] - xi[2] );
    typeRNum dz = ( xj[4] - xi[4] );
    out[0] += 0.0;
    out[2] += 0.0;
    out[4] += 0.0;
    if( abs(d) < 1e-12 )
    {
        out[1] += 0.0;
        out[3] += 0.0;
        out[5] += 0.0;
    }
    else
    {
        out[1] += p2_ / p1_ * ( 1 - d0_ / d ) * dx;
        out[3] += p2_ / p1_ * ( 1 - d0_ / d ) * dy;
        out[5] += p2_ / p1_ * ( 1 - d0_ / d ) * dz;
    }

}

void SSMS3DCouplingModel::dfdxi_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
    typeRNum d = sqrt( ( xi[0] - xj[0] )*( xi[0] - xj[0] ) + ( xi[2] - xj[2] )*( xi[2] - xj[2] ) + ( xi[4] - xj[4] )*( xi[4] - xj[4] ) );
    typeRNum dx = ( xj[0] - xi[0] );
    typeRNum dy = ( xj[2] - xi[2] );
    typeRNum dz = ( xj[4] - xi[4] );
    out[1] += 0;
    out[3] += 0;
    out[5] += 0;
    if( abs(d) < 1e-12 )
    {
        out[0] += 0;
        out[2] += 0;
        out[4] += 0;
    }
    else
    {
        out[0] += p2_ / p1_ * ( d0_ / d - 1 - d0_ / ( d*d*d ) * dx*dx ) * vec[1] - p2_ / p1_ * d0_ / ( d*d*d ) * dx * dy * vec[3] - p2_ / p1_ * d0_ / ( d*d*d ) * dx * dz * vec[5];
        out[2] += -p2_ / p1_ * d0_ / ( d*d*d ) * dx * dy * vec[1] + p2_ / p1_ * ( d0_ / d - 1 - d0_ / ( d*d*d ) * dy*dy ) * vec[3] - p2_ / p1_ * d0_ / ( d*d*d ) * dy * dz * vec[5];
        out[4] += -p2_ / p1_ * d0_ / ( d*d*d ) * dx * dz * vec[1] - p2_ / p1_ * d0_ / ( d*d*d ) * dy * dz * vec[3] + p2_ / p1_ * ( d0_ / d - 1 - d0_ / ( d*d*d ) * dz*dz ) * vec[5];
    }

}

void SSMS3DCouplingModel::dfdui_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
    out[0] += 0.0;
    out[1] += 0.0;
    out[2] += 0.0;
}

void SSMS3DCouplingModel::dfdxj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
    typeRNum d = sqrt( ( xi[0] - xj[0] )*( xi[0] - xj[0] ) + ( xi[2] - xj[2] )*( xi[2] - xj[2] ) + ( xi[4] - xj[4] )*( xi[4] - xj[4] ) );
    typeRNum dx = ( xj[0] - xi[0] );
    typeRNum dy = ( xj[2] - xi[2] );
    typeRNum dz = ( xj[4] - xi[4] );
    out[1] += 0;
    out[3] += 0;
    out[5] += 0;
    if( abs(d) < 1e-12 )
    {
        out[0] += 0;
        out[2] += 0;
        out[4] += 0;
    }
    else
    {
        out[0] += p2_ / p1_ * ( d0_ / ( d*d*d ) * dx*dx - d0_ / d + 1 ) * vec[1] + p2_ / p1_ * d0_ / ( d*d*d ) * dx * dy * vec[3] + p2_ / p1_ * d0_ / ( d*d*d ) * dx * dz * vec[5];
        out[2] += p2_ / p1_ * d0_ / ( d*d*d ) * dx * dy * vec[1] + p2_ / p1_ * ( d0_ / ( d*d*d ) * dy*dy - d0_ / d + 1 ) * vec[3] + p2_ / p1_ * d0_ / ( d*d*d ) * dy * dz * vec[1];
        out[4] += p2_ / p1_ * d0_ / ( d*d*d ) * dx * dz * vec[1] + p2_ / p1_ * d0_ / ( d*d*d ) * dy * dz * vec[3] + p2_ / p1_ * ( d0_ / ( d*d*d ) * dz*dz - d0_ / d + 1 ) * vec[5];
    }

}

void SSMS3DCouplingModel::dfduj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
    out[0] += 0.0;
    out[1] += 0.0;
    out[2] += 0.0;
}
