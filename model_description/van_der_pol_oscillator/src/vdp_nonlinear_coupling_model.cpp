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

#include "../include/vdp_nonlinear_coupling_model.hpp"

VDPNonlinearCouplingModel::VDPNonlinearCouplingModel(const std::vector<typeRNum>& model_parameters, const std::string& name)
	: CouplingModel(2, 1, 2, 1, 0, 0,
		model_parameters,
		name)
{
    p1_ = model_parameters[0];
}

dmpc::CouplingModelPtr VDPNonlinearCouplingModel::create(const std::vector<typeRNum>& model_parameters, const std::string& name)
{
    return std::shared_ptr<CouplingModel>(new VDPNonlinearCouplingModel(model_parameters, name));
}

void VDPNonlinearCouplingModel::ffct(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{
    out[0] += 0.0;
    out[1] += p1_ * xj[0] * xj[1];
}

void VDPNonlinearCouplingModel::dfdxi_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
    out[0] += 0.0;
    out[1] += 0.0;
}

void VDPNonlinearCouplingModel::dfdui_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
    out[0] += 0.0;
}

void VDPNonlinearCouplingModel::dfdxj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
    out[0] += vec[1] * p1_ * xj[1];
    out[1] += vec[1] * p1_ * xj[0];
}

void VDPNonlinearCouplingModel::dfduj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
    out[0] += 0.0;
}
