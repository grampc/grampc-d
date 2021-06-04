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

#include "../include/vdp_linear_coupling_model.hpp"

VDPLinearCouplingModel::VDPLinearCouplingModel(const std::vector<typeRNum>& model_parameters, const std::string& name)
	: CouplingModel(2, 1, 2, 1, 0, 0,
		model_parameters,
		name)
{
	p1_ = model_parameters[0];
}

dmpc::CouplingModelPtr VDPLinearCouplingModel::create(const std::vector<typeRNum>& model_parameters, const std::string& name)
{
	return std::shared_ptr<CouplingModel>(new VDPLinearCouplingModel(model_parameters, name));
}

void VDPLinearCouplingModel::ffct(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{
	out[0] += 0.0;
	out[1] += p1_ * (xj[0] - xi[0]);
}

void VDPLinearCouplingModel::dfdxi_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
	out[0] += -vec[1] * p1_;
	out[1] += 0.0;
}

void VDPLinearCouplingModel::dfdui_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
	out[0] += 0.0;
}

void VDPLinearCouplingModel::dfdxj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
	out[0] += vec[1] * p1_;
	out[1] += 0.0;
}

void VDPLinearCouplingModel::dfduj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
	out[0] += 0.0;
}