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

#include "../include/ssms_coupling_model.hpp"
#include <cmath>

SSMSCouplingModel::SSMSCouplingModel
(
	const std::vector<typeRNum>& model_parameters,
	const std::vector<typeRNum>& cost_parameters,
	const std::string& name
)
	: CouplingModel(2, 1, 2, 1, 0, 0,
		model_parameters,
		cost_parameters,
		name)
{
	p1_ = model_parameters[0]; // m_i
	p2_ = model_parameters[1]; // c
}

grampcd::CouplingModelPtr SSMSCouplingModel::create
(
	const std::vector<typeRNum>& model_parameters,
	const std::vector<typeRNum>& cost_parameters,
	const std::string& name
)
{
	return std::shared_ptr<CouplingModel>(new SSMSCouplingModel(model_parameters, cost_parameters, name));
}

void SSMSCouplingModel::ffct(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{
	const typeRNum dx = std::fabs(xj[0] - xi[0]);
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

void SSMSCouplingModel::lfct(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{

}

void SSMSCouplingModel::dldxi(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{

}

void SSMSCouplingModel::dldui(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{

}

void SSMSCouplingModel::dldxj(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{

}

void SSMSCouplingModel::dlduj(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{

}

void SSMSCouplingModel::Vfct(typeRNum* out, ctypeRNum T, ctypeRNum* xi, ctypeRNum* xj)
{

}

void SSMSCouplingModel::dVdxi(typeRNum* out, ctypeRNum T, ctypeRNum* xi, ctypeRNum* xj)
{

}

void SSMSCouplingModel::dVdxj(typeRNum* out, ctypeRNum T, ctypeRNum* xi, ctypeRNum* xj)
{

}