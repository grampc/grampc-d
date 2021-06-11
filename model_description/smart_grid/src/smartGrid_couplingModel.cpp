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

#include "../include/smartGrid_couplingModel.hpp"
#include <math.h>

 // https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}

SmartGridCouplingModel::SmartGridCouplingModel
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
	I_ = model_parameters[0]; // inertia
	Omega_ = model_parameters[1]; // frequency
	P_max_ij_ = model_parameters[2];
}

dmpc::CouplingModelPtr SmartGridCouplingModel::create
(
	const std::vector<typeRNum>& model_parameters,
	const std::vector<typeRNum>& cost_parameters,
	const std::string& name
)
{
	return dmpc::CouplingModelPtr(new SmartGridCouplingModel(model_parameters, cost_parameters, name));
}

void SmartGridCouplingModel::ffct(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{
	out[0] += 0.0;
	out[1] += P_max_ij_ / (I_ * Omega_) * sin(xi[0] - xj[0]);
}

void SmartGridCouplingModel::dfdxi_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
	out[0] += P_max_ij_ / (I_ * Omega_) * vec[1] * cos(xi[0] - xj[0]);
	out[1] += 0.0;
}

void SmartGridCouplingModel::dfdui_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
	out[0] += 0.0;
}

void SmartGridCouplingModel::dfdxj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
	out[0] += -P_max_ij_ / (I_ * Omega_) * vec[1] * cos(xi[0] - xj[0]);
	out[1] += 0.0;
}

void SmartGridCouplingModel::dfduj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
	out[0] += 0.0;
}

void SmartGridCouplingModel::lfct(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{

}

void SmartGridCouplingModel::dldxi(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{

}

void SmartGridCouplingModel::dldui(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{

}

void SmartGridCouplingModel::dldxj(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{

}

void SmartGridCouplingModel::dlduj(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{

}

void SmartGridCouplingModel::Vfct(typeRNum* out, ctypeRNum T, ctypeRNum* xi, ctypeRNum* xj)
{

}

void SmartGridCouplingModel::dVdxi(typeRNum* out, ctypeRNum T, ctypeRNum* xi, ctypeRNum* xj)
{

}

void SmartGridCouplingModel::dVdxj(typeRNum* out, ctypeRNum T, ctypeRNum* xi, ctypeRNum* xj)
{

}