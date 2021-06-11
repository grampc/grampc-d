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

#include "../include/oscillators_coupling_model.hpp"

OscillatorsCouplingModel::OscillatorsCouplingModel
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
	m_ = model_parameters[0]; // tank area
	c_ = model_parameters[1]; // pipe area
}

dmpc::CouplingModelPtr OscillatorsCouplingModel::create
(
	const std::vector<typeRNum>& model_parameters,
	const std::vector<typeRNum>& cost_parameters,
	const std::string& name
)
{
	return std::shared_ptr<CouplingModel>(new OscillatorsCouplingModel(model_parameters, cost_parameters, name));
}

void OscillatorsCouplingModel::ffct(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{
	out[0] += 0;
	out[1] += c_ / m_ * xj[0] - c_ / m_ * xi[0];
}

void OscillatorsCouplingModel::dfdxi_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
	out[0] += -c_ / m_ * vec[1];
	out[1] += 0;
}

void OscillatorsCouplingModel::dfdui_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
	out[0] += 0.0;
}

void OscillatorsCouplingModel::dfdxj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
	out[0] += c_ / m_ * vec[1];
}

void OscillatorsCouplingModel::dfduj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
	out[0] += 0.0;
}

void OscillatorsCouplingModel::lfct(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{

}

void OscillatorsCouplingModel::dldxi(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{

}

void OscillatorsCouplingModel::dldui(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{

}

void OscillatorsCouplingModel::dldxj(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{

}

void OscillatorsCouplingModel::dlduj(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{

}

void OscillatorsCouplingModel::Vfct(typeRNum* out, ctypeRNum T, ctypeRNum* xi, ctypeRNum* xj)
{

}

void OscillatorsCouplingModel::dVdxi(typeRNum* out, ctypeRNum T, ctypeRNum* xi, ctypeRNum* xj)
{

}

void OscillatorsCouplingModel::dVdxj(typeRNum* out, ctypeRNum T, ctypeRNum* xi, ctypeRNum* xj)
{

}