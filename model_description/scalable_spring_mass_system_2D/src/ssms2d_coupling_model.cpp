/* This file is part of GRAMPC-D - (https://github.com/grampc-d/grampc-d.git)
 *
 * GRAMPC-D -- A software framework for distributed model predictive control (DMPC)
 * 
 *
 * Copyright 2023 by Daniel Burk, Maximilian Pierer von Esch, Andreas Voelz, Knut Graichen
 * All rights reserved.
 *
 * GRAMPC-D is distributed under the BSD-3-Clause license, see LICENSE.txt
 *
 */

#include "../include/ssms2d_coupling_model.hpp"
#include <cmath>

SSMS2DCouplingModel::SSMS2DCouplingModel
(
	const std::vector<typeRNum>& model_parameters,
	const std::vector<typeRNum>& cost_parameters,
	const std::string& name
)
	: CouplingModel(4, 2, 4, 2, 0, 0,
		model_parameters,
		cost_parameters,
		name)
{
	p1_ = model_parameters[0]; // m_i
	p2_ = model_parameters[1]; // c

	d0_ = 1;
	dmin = 0.2;
}

grampcd::CouplingModelPtr SSMS2DCouplingModel::create
(
	const std::vector<typeRNum>& model_parameters,
	const std::vector<typeRNum>& cost_parameters,
	const std::string& name
)
{
	return std::shared_ptr<CouplingModel>(new SSMS2DCouplingModel(model_parameters, cost_parameters, name));
}

void SSMS2DCouplingModel::ffct(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{
	typeRNum d = std::sqrt((xi[0] - xj[0]) * (xi[0] - xj[0]) + (xi[2] - xj[2]) * (xi[2] - xj[2]));
	typeRNum dx = (xj[0] - xi[0]);
	typeRNum dy = (xj[2] - xi[2]);
	if (d < dmin)
	{
		out[0] += 0.0;
		out[1] += 0.0;
		out[2] += 0.0;
		out[3] += 0.0;
	}
	else
	{
		out[0] += 0.0;
		out[1] += p2_ / p1_ * (1 - d0_ / d) * dx;
		out[2] += 0.0;
		out[3] += p2_ / p1_ * (1 - d0_ / d) * dy;
	}

}

void SSMS2DCouplingModel::dfdxi_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
	typeRNum d = std::sqrt((xi[0] - xj[0]) * (xi[0] - xj[0]) + (xi[2] - xj[2]) * (xi[2] - xj[2]));
	typeRNum dx = (xj[0] - xi[0]);
	typeRNum dy = (xj[2] - xi[2]);
	if (d < dmin)
	{
		out[0] += 0.0;
		out[1] += 0.0;
		out[2] += 0.0;
		out[3] += 0.0;
	}
	else
	{
		out[0] += p2_ / p1_ * (d0_ / d - 1 - d0_ / (d * d * d) * dx * dx) * vec[1] - p2_ / p1_ * d0_ / (d * d * d) * dx * dy * vec[3];
		out[1] += 0;
		out[2] += -p2_ / p1_ * d0_ / (d * d * d) * dx * dy * vec[1] + p2_ / p1_ * (d0_ / d - 1 - d0_ / (d * d * d) * dy * dy) * vec[3];
		out[3] += 0;
	}
}

void SSMS2DCouplingModel::dfdui_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
	out[0] += 0.0;
	out[1] += 0.0;
}

void SSMS2DCouplingModel::dfdxj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
	typeRNum d = std::sqrt((xi[0] - xj[0]) * (xi[0] - xj[0]) + (xi[2] - xj[2]) * (xi[2] - xj[2]));
	typeRNum dx = (xj[0] - xi[0]);
	typeRNum dy = (xj[2] - xi[2]);
	if (d < dmin)
	{
		out[0] += 0.0;
		out[1] += 0.0;
		out[2] += 0.0;
		out[3] += 0.0;
	}
	else
	{
		out[0] += p2_ / p1_ * (d0_ / (d * d * d) * dx * dx - d0_ / d + 1) * vec[1] + p2_ / p1_ * d0_ / (d * d * d) * dx * dy * vec[3];
		out[1] += 0;
		out[2] += p2_ / p1_ * d0_ / (d * d * d) * dx * dy * vec[1] + p2_ / p1_ * (d0_ / (d * d * d) * dy * dy - d0_ / d + 1) * vec[3];
		out[3] += 0;
	}
}

void SSMS2DCouplingModel::dfduj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
	out[0] += 0.0;
	out[1] += 0.0;
}

void SSMS2DCouplingModel::lfct(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{

}

void SSMS2DCouplingModel::dldxi(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{

}

void SSMS2DCouplingModel::dldui(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{

}

void SSMS2DCouplingModel::dldxj(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{

}

void SSMS2DCouplingModel::dlduj(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{

}

void SSMS2DCouplingModel::Vfct(typeRNum* out, ctypeRNum T, ctypeRNum* xi, ctypeRNum* xj)
{

}

void SSMS2DCouplingModel::dVdxi(typeRNum* out, ctypeRNum T, ctypeRNum* xi, ctypeRNum* xj)
{

}

void SSMS2DCouplingModel::dVdxj(typeRNum* out, ctypeRNum T, ctypeRNum* xi, ctypeRNum* xj)
{

}