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

#include "../include/mobile_robot_coupling_model.hpp"
#include <cmath>

// Utility macros for integer powers
#define POW2(x) (x)*(x)
#define POW3(x) (x)*(x)*(x)
#define POW4(x) (x)*(x)*(x)*(x)

MobileRobotCouplingModel::MobileRobotCouplingModel
(
	const std::vector<typeRNum>& model_parameters,
	const std::vector<typeRNum>& cost_parameters,
	const std::string& name
)
	: CouplingModel(3, 2, 4, 3, 0, 0,
		model_parameters,
		cost_parameters,
		name)
{
	// Agent i should be follower with Nxi=3, Nui=2
	// Neighbor j should be leader with Nxj=4, Nuj=3
}

grampcd::CouplingModelPtr MobileRobotCouplingModel::create
(
	const std::vector<typeRNum>& model_parameters,
	const std::vector<typeRNum>& cost_parameters,
	const std::string& name
)
{
	return std::shared_ptr<CouplingModel>(new MobileRobotCouplingModel(model_parameters, cost_parameters, name));
}

void MobileRobotCouplingModel::ffct(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{
}

void MobileRobotCouplingModel::dfdxi_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
}

void MobileRobotCouplingModel::dfdui_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
}

void MobileRobotCouplingModel::dfdxj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
}

void MobileRobotCouplingModel::dfduj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
}

void MobileRobotCouplingModel::lfct(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{
    // longitudinal and lateral errors
    typeRNum elon =  cos(xj[2]) * (xi[0] - xj[0]) + sin(xj[2]) * (xi[1] - xj[1]) - model_parameters_[0];
    typeRNum elat = -sin(xj[2]) * (xi[0] - xj[0]) + cos(xj[2]) * (xi[1] - xj[1]) - model_parameters_[1];
    // non-quadratic cost function
    out[0] += cost_parameters_[0] * POW4(elon)
            + cost_parameters_[1] * POW2(elat)
            + cost_parameters_[2] * POW4(sin((xi[2] - xj[2]) / 2))
			+ cost_parameters_[3] * POW4(ui[0] - uj[0])
            + cost_parameters_[4] * POW4(ui[1] - uj[1]);
}

void MobileRobotCouplingModel::dldxi(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{
    // longitudinal and lateral errors
	typeRNum elon =  cos(xj[2]) * (xi[0] - xj[0]) + sin(xj[2]) * (xi[1] - xj[1]) - model_parameters_[0];
    typeRNum elat = -sin(xj[2]) * (xi[0] - xj[0]) + cos(xj[2]) * (xi[1] - xj[1]) - model_parameters_[1];
    // gradient of non-quadratic cost function
    out[0] += cost_parameters_[0] * 4 * POW3(elon) * cos(xj[2]) + cost_parameters_[1] * 2 * elat * -sin(xj[2]);
    out[1] += cost_parameters_[0] * 4 * POW3(elon) * sin(xj[2]) + cost_parameters_[1] * 2 * elat * cos(xj[2]);
    out[2] += cost_parameters_[2] * 4 * POW3(sin(xi[2] - xj[2]) / 2) * (cos((xi[2] - xj[2]) / 2) / 2);
}

void MobileRobotCouplingModel::dldui(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{
	out[0] += cost_parameters_[3] * 4 * POW3(ui[0] - uj[0]);
	out[1] += cost_parameters_[4] * 4 * POW3(ui[1] - uj[1]);
}

void MobileRobotCouplingModel::dldxj(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{
    // longitudinal and lateral errors
	typeRNum elon =  cos(xj[2]) * (xi[0] - xj[0]) + sin(xj[2]) * (xi[1] - xj[1]) - model_parameters_[0];
    typeRNum elat = -sin(xj[2]) * (xi[0] - xj[0]) + cos(xj[2]) * (xi[1] - xj[1]) - model_parameters_[1];
    // gradient of non-quadratic cost function
    out[0] += cost_parameters_[0] * 4 * POW3(elon) * -cos(xj[2]) + cost_parameters_[1] * 2 * elat * sin(xj[2]);
    out[1] += cost_parameters_[0] * 4 * POW3(elon) * -sin(xj[2]) + cost_parameters_[1] * 2 * elat * -cos(xj[2]);
    out[2] += cost_parameters_[0] * 4 * POW3(elon) * (-sin(xj[2]) * (xi[0] - xj[0]) + cos(xj[2]) * (xi[1] - xj[1]))
	        + cost_parameters_[1] * 2 *      elat  * (-cos(xj[2]) * (xi[0] - xj[0]) - sin(xj[2]) * (xi[1] - xj[1]))
			+ cost_parameters_[2] * 4 * POW3(sin(xi[2] - xj[2]) / 2) * -(cos((xi[2] - xj[2]) / 2) / 2);
}

void MobileRobotCouplingModel::dlduj(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{
	out[0] += cost_parameters_[3] * 4 * POW3(ui[0] - uj[0]) * (-1);
	out[1] += cost_parameters_[4] * 4 * POW3(ui[1] - uj[1]) * (-1);
}

void MobileRobotCouplingModel::Vfct(typeRNum* out, ctypeRNum T, ctypeRNum* xi, ctypeRNum* xj)
{
}

void MobileRobotCouplingModel::dVdxi(typeRNum* out, ctypeRNum T, ctypeRNum* xi, ctypeRNum* xj)
{
}

void MobileRobotCouplingModel::dVdxj(typeRNum* out, ctypeRNum T, ctypeRNum* xi, ctypeRNum* xj)
{
}

void MobileRobotCouplingModel::hfct(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj)
{
}

void MobileRobotCouplingModel::dhdxi_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
}

void MobileRobotCouplingModel::dhdui_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
}

void MobileRobotCouplingModel::dhdxj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
}

void MobileRobotCouplingModel::dhduj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec)
{
}