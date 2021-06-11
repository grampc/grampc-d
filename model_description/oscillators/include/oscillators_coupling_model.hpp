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

#pragma once

#include "grampcd/model/coupling_model.hpp"

class OscillatorsCouplingModel : public grampcd::CouplingModel
{
public:
	OscillatorsCouplingModel
	(
		const std::vector<typeRNum>& model_parameters,
		const std::vector<typeRNum>& cost_parameters,
		const std::string& name
	);

	static grampcd::CouplingModelPtr create
	(
		const std::vector<typeRNum>& model_parameters,
		const std::vector<typeRNum>& cost_parameters,
		const std::string& name
	);

	void ffct(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj) override;
	void dfdxi_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec) override;
	void dfdui_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec) override;
	void dfdxj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec) override;
	void dfduj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec) override;

	void lfct(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj) override;
	void dldxi(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj) override;
	void dldui(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj) override;
	void dldxj(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj) override;
	void dlduj(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj) override;

	void Vfct(typeRNum* out, ctypeRNum T, ctypeRNum* xi, ctypeRNum* xj) override;
	void dVdxi(typeRNum* out, ctypeRNum T, ctypeRNum* xi, ctypeRNum* xj) override;
	void dVdxj(typeRNum* out, ctypeRNum T, ctypeRNum* xi, ctypeRNum* xj) override;

private:
	typeRNum m_;
	typeRNum c_;
};