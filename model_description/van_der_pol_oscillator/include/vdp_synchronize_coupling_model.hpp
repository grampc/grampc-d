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

#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/archives/binary.hpp>

class VDPSynchronizeCouplingModel : public grampcd::CouplingModel
{
public:
	VDPSynchronizeCouplingModel
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

	// model parameters
	typeRNum p1_;

	// cost parameters
	typeRNum Q_;

	/*
	* The following functions enable serializing the object.
	*/

	// A default constructor is required.
	VDPSynchronizeCouplingModel() {};

	// The serialize function is required.
	template<class Archive>
	void serialize(Archive& ar)
	{
		ar(
			// serialize member variables of this specific coupling model
			p1_, Q_,
			//serialize member variables of the general coupling model
			Nxi_, Nui_, Nxj_, Nuj_, Ngij_, Nhij_, model_parameters_, cost_parameters_, model_name_
		);
	}
};

// Register coupling model
CEREAL_REGISTER_TYPE(VDPSynchronizeCouplingModel);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::CouplingModel, VDPSynchronizeCouplingModel)