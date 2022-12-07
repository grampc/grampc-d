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

#pragma once

#include "grampcd/model/agent_model.hpp"

#include "grampcd/util/logging.hpp"

#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/archives/binary.hpp>

class SSMS2DAgentModel : public grampcd::AgentModel
{
public:
	SSMS2DAgentModel(
		const std::vector<typeRNum>& model_parameters,
		const std::vector<typeRNum>& cost_parameters,
		const std::string& name,
		const grampcd::LoggingPtr& log);

	static grampcd::AgentModelPtr create(
		const std::vector<typeRNum>& model_parameters,
		const std::vector<typeRNum>& cost_parameters,
		const std::string& name,
		const grampcd::LoggingPtr& log);

	virtual void ffct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u) override;

	virtual void dfdx_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec) override;

	virtual void dfdu_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec) override;

	virtual void lfct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes) override;

	virtual void dldx(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes) override;

	virtual void dldu(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes) override;

	virtual void Vfct(typeRNum* out, ctypeRNum T, ctypeRNum* x, ctypeRNum* xdes) override;

	virtual void dVdx(typeRNum* out, ctypeRNum T, ctypeRNum* x, ctypeRNum* xdes) override;

	// model parameters
	typeRNum p1_;
	typeRNum p2_;

	// cost parameters
	std::vector<typeRNum> P_;
	std::vector<typeRNum> Q_;
	std::vector<typeRNum> R_;

	/*
	* The following functions enable serializing the object.
	*/

	// A default constructor is required.
	SSMS2DAgentModel() {};

	// The serialize function is required.
	template<class Archive>
	void serialize(Archive& ar)
	{
		ar(
			// serialize member variables of this specific agent model
			p1_, p2_, P_, Q_, R_,
			//serialize member variables of the general agent model
			Nxi_, Nui_, Ngi_, Nhi_, umin_, umax_, model_parameters_, cost_parameters_, model_name_, log_
		);
	}
};

// Register agent model
CEREAL_REGISTER_TYPE(SSMS2DAgentModel);
CEREAL_REGISTER_POLYMORPHIC_RELATION(grampcd::AgentModel, SSMS2DAgentModel)