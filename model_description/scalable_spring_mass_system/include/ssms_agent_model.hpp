/* This file is part of GRAMPC-D - (https://github.com/DanielBurk/GRAMPC-D.git)
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

#ifndef SSMS_agentModel_HPP
#define SSMS_agentModel_HPP

#include "dmpc/model/agent_model.hpp"

class SSMSAgentModel : public dmpc::AgentModel
{
public:
    SSMSAgentModel(
        const std::vector<typeRNum>& model_parameters,
        const std::vector<typeRNum>& cost_parameters,
		const std::string& name,
		const LoggingPtr& log);

	static dmpc::AgentModelPtr create(
		const std::vector<typeRNum>& model_parameters,
		const std::vector<typeRNum>& cost_parameters,
		const std::string& name,
		const LoggingPtr& log);

    virtual void ffct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u) override;

    virtual void dfdx_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec) override;

    virtual void dfdu_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec) override;

    virtual void lfct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes) override;

    virtual void dldx(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes) override;

    virtual void dldu(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes) override;

    virtual void Vfct(typeRNum* out, ctypeRNum T, ctypeRNum* x, ctypeRNum* xdes) override;

    virtual void dVdx(typeRNum* out, ctypeRNum T, ctypeRNum* x, ctypeRNum* xdes) override;

private:
    typeRNum m_;
    typeRNum c_;
    typeRNum d_;
    std::vector<typeRNum> P_;
    std::vector<typeRNum> Q_;
    std::vector<typeRNum> R_;
};

#endif // SSMS_agentModel_HPP
