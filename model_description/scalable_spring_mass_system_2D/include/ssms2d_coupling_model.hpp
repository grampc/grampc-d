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

#ifndef SSMS_2D_couplingModel_HPP
#define SSMS_2D_couplingModel_HPP

#include "dmpc/model/coupling_model.hpp"

class SSMS2DCouplingModel : public dmpc::CouplingModel
{
public:
	SSMS2DCouplingModel(const std::vector<typeRNum>& model_parameters, const std::string& name);

	static dmpc::CouplingModelPtr create(const std::vector<typeRNum>& model_parameters, const std::string& name);

    virtual void ffct(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj) override;

    virtual void dfdxi_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec) override;

    virtual void dfdui_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec) override;

    virtual void dfdxj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec) override;

    virtual void dfduj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec) override;

private:
    typeRNum p1_;
    typeRNum p2_;
    typeRNum d0_;
    typeRNum dmin;
};

#endif // SSMS_2D_couplingModel_HPP
