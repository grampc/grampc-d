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

#ifndef COUPLING_MODEL_HPP
#define COUPLING_MODEL_HPP

#include "dmpc/util/types.hpp"

namespace dmpc
{

/*@brief Model of a coupling between two agents*/
class CouplingModel
{
public:

    CouplingModel(unsigned int Nxi, unsigned int Nui, unsigned int Nxj, unsigned int Nuj, unsigned int Ngij, unsigned int Nhij,
                  const std::vector<typeRNum>& model_parameters,
                  const std::string& model_name);

    virtual ~CouplingModel();

    /*Returns number of states of the agent.*/
    const unsigned int get_Nxi() const;
    /*Returns number of controls of the agent.*/
    const unsigned int get_Nui() const;
    /*Returns number of states of the neighbor.*/
    const unsigned int get_Nxj() const;
    /*Returns number of controls of the neighbor.*/
    const unsigned int get_Nuj() const;
    /*Returns number of equality constraints.*/
    const unsigned int get_Ngij() const;
    /*Returns number of inequality constraints.*/
    const unsigned int get_Nhij() const;

    /*Returns the model parameters.*/
    const std::vector<typeRNum> get_modelParameters() const;
    /*Returns the model name.*/
    const std::string get_modelName() const;

    /*Coupling dynamics f_{ij}(x_i, u_i, x_j, u_i, t)*/
    virtual void ffct(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj) = 0;
    /*Partial derivative of the coupling dynamics with respect to the agents states multiplied with Lagrangian multipliers.*/
	virtual void dfdxi_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec) = 0;
	/*Partial derivative of the coupling dynamics with respect to the agents controls multiplied with Lagrangian multipliers.*/
	virtual void dfdui_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec) = 0;
	/*Partial derivative of the coupling dynamics with respect to the neighbors states multiplied with Lagrangian multipliers.*/
	virtual void dfdxj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec) = 0;
	/*Partial derivative of the coupling dynamics with respect to the neighbors controls multiplied with Lagrangian multipliers.*/
    virtual void dfduj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec) = 0;

    /*Equality constraints g_{ij}(x_i, u_i, x_j, u_j, t) = 0*/
	virtual void gfct(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj) {}
	/*Partial derivative of the equality constraints with respect to the agents states multiplied with Lagrangian multipliers.*/
	virtual void dgdxi_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec) {}
	/*Partial derivative of the equality constraints with respect to the agents controls multiplied with Lagrangian multipliers.*/
	virtual void dgdui_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec) {}
	/*Partial derivative of the equality constraints with respect to the neighbors states multiplied with Lagrangian multipliers.*/
	virtual void dgdxj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec) {}
	/*Partial derivative of the equality constraints with respect to the neighbors controls multiplied with Lagrangian multipliers.*/
    virtual void dgduj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec) {}

	/*Inequality constraints h_{ij}(x_i, u_i, x_j, u_j, t) <= 0*/
	virtual void hfct(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj) {}
	/*Partial derivative of the inequality constraints with respect to the agents states multiplied with Lagrangian multipliers.*/
	virtual void dhdxi_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec) {}
	/*Partial derivative of the inequality constraints with respect to the agents controls multiplied with Lagrangian multipliers.*/
	virtual void dhdui_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec) {}
	/*Partial derivative of the inequality constraints with respect to the neighbors states multiplied with Lagrangian multipliers.*/
	virtual void dhdxj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec) {}
	/*Partial derivative of the inequality constraints with respect to the neighbors controls multiplied with Lagrangian multipliers.*/
    virtual void dhduj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec) {}

private:
    unsigned int Nxi_;
    unsigned int Nui_;
    unsigned int Nxj_;
    unsigned int Nuj_;
    unsigned int Ngij_;
    unsigned int Nhij_;

    std::vector<typeRNum> model_parameters_;
    std::string model_name_;
};

typedef std::shared_ptr<CouplingModel> CouplingModelPtr;

}

#endif // COUPLING_MODEL_HPP
