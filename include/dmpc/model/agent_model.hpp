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

#ifndef AGENT_MODEL_HPP
#define AGENT_MODEL_HPP

#include "dmpc/util/types.hpp"
#include "dmpc/util/logging.hpp"

namespace dmpc
{

/*@brief Agent models describe a part of the OCP regarding a specific agent.*/
class AgentModel
{
public:

    AgentModel(unsigned int Nxi, unsigned int Nui, unsigned int Ngi, unsigned int Nhi,
		const std::vector<double>& umin, const std::vector<double>& umax,
		const std::vector<typeRNum>& model_parameters,
		const std::vector<typeRNum>& cost_parameters,
        const std::string& model_name,
        const LoggingPtr& log);

    virtual ~AgentModel();

    /*Returns number of states of the agent.*/
    const unsigned int get_Nxi() const;
    /*Returns number of controls of the agent.*/
    const unsigned int get_Nui() const;
    /*Returns number of equality constraints.*/
    const unsigned int get_Ngi() const;
    /*Returns number of inequality constraints.*/
    const unsigned int get_Nhi() const;

    /*Returns the control limits.*/
    const void get_controlLimits(std::vector<typeRNum>& umin, std::vector<typeRNum>& umax);

    /*Agent dynamics f(x_i, u_i, t)*/
    virtual void ffct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u) = 0;
    /*Partial derivate of the agent dynamics with respect to states multiplied with Lagrangian multipliers.*/
	virtual void dfdx_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec) = 0;
	/*Partial derivate of the agent dynamics with respect to controls multiplied with Lagrangian multipliers.*/
    virtual void dfdu_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec) = 0;

    /*Cost function l_i(x_i, u_i, t)*/
	virtual void lfct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes) = 0;
	/*Partial derivate of the cost function with respect to states multiplied with Lagrangian multipliers.*/
	virtual void dldx(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes) = 0;
	/*Partial derivate of the cost function with respect to controls multiplied with Lagrangian multipliers.*/
    virtual void dldu(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* xdes) = 0;

    /*Terminating cost V_i(x_i(T))*/
	virtual void Vfct(typeRNum* out, ctypeRNum T, ctypeRNum* x, ctypeRNum* xdes) = 0;
	/*Partial derivate of the terminating cost with respect to states multiplied with Lagrangian multipliers.*/
    virtual void dVdx(typeRNum* out, ctypeRNum T, ctypeRNum* x, ctypeRNum* xdes) = 0;

    /*Equality constraints g_i(x_i, u_i, t) = 0*/
	virtual void gfct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u) {}
	/*Partial derivate of the equality constraints with respect to states multiplied with Lagrangian multipliers.*/
	virtual void dgdx_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec) {}
	/*Partial derivate of the equality constraints with respect to controls multiplied with Lagrangian multipliers.*/
    virtual void dgdu_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec) {}

	/*Inequality constraints h_i(x_i, u_i, t) <= 0*/
	virtual void hfct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u) {}
	/*Partial derivate of the inequality constraints with respect to states multiplied with Lagrangian multipliers.*/
	virtual void dhdx_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec) {}
	/*Partial derivate of the inequality constraints with respect to controls multiplied with Lagrangian multipliers.*/
    virtual void dhdu_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* vec) {}

    /*Returns model parameters.*/
    const std::vector<typeRNum> get_modelParameters() const;
    /*Returns cost parameters.*/
    const std::vector<typeRNum> get_costParameters() const;
    /*Returns the model name.*/
    const std::string get_modelName() const;

private:
    unsigned int Nxi_;
    unsigned int Nui_;
    unsigned int Ngi_;
    unsigned int Nhi_;
    std::vector<typeRNum> umin_;
    std::vector<typeRNum> umax_;

    // These variables are for being able to recreate the model
    std::vector<typeRNum> model_parameters_;
    std::vector<typeRNum> cost_parameters_;
    std::string model_name_;

    LoggingPtr log_;
};

typedef std::shared_ptr<AgentModel> AgentModelPtr;

}

#endif // AGENT_MODEL_HPP
