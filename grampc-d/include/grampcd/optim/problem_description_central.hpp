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

#include <problem_description.hpp>

#include "grampcd/util/class_forwarding.hpp"

#include "grampcd/state/agent_state.hpp"

namespace grampcd
{

    /**
     * @brief Description of the centralized optimization problem for solution with GRAMPC.
     */
    class ProblemDescriptionCentral : public grampc::ProblemDescription
    {
    public:
        ProblemDescriptionCentral(const std::vector<AgentPtr>& agents);

        /*Returns the mapping for x_indices*/
	    const int get_x_index(int agent_id) const;
	    /*Returns the mapping for u_indices*/
        const int get_u_index( int agent_id ) const;

	    /*Set dimensions of the OCP*/
        virtual void ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT) override;

        /*Dynamics*/
        virtual void ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) override;
        /*Partial derivate of the dynamics with respect to states*/
	    virtual void dfdx_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* vec, ctypeRNum* u, ctypeRNum* p) override;
	    /*Partial derivate of the dynamics with respect to controls*/
        virtual void dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p) override;

        /*Cost function*/
	    virtual void lfct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* p, ctypeRNum* xdes, ctypeRNum* udes) override;
	    /*Partial derivate of the cost function with respect to states*/
	    virtual void dldx(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* p, ctypeRNum* xdes, ctypeRNum* udes) override;
	    /*Partial derivate of the cost function with respect to controls*/
        virtual void dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes) override;

        /*Terminal cost*/
	    virtual void Vfct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* p, ctypeRNum* xdes) override;
	    /*Partial derivate of the terminal cost with respect to states*/
        virtual void dVdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes) override;

        /*Equality constraints*/
	    virtual void gfct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* p) override;
	    /*Partial derivate of the equality constraints with respect to states*/
	    virtual void dgdx_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* p, ctypeRNum* vec) override;
	    /*Partial derivate of the equality constraints with respect to controls*/
        virtual void dgdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec) override;

	    /*Inequality constraints*/
	    virtual void hfct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* p) override;
	    /*Partial derivate of the inequality constraints with respect to states*/
	    virtual void dhdx_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* p, ctypeRNum* vec) override;
	    /*Partial derivate of the inequality constraints with respect to controls*/
        virtual void dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec) override;

    private:
        std::vector<AgentPtr> agents_;

        std::vector<int> x_index_ = std::vector<int>(0, 0);
        std::vector<int> u_index_ = std::vector<int>(0, 0);

        int Nx_ = 0;
        int Nu_ = 0;
        int Ng_ = 0;
        int Nh_ = 0;

        AgentState desired_;
    };

}