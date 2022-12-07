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

#include "grampcd/util/class_forwarding.hpp"

namespace grampcd
{

    /*@brief Model of a coupling between two agents*/
    class CouplingModel
    {
    public:

        CouplingModel
        (
            unsigned int Nxi, unsigned int Nui, 
            unsigned int Nxj, unsigned int Nuj, 
            unsigned int Ngij, unsigned int Nhij,
			const std::vector<typeRNum>& model_parameters,
			const std::vector<typeRNum>& cost_parameters,
            const std::string& model_name
        );

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

		/*Second-order partial derivative of the coupling dynamics w.r.t. to the neighbors states multiplied with the Lagrangian multipliers*/
		virtual void dfdxjdxj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec){}
		/*Second-order partial derivative of the coupling dynamics w.r.t. to the neighbors controls multiplied with the Lagrangian multipliers*/
		virtual void dfdujduj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec) {}
		/*Second-order mixed partial derivative of the coupling dynamics w.r.t. to the neighbors states and controls multiplied with the Lagrangian multipliers*/
		virtual void dfdxjduj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec) {}


		/*Coupling cost function l_{ij}(x_i, u_i, x_j, u_i, t)*/
		virtual void lfct(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj) = 0;
		/*Partial derivative of the cost function with respect to the agents states.*/
		virtual void dldxi(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj) = 0;
		/*Partial derivative of the cost function with respect to the agents controls.*/
		virtual void dldui(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj) = 0;
		/*Partial derivative of the cost function with respect to the neighbors states.*/
		virtual void dldxj(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj) = 0;
		/*Partial derivative of the cost function with respect to the neighbors controls.*/
		virtual void dlduj(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj) = 0;

		/*Second-order partial derivative of the coupling cost function w.r.t.to the neighbors states */
		virtual void dldxjdxj(typeRNum * out, typeRNum t, ctypeRNum * xi, ctypeRNum * ui, ctypeRNum * xj, ctypeRNum * uj) {}
		/*Second-order partial derivative of the coupling function w.r.t. to the neighbors controls*/
		virtual void dldujduj(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj) {}
		/*Second-order mixed partial derivative of the coupling function w.r.t. to the neighbors states and controls*/
		virtual void dldxjduj(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj) {}

		/*Coupled terminating cost V_i(x_i(T))*/
		virtual void Vfct(typeRNum* out, ctypeRNum T, ctypeRNum* xi, ctypeRNum* xj) = 0;
		/*Partial derivate of the coupled terminating cost with respect to the agents states.*/
		virtual void dVdxi(typeRNum* out, ctypeRNum T, ctypeRNum* xi, ctypeRNum* xj) = 0;
		/*Partial derivate of the coupled terminating cost with respect to the neighbors states.*/
		virtual void dVdxj(typeRNum* out, ctypeRNum T, ctypeRNum* xi, ctypeRNum* xj) = 0;
		/*Partial derivate of the coupled terminating cost with respect to the neighbors states.*/
		virtual void dVdxjdxj(typeRNum* out, ctypeRNum T, ctypeRNum* xi, ctypeRNum* xj) {}

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

		/*Second-order partial derivative of the coupled equality constraints w.r.t.to the neighbors states multiplied with the Lagrangian multipliers */
		virtual void dgdxjdxj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec) {}
		/*Second-order partial derivative of the coupled equality constraints w.r.t. to the neighbors controls multiplied with the Lagrangian multipliers*/
		virtual void dgdujduj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec) {}
		/*Second-order mixed partial derivative of the coupled equality constraints w.r.t. to the neighbors states and controls multiplied with the Lagrangian multipliers*/
		virtual void dgdxjduj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec) {}

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

		/*Second-order partial derivative of the coupled inequality constraints w.r.t.to the neighbors states multiplied with the Lagrangian multipliers */
		virtual void dhdxjdxj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec) {}
		/*Second-order partial derivative of the coupled inequality constraints w.r.t. to the neighbors controls multiplied with the Lagrangian multipliers*/
		virtual void dhdujduj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec) {}
		/*Second-order mixed partial derivative of the coupled inequality constraints w.r.t. to the neighbors states and controls multiplied with the Lagrangian multipliers*/
		virtual void dhdxjduj_vec(typeRNum* out, typeRNum t, ctypeRNum* xi, ctypeRNum* ui, ctypeRNum* xj, ctypeRNum* uj, ctypeRNum* vec) {}

        unsigned int Nxi_;
        unsigned int Nui_;
        unsigned int Nxj_;
        unsigned int Nuj_;
        unsigned int Ngij_;
        unsigned int Nhij_;

		std::vector<typeRNum> model_parameters_;
		std::vector<typeRNum> cost_parameters_;
		std::string model_name_;

		CouplingModel() {};
    };

}