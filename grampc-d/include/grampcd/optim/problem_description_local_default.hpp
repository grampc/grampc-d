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
#include "grampcd/util/types.hpp"

#include "grampcd/state/agent_state.hpp"
#include "grampcd/state/coupling_state.hpp"
#include "grampcd/state/multiplier_state.hpp"
#include "grampcd/state/penalty_state.hpp"

namespace grampcd
{

	/**
	 * @brief Description of the local optimization problem for solution with GRAMPC.
	 */
	class ProblemDescriptionLocalDefault : public grampc::ProblemDescription
	{
	public:
		ProblemDescriptionLocalDefault(Agent* agent);

		/*Returns the mapping of local copies x_{ji} to vector of controls.*/
		const std::vector<int>& get_u_index_xji() const;
		/*Returns the mapping of local copies u_{ji} to vector of controls.*/
		const std::vector<int>& get_u_index_uji() const;

		/*Returns the position of the local copy x_{ji} inside the vector of controls.*/
		const int get_u_index_xji(int agent_id) const;
		/*Returns the position of the local copy u_{ji} inside the vector of controls.*/
		const int get_u_index_uji(int agent_id) const;
		/*Returns the position of the local copy v_{ji} inside the vector of controls.*/
		const int get_u_index_vji(int agent_id) const;
		/*Returns the position of the local copy x_{ji} inside the vector of states.*/
		const int get_x_index_xji(int agent_id) const;

		/*Set dimensions of the OCP*/
		virtual void ocp_dim(typeInt* Nx, typeInt* Nu, typeInt* Np, typeInt* Ng, typeInt* Nh, typeInt* NgT, typeInt* NhT) override;

		/*Dynamics*/
		virtual void ffct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* p) override;
		/*Partial derivate of the dynamics with respect to states*/
		virtual void dfdx_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* vec, ctypeRNum* u, ctypeRNum* p) override;
		/*Partial derivate of the dynamics with respect to controls*/
		virtual void dfdu_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* vec, ctypeRNum* u, ctypeRNum* p) override;

		/*Cost function*/
		virtual void lfct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* p, ctypeRNum* xdes, ctypeRNum* udes) override;
		/*Partial derivate of the cost function with respect to states*/
		virtual void dldx(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* p, ctypeRNum* xdes, ctypeRNum* udes) override;
		/*Partial derivate of the cost function with respect to controls*/
		virtual void dldu(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* p, ctypeRNum* xdes, ctypeRNum* udes) override;

		/*Terminal cost*/
		virtual void Vfct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* p, ctypeRNum* xdes) override;
		/*Partial derivate of the terminal cost with respect to states*/
		virtual void dVdx(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* p, ctypeRNum* xdes) override;

		/*Equality constraints*/
		virtual void gfct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* p) override;
		/*Partial derivate of the equality constraints with respect to states*/
		virtual void dgdx_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* p, ctypeRNum* vec) override;
		/*Partial derivate of the equality constraints with respect to controls*/
		virtual void dgdu_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* p, ctypeRNum* vec) override;

		/*Inequality constraints*/
		virtual void hfct(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* p) override;
		/*Partial derivate of the inequality constraints with respect to states*/
		virtual void dhdx_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* p, ctypeRNum* vec) override;
		/*Partial derivate of the inequality constraints with respect to controls*/
		virtual void dhdu_vec(typeRNum* out, ctypeRNum t, ctypeRNum* x, ctypeRNum* u, ctypeRNum* p, ctypeRNum* vec) override;

	private:
		Agent* agent_;
		std::vector<int> u_index_uji_;
		std::vector<int> u_index_xji_;
		std::vector<int> u_index_vji_;
		std::vector<int> x_index_xji_;

		unsigned int Nx_;
		unsigned int Nu_;
		unsigned int Ng_;
		unsigned int Nh_;

		AgentState desired_state_;
		CouplingState couplingState_;
		MultiplierState multiplierState_;
		PenaltyState penaltyState_;
	};

}