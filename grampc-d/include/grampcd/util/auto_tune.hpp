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
#include "grampcd/info/optimization_info.hpp"

#include <mutex>

namespace grampcd
{

	class Autotune
	{
	public:
		Autotune
		(
			const LoggingPtr& log,
			const OptimizationInfo& optimization_info,
			const std::vector<AgentPtr>& agents,
			const CoordinatorPtr& coordinator 
		);

		/*Calling this functions starts the auto tuning.*/
		const OptimizationInfo auto_tune_parameters
		(
			const TuningInfo& tuning_info,
			ctypeRNum convergence_tolerance,
			const std::string& type,
			const int size_of_population,
			const int number_of_generations
		);

	private:
		/*This mutex is used to serialize the */
		std::mutex mutex_tuning_;

		/*This functions measures the required CPU time in case of MPC*/
		ctypeRNum run_MPC_tuning(const std::vector<typeRNum>& x);

		/*This functions measures the required CPU time in case of DMPC*/
		ctypeRNum run_DMPC_tuning(const std::vector<typeRNum>& x);

		/*Function pointer to the run_tuning functions.*/
		typedef  const typeRNum(Autotune::* FunctionRunTuning)(const std::vector<typeRNum>&);
		FunctionRunTuning function_run_tuning_ = nullptr;

		/*This template function for the objective is a requirement due to Galgo.*/
		template <typename T> static std::vector<T> Objective(const std::vector<T>& x);

		/*Chosen convergence tolerance*/
		typeRNum convergence_tolerance_ = 0.0;

		/*Initial optimization info that parameterizes the OCP.*/
		const OptimizationInfo optimization_info_;

		/*Pointer to agents*/
		const std::vector<AgentPtr> agents_;

		/*Pointer to coordinator*/
		const CoordinatorPtr coordinator_;

		/*Pointer to logging class*/
		const LoggingPtr log_;
	};
}