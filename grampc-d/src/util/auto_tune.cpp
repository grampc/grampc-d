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

#include "grampcd/optim/solver_central.hpp"

#include "grampcd/agent/agent.hpp"

#include "grampcd/coord/coordinator.hpp"

#include "grampcd/util/logging.hpp"
#include "grampcd/util/auto_tune.hpp"

#include "grampcd/info/tuning_info.hpp"

#ifdef AUTO_TUNE
#include "galgo.hpp"
#endif

namespace grampcd
{
	/*This global variables allows to access member variables from static functions.*/
	Autotune* global_autotune_;

	Autotune::Autotune
	(
		const LoggingPtr& log,
		const OptimizationInfo& optimization_info,
		const std::vector<AgentPtr>& agents,
		const CoordinatorPtr& coordinator
	) :
		log_(log),
		optimization_info_(optimization_info),
		agents_(agents),
		coordinator_(coordinator)
	{
		global_autotune_ = this;
	}

	template<typename T> std::vector<T> Autotune::Objective(const std::vector<T>& x)
	{
		// This template function simply calls the corresponding run_tuning function.
		return { -((*global_autotune_).*(global_autotune_->function_run_tuning_))(x) };
	}

	const OptimizationInfo Autotune::auto_tune_parameters
	(
		const TuningInfo& tuning_info,
		ctypeRNum convergence_tolerance,
		const std::string& type,
		const int size_of_population,
		const int number_of_generations
	)
	{
#ifndef AUTO_TUNE
		log_->print(DebugType::Error) << "[Autotune::auto_tune_parameters]: "
			<< "The auto-tuning functionality is not available at the moment. "
			<< "Please set AUTO_TUNE in the top CMakeLists.txt to turn it on." << std::endl;

		return optimization_info_;
#else
		convergence_tolerance_ = convergence_tolerance;

		if (type == "MPC")
		{
			// set corresponding run_tuning function
			function_run_tuning_ = &Autotune::run_MPC_tuning;

			// initialize Galgo
			galgo::Parameter<double, 6> par0({ tuning_info.COMMON_Nhor_[0], tuning_info.COMMON_Nhor_[1] });
			galgo::Parameter<double, 6> par1({ tuning_info.GRAMPC_MaxGradIter_[0], tuning_info.GRAMPC_MaxGradIter_[1] });
			galgo::Parameter<double, 3> par2({ tuning_info.GRAMPC_MaxMultIter_[0], tuning_info.GRAMPC_MaxMultIter_[1] });

			galgo::GeneticAlgorithm<double> ga(Objective, size_of_population, number_of_generations, true, par0, par1, par2);

			// Print each generation
			ga.genstep = 1;
			ga.precision = 0;

			// run Galgo
			ga.run();

			// get tuning result
			const std::vector<double> tuning_result = ga.result()->getParam();

			// Set tuned parameters in initial optimization info.
			OptimizationInfo optimization_info = optimization_info_;
			optimization_info.COMMON_Nhor_ = static_cast<unsigned int>(std::round(tuning_result[0]));
			optimization_info.GRAMPC_MaxGradIter_ = static_cast<unsigned int>(std::round(tuning_result[1]));
			optimization_info.GRAMPC_MaxMultIter_ = static_cast<unsigned int>(std::round(tuning_result[2]));

			return optimization_info;
		}
		else if (type == "DMPC")
		{
			// set corresponding run_tuning function
			function_run_tuning_ = &Autotune::run_DMPC_tuning;

			// initialize Galgo
			galgo::Parameter<double, 6> par0({ tuning_info.COMMON_Nhor_[0], tuning_info.COMMON_Nhor_[1] });
			galgo::Parameter<double, 6> par1({ tuning_info.GRAMPC_MaxGradIter_[0], tuning_info.GRAMPC_MaxGradIter_[1] });
			galgo::Parameter<double, 3> par2({ tuning_info.GRAMPC_MaxMultIter_[0], tuning_info.GRAMPC_MaxMultIter_[1] });
			galgo::Parameter<double, 6> par3({ tuning_info.ADMM_maxIterations_[0], tuning_info.ADMM_maxIterations_[1] });
			galgo::Parameter<double, 3> par4({ tuning_info.ADMM_innerIterations_[0], tuning_info.ADMM_innerIterations_[1] });

			galgo::GeneticAlgorithm<double> ga(Objective, size_of_population, number_of_generations,
				true, par0, par1, par2, par3, par4);

			// Print each generation
			ga.genstep = 1;
			ga.precision = 0;

			// run Galgo
			ga.run();

			// get tuning result
			std::vector<double> tuning_result = ga.result()->getParam();

			// Galgo does not find the minimum number of ADMM iterations, which is fixed by this small algorithm.
			while (run_DMPC_tuning(tuning_result) < 5e5)
				--tuning_result[3];
			++tuning_result[3];

			// Set tuned parameters in initial optimization info.
			OptimizationInfo optimization_info = optimization_info_;
			optimization_info.COMMON_Nhor_ = static_cast<unsigned int>(std::round(tuning_result[0]));
			optimization_info.GRAMPC_MaxGradIter_ = static_cast<unsigned int>(std::round(tuning_result[1]));
			optimization_info.GRAMPC_MaxMultIter_ = static_cast<unsigned int>(std::round(tuning_result[2]));
			optimization_info.ADMM_maxIterations_ = static_cast<unsigned int>(std::round(tuning_result[3]));
			optimization_info.ADMM_innerIterations_ = static_cast<unsigned int>(std::round(tuning_result[4]));

			return optimization_info;
		}
		else
		{
			log_->print(DebugType::Error) << "[DmpcInterface::auto_tune_parameters]: "
				<< "Unknown type " << type << ". Allowed types are MPC and DMPC." << std::endl;
			return optimization_info_;
		}
#endif
	}

	ctypeRNum Autotune::run_MPC_tuning(const std::vector<typeRNum>& x)
	{
		// Serialize the tuning process
		std::lock_guard<std::mutex> guard(mutex_tuning_);

		// Generate local copy of current candidate
		OptimizationInfo optimization_info = optimization_info_;
		optimization_info.GRAMPC_ConvergenceCheck_ = "on";
		optimization_info.GRAMPC_ConvergenceGradientRelTol_ = convergence_tolerance_;
		optimization_info.COMMON_Nhor_ = static_cast<int>(std::round(x[0]));
		optimization_info.GRAMPC_MaxGradIter_ = static_cast<unsigned int>(std::round(x[1]));
		optimization_info.GRAMPC_MaxMultIter_ = static_cast<unsigned int>(std::round(x[2]));

		// Initialize central solver
		SolverCentral solver(agents_, optimization_info);

		// Initialize agents
		for (const auto& agent : agents_)
			agent->initialize(optimization_info);

		// Run optimization and measure required time
		const auto tstart = std::chrono::high_resolution_clock::now();
		solver.solve();
		const auto elapsed = std::chrono::high_resolution_clock::now() - tstart;
		const long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
		ctypeRNum measured_time = static_cast<typeRNum>(microseconds);

		if (solver.is_converged())
			return measured_time;
		else
			return 1e6;
	}

	ctypeRNum Autotune::run_DMPC_tuning(const std::vector<typeRNum>& x)
	{
		// Serialize the tuning process
		std::lock_guard<std::mutex> guard(mutex_tuning_);

		// Generate local copy of current candidate
		OptimizationInfo optimization_info = optimization_info_;
		optimization_info.ADMM_ConvergenceTolerance_ = convergence_tolerance_;
		optimization_info.COMMON_Nhor_ = static_cast<unsigned int>(std::round(x[0]));
		optimization_info.GRAMPC_MaxGradIter_ = static_cast<unsigned int>(std::round(x[1]));
		optimization_info.GRAMPC_MaxMultIter_ = static_cast<unsigned int>(std::round(x[2]));
		optimization_info.ADMM_maxIterations_ = static_cast<unsigned int>(std::round(x[3]));
		optimization_info.ADMM_innerIterations_ = static_cast<unsigned int>(std::round(x[4]));

		// Initialize local solver
		coordinator_->initialize_ADMM(optimization_info);

		// Run optimization and measure required time
		const auto tstart = std::chrono::high_resolution_clock::now();
		const bool converged = coordinator_->solve_ADMM(optimization_info.ADMM_maxIterations_, optimization_info.ADMM_innerIterations_);
		const auto elapsed = std::chrono::high_resolution_clock::now() - tstart;
		const long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
		typeRNum measured_time = static_cast<typeRNum>(microseconds) / static_cast<typeRNum>(agents_.size());

		// This adaption of the time tells the algorithm to increase to computational effort
		// if the algorithm does not converge
		if(!converged)
			measured_time = 1e6 + 1e6 * std::exp(-1e-5 * measured_time);

		return measured_time;
	}
}