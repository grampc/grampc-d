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

#include "dmpc/interface/dmpc_interface.hpp"


int main(int argc, char** argv)
{
	// create interface
	const auto interface = std::make_shared<dmpc::DmpcInterface>();

	// show logging
	interface->set_print_message(true);
	interface->set_print_warning(true);
	interface->set_print_error(true);

	// initialize communication interface
	interface->initialize_local_communicationInterface_as_coordinator(7777);
	
	// set optimization info
	auto optimization_info = interface->optimizationInfo();
	optimization_info.COMMON_Nhor_ = 21;
	optimization_info.COMMON_Thor_ = 5;
	optimization_info.COMMON_dt_ = 0.1;
	optimization_info.GRAMPC_MaxGradIter_ = 15;
	optimization_info.GRAMPC_MaxMultIter_ = 1;
	optimization_info.ADMM_maxIterations_ = 20;
	optimization_info.ADMM_ConvergenceTolerance_ = 0.05;

	const bool approx = true;
	optimization_info.APPROX_ApproximateConstraints_ = approx;
	optimization_info.APPROX_ApproximateCost_ = approx;
	optimization_info.APPROX_ApproximateDynamics_ = approx;

	interface->set_optimizationInfo(optimization_info);

	interface->simulate_realtime(true);

	// run distributed controller
	interface->run_DMPC();

	return 0;
}
