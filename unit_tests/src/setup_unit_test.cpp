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

#include "grampcd/util/logging.hpp"

#include "unit_tests/checksum_handler.hpp"
#include "unit_tests/file_handler.hpp"

int main(int argc, char* argv[]) 
{
	const auto log = std::make_shared<grampcd::Logging>();
	unit_test::ChecksumHandler checksumHandler(log);

	checksumHandler.add_unit_test("coupled_watertanks");
	checksumHandler.add_unit_test("evaluate_neighborApproximation");
	checksumHandler.add_unit_test("high_scaled_system");
	checksumHandler.add_unit_test("plug-and-play");
}