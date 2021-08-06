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

const bool evaluate
(
	const std::shared_ptr<unit_test::SimulationResult> simulation_result_reference,
	const std::shared_ptr<grampcd::Logging>& log
) 
{

	// check number of files
	const auto desired_number_of_files = simulation_result_reference->simulationFile_References_.size();
	const auto number_of_files = unit_test::FileHandler::get_number_of_text_files(".");
	if (number_of_files != desired_number_of_files)
	{
		log->print(grampcd::DebugType::Error) << "[Unit_Test]: "
			<< "Test failed. Number of solution files should be "
			<< desired_number_of_files << " but is " << number_of_files
			<< "." << std::endl;
		return false;
	}

	// check checksum
	for (const auto simulation_file_reference : simulation_result_reference->simulationFile_References_)
	{
		const auto& filename = simulation_file_reference->filename_;
		const auto& desired_checksum = simulation_file_reference->hash_code_;
		const auto checksum = unit_test::ChecksumHandler::generate_checksum(filename, log);

		if (desired_checksum != checksum)
		{
			log->print(grampcd::DebugType::Error) << "[Unit test]: "
				<< "Test failed. Checksum does not match." << std::endl;
			return false;
		}
	}

	return true;
}

int main(int argc, char* argv[])
{
	//initialize Logging and specify which MessageTypes are printed
	const auto log = std::make_shared<grampcd::Logging>();
	log->set_print_error(true);
	log->set_print_message(true);

	//obtain the names of executable simulation files from checksumHandler
	const auto results_of_all_unittests = unit_test::ChecksumHandler::get_results_of_all_unittests(log);
	
	// catch case that referenceChecksums.txt could not be read
	if (results_of_all_unittests == nullptr)
		return unit_test::ChecksumHandler::REFERENCE_CHECKSUM_FILE_NOT_READABLE_;
	
	// flag that indicates, whether the unit test succeeded of failed
	bool whole_unit_test_suceeded = true;

	for(const auto& simulation_result_reference : results_of_all_unittests->simulationResults)
	{
		if (!simulation_result_reference->runtest_)
			continue;
	
		// delete old solution files
		unit_test::FileHandler::remove_text_files(".", log);
		
		//start simulation executable
		log->print(grampcd::DebugType::Message) << "[Unit test]: unit test starts for " 
			<< simulation_result_reference->executable_name_ << " ..." << std::endl;

		// call executable
		unit_test::FileHandler::run_simulation_executable_in_working_directory(simulation_result_reference->executable_name_, log);

		// evaluate solution files
		const bool individual_test_succeeded = evaluate(simulation_result_reference, log);

		if(individual_test_succeeded)
			log->print(grampcd::DebugType::Error) << "[Unit_Test]: Passed" << std::endl;

		whole_unit_test_suceeded &= individual_test_succeeded;

		log->print(grampcd::DebugType::Message) << std::endl;
	}

	//print test result
	if (whole_unit_test_suceeded)
		log->print(grampcd::DebugType::Message) << "[Unit test] : Unit test passed." << std::endl;
	else
		log->print(grampcd::DebugType::Error) << "[Unit test]: Unit test failed." << std::endl;

	return 0;
}
