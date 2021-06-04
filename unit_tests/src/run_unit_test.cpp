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

#include "dmpc/util/logging.hpp"

#include "unit_tests/checksum_handler.hpp"
#include "unit_tests/file_handler.hpp"

const bool evaluate
(
	const std::string& name_of_simulation, 
	const unit_test::ChecksumHandler& checksumHandler, 
	const std::shared_ptr<dmpc::Logging>& log
) 
{
	//obtain reference values of checksums for the current simulation from checksumHandler
	const auto referenceValues = checksumHandler.get_reference_checksums(name_of_simulation);
	if (referenceValues == nullptr)
	{
		log->print(dmpc::DebugType::Error) << "[Unit_Test]: " 
			<< "Reference values not found for simulation "
			<< name_of_simulation << std::endl;
		return false;
	}

	// check number of files
	const auto desired_number_of_files = referenceValues->size();
	const auto number_of_files = unit_test::FileHandler::get_number_of_text_files(".");
	if (number_of_files != desired_number_of_files)
	{
		log->print(dmpc::DebugType::Error) << "[Unit_Test]: "
			<< "Test failed. Number of solution files should be "
			<< desired_number_of_files << " but is " << number_of_files
			<< "." << std::endl;
		return false;
	}

	// check checksum
	for (auto iter = referenceValues->begin(); iter != referenceValues->end(); ++iter)
	{
		const auto& filename = iter->first;
		const auto& desired_checksum = iter->second;
		const auto checksum = checksumHandler.generate_checksum(filename);

		if (desired_checksum != checksum)
		{
			log->print(dmpc::DebugType::Error) << "[Unit test]: "
				<< "Test failed. Checksum does not match." << std::endl;
			return false;
		}
	}

	return true;
}

int main(int argc, char* argv[])
{
	//initialize Logging and specify which MessageTypes are printed
	const auto log = std::make_shared<dmpc::Logging>();
	log->set_print_error(true);
	log->set_print_message(true);

	//setup checksumHandler and specify the output stream that should be used to print errors
	const unit_test::ChecksumHandler checksumHandler(log);

	//obtain the names of executable simulation files from checksumHandler
	const auto names_of_executables = checksumHandler.get_names_of_executables();
	
	// catch nullptr
	if (names_of_executables == nullptr)
		return -1;
	
	// flag that indicates, whether the unit test succeeded of failed
	bool whole_unit_test_suceeded = true;

	for(const auto& filename : *names_of_executables)
	{
		// delete old solution files
		unit_test::FileHandler::remove_text_files(".");
		
		//start simulation executable
		log->print(dmpc::DebugType::Message) << "[Unit test]: unit test starts for " 
			<< filename << " ..." << std::endl;

		// call .exe
		system(filename.c_str());

		// evaluate solution files
		const bool individual_test_succeeded = evaluate(filename, checksumHandler, log);

		if(individual_test_succeeded)
			log->print(dmpc::DebugType::Error) << "[Unit_Test]: Passed" << std::endl;

		whole_unit_test_suceeded &= individual_test_succeeded;

		log->print(dmpc::DebugType::Message) << std::endl;
	}

	//print test result
	if (whole_unit_test_suceeded)
		log->print(dmpc::DebugType::Message) << "[Unit test] : Unit test passed." << std::endl;
	else
		log->print(dmpc::DebugType::Error) << "[Unit test]: Unit test failed." << std::endl;

	return 0;
}
