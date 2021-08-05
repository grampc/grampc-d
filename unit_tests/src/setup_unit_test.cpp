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

std::shared_ptr<unit_test::SimulationResult> generate_SimulationResult_from_simulation_run(const std::string& executableName, const grampcd::LoggingPtr log_)
{
	// delete old solution files
	unit_test::FileHandler::remove_text_files(".", log_);

	// call executable
	unit_test::FileHandler::run_simulation_executable_in_working_directory(executableName, log_);

	// get all simulation file names
	const auto simulation_file_names = unit_test::FileHandler::get_list_of_text_files(".");

	//create vector for all simulation files with corresponding checksums
	std::vector<std::shared_ptr<unit_test::SimulationFile_Reference>> simulation_data_vector(simulation_file_names->size());

	//iterate over all simulation files
	for (int i = 0; i < simulation_file_names->size(); ++i)
	{
		// obtain filename
		std::string simulation_file_name = (*simulation_file_names)[i];

		//create a simulationFile_reference and store the data of the current simulation file
		unit_test::SimulationFile_Reference simulation_file_reference;

		simulation_file_reference.filename_ = simulation_file_name;
		simulation_file_reference.hash_code_ = unit_test::ChecksumHandler::generate_checksum(simulation_file_name, log_);
		simulation_data_vector[i] = std::make_shared< unit_test::SimulationFile_Reference>(simulation_file_reference);
	}

	//put together the simulation Result
	const auto simulationResult = std::make_shared<unit_test::SimulationResult>();
	simulationResult->executable_name_ = executableName;
	simulationResult->runtest_ = true;
	simulationResult->simulationFile_References_ = simulation_data_vector;
	return simulationResult;
}

void add_unit_test(const std::string& executable_name_to_add, const grampcd::LoggingPtr log_)
{
	// get results of all unittests
	const auto results_of_all_unittests = unit_test::ChecksumHandler::get_results_of_all_unittests(log_);

	//check if file was successfully read
	if (results_of_all_unittests == nullptr)
		exit(unit_test::ChecksumHandler::REFERENCE_CHECKSUM_FILE_NOT_READABLE_);

	//iterate over reference data to check if there are already entrys for requested simulation
	for (int i = 0; i < results_of_all_unittests->simulationResults.size(); ++i)
	{
		const bool name_already_exists = results_of_all_unittests->simulationResults[i]->executable_name_ == executable_name_to_add;
		if (!name_already_exists ) 
			continue;

		log_->print(grampcd::DebugType::Base) << "There is already a entry for \"" << executable_name_to_add
			<< "\"\nWhat do you want to do?\nc = cancel (will stop the whole executable)\n"
			<< "s = skip (will skip this operation, but continue the program)\nu = update entry"
			<< " (will run the simulation and update the entry files with the new data)" << std::endl;
		char input_option;

		while (std::cin >> input_option, input_option != 'c' && input_option != 's' && input_option != 'u')
			log_->print(grampcd::DebugType::Base) << "unknown option. Please choose c (cancel), s (skip) or u (update)" << std::endl;

		if (input_option == 'c')
			std::exit(EXIT_FAILURE);
		else if (input_option == 's')
			return;
		else /* input_option == 'u' */
		{
			const auto simulationResult = generate_SimulationResult_from_simulation_run(executable_name_to_add, log_);
			results_of_all_unittests->simulationResults[i] = simulationResult;
			unit_test::ChecksumHandler::serialize_reference_data(results_of_all_unittests, log_);
			return;
		}
	}

	std::shared_ptr<unit_test::SimulationResult> simulationResult = generate_SimulationResult_from_simulation_run(executable_name_to_add, log_);
	results_of_all_unittests->simulationResults.push_back(simulationResult);
	unit_test::ChecksumHandler::serialize_reference_data(results_of_all_unittests, log_);

}

int main(int argc, char* argv[]) 
{
	const auto log = std::make_shared<grampcd::Logging>();
	log->set_print_error(true);

	add_unit_test("coupled_watertanks", log);
	add_unit_test("evaluate_neighborApproximation", log);
	add_unit_test("high_scaled_system", log);
	add_unit_test("coupled_cost_functions", log);
}