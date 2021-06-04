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

#pragma once

#include <vector>
#include <map>
#include <string>
#include <memory>
#include "dmpc/util/types.hpp"

#include "dmpc/util/logging.hpp"

namespace unit_test 
{

	class ChecksumHandler 
	{
	public:
		enum error_codes 
		{
			referenceChecksum_file_could_not_be_opened = -1,
			no_match_for_requested_simulation_found = -2
		};

		// constructor initializes Logging
		ChecksumHandler(const std::shared_ptr<dmpc::Logging>&);

		// returns a map of reference data; the keys are the solution file names of the given simulation, 
		// the values are the reference checksums of the according simulation files
		const std::shared_ptr<std::map<std::string, std::string>> get_reference_checksums(const std::string& name) const;

		// returns a shared pointer to a string holding the content of the given file
		const std::shared_ptr<std::string> read_file_to_string(const std::string& filename) const;

		// creates a checksum from the given file and returns it
		const std::string generate_checksum(const std::string& filename) const;

		// searches the referenceChecksums.csv file for the given entry_name and returns a pair int values, 
		// where the first value is -2, if the referenceChecksums.csv was not opened successfully and -1, if the entry name was not found
		// otherwise the first value corresponds to the first line of the given entry and the second value corresponds to the last line 
		// that is part of the given entry
		const std::shared_ptr<std::pair<int, int>> evaluate_row_in_checksum_file(const std::string& entry_name) const;

		// returns the names of every simEntry that is not commented out
		const std::shared_ptr<std::vector<std::string>> get_names_of_executables() const;

		// adds a unitTest for the simulation which executable is called <executableName>.exe
		void add_unit_test(const std::string& executableName) const;

	private:
		const std::shared_ptr<dmpc::Logging> log_;

		//path relative to bin
		const std::string PATH_REFERENCE_CHECKSUMS_ = "../unit_tests/referenceChecksums.csv";
		const std::string SIM_DIR_ = "./";

		void strip(std::string& s) const;
		const std::shared_ptr<std::ifstream> open_file_including_checksums() const;
		void update_file_including_checksums(const unsigned int startDel, const unsigned int endDel, const std::string& textToInsert) const;
		const std::vector<std::string> split(const std::string& inputString, const char delimiter) const;
	};

}
