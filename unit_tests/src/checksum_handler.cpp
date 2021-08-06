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

#include "unit_tests/checksum_handler.hpp"
#include "grampcd/util/logging.hpp"

#include <sha1.hpp>
#include <fstream>

const std::string unit_test::ChecksumHandler::PATH_REFERENCE_CHECKSUMS_ = "../unit_tests/referenceChecksums.txt";
const int unit_test::ChecksumHandler::REFERENCE_CHECKSUM_FILE_NOT_READABLE_ = -1;

const std::shared_ptr<std::string> unit_test::ChecksumHandler::read_file_to_string(const std::string& filename, grampcd::LoggingPtr log_)
{
	std::ifstream file_to_read_to_string(filename);
	if (!file_to_read_to_string.is_open())
	{
		log_->print(grampcd::DebugType::Error) << "[ChecksumHandler::read_file_to_string]: "
			<< "Error: file '" << filename << "' cannot be opened." << std::endl;
		return nullptr;
	}

	// prepare buffer
	const auto buffer = std::make_shared<std::string>("");

	std::string line;
	while (std::getline(file_to_read_to_string, line))
		*buffer += line + "\n";

	// close file
	file_to_read_to_string.close();

	return buffer;
}

const std::string unit_test::ChecksumHandler::generate_checksum(const std::string& filename, grampcd::LoggingPtr log_)
{

	auto file_content = read_file_to_string(filename, log_);
	if (file_content == nullptr)
	{
		log_->print(grampcd::DebugType::Error) << "[ChecksumHandler::generate_checksum]: "
			<< "Error: cannot generate checksum for file '" << filename << "' because reading the file failed." << std::endl;
		return "";
	}

	// generate SHA1 checksum
	char checksum[SHA1_HEX_SIZE];
	sha1(file_content->c_str()).finalize().print_hex(checksum);

	return checksum;
}

const std::shared_ptr<unit_test::Results_of_UnitTests> unit_test::ChecksumHandler::get_results_of_all_unittests(const grampcd::LoggingPtr log_)
{
	std::ifstream input_file_stream;
	input_file_stream.open(PATH_REFERENCE_CHECKSUMS_, std::ios::binary);
	if (!input_file_stream.is_open())
	{
		log_->print(grampcd::DebugType::Error) << "[ChecksumHandler::get_results_of_all_unittests]: "
			<< "referenceChecksums.txt cannot be opened." << std::endl;
		return nullptr;
	}


	std::shared_ptr<Results_of_UnitTests> all_results;

	if (input_file_stream.peek() != std::ifstream::traits_type::eof())
	{
		// deserialize the reference data for the unit test results
		cereal::BinaryInputArchive iarchive(input_file_stream);
		iarchive(all_results);
	}
	else
	{
		//return a shared_pointer to an empty vector of SimulationResults
		all_results = std::make_shared<Results_of_UnitTests>();
	}

	return all_results;
}

void unit_test::ChecksumHandler::serialize_reference_data(std::shared_ptr<unit_test::Results_of_UnitTests> results_of_all_unittests, grampcd::LoggingPtr log_)
{
	//serialize data
	std::ofstream output_file_stream;
	output_file_stream.open(PATH_REFERENCE_CHECKSUMS_, std::ios::binary);
	if (output_file_stream.is_open())
	{
		cereal::BinaryOutputArchive oarchive(output_file_stream);
		oarchive(results_of_all_unittests);
		output_file_stream.close();
	}
	else
	{
		log_->print(grampcd::DebugType::Error) << "[ChecksumHandler::serialize_reference_data]: "
			<< "refereceChecksums.txt cannot be opened." << std::endl;
		exit(REFERENCE_CHECKSUM_FILE_NOT_READABLE_);
	}
}