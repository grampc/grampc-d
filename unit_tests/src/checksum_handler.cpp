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

#include "unit_tests/file_handler.hpp"
#include "unit_tests/checksum_handler.hpp"

#include "sha1.hpp"

#include <fstream>
#include <sstream>
#include <cstdlib>

namespace unit_test
{
	ChecksumHandler::ChecksumHandler(const std::shared_ptr<grampcd::Logging>& log) : log_(log) 
	{
	};

	void ChecksumHandler::strip(std::string& s) const
	{
		s.erase(std::remove(s.begin(), s.end(), ' '), s.end());
		s.erase(std::remove(s.begin(), s.end(), '\t'), s.end());
	}

	const std::vector<std::string> ChecksumHandler::split(const std::string& inputString, const char delimiter = ';') const
	{
		std::vector<std::string> tokens;
		std::string token;
		std::istringstream tokenStream(inputString);

		while (std::getline(tokenStream, token, delimiter))
			tokens.push_back(token);

		return tokens;
	}

	const std::shared_ptr<std::ifstream> ChecksumHandler::open_file_including_checksums() const
	{
		const std::shared_ptr<std::ifstream> file_including_checksums = std::make_shared<std::ifstream>(std::ifstream());

		// open file with checksums
		file_including_checksums->open(PATH_REFERENCE_CHECKSUMS_);

		// check if file is successfully opened
		if (!file_including_checksums->is_open())
		{
			log_->print(grampcd::DebugType::Error) << "[ChecksumHandler::open_file_including_checksums]: "
				<< "Error: file 'referenceChecksums.csv' not found at " << PATH_REFERENCE_CHECKSUMS_ << std::endl;
			return nullptr;
		}
		return file_including_checksums;
	}

	const std::shared_ptr<std::string> ChecksumHandler::read_file_to_string(const std::string& filename) const
	{
		std::ifstream file_to_read_to_string(SIM_DIR_ + filename);
		if (!file_to_read_to_string.is_open())
		{
			log_->print(grampcd::DebugType::Error) << "[ChecksumHandler::read_file_to_string]: "
				<< "Error: file '" << filename << "' cannot be opened at " << SIM_DIR_ + filename << std::endl;
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

	const std::string ChecksumHandler::generate_checksum(const std::string& filename) const
	{

		auto file_content = read_file_to_string(filename);
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

	const std::shared_ptr<std::vector<std::string>> ChecksumHandler::get_names_of_executables() const
	{
		const auto file_including_checksums = open_file_including_checksums();

		// check if file was found
		if (file_including_checksums == nullptr)
		{
			log_->print(grampcd::DebugType::Error) << "File including the reference checksums "
				<< "could not be opened." << std::endl;
			return nullptr;
		}

		// prepare string with names of executables
		auto simExecutableNames = std::make_shared<std::vector<std::string>>();

		// read names of executables from file
		std::string line;
		while (std::getline(*file_including_checksums, line))
		{
			strip(line);
			if (line != "" && split(line, ';').size() == 1 && line[0] != '#')
				simExecutableNames->push_back(line);
		}

		// close file
		file_including_checksums->close();

		return simExecutableNames;
	}

	const std::shared_ptr<std::map<std::string, std::string>> ChecksumHandler::get_reference_checksums(const std::string& name) const
	{
		//create shared pointer to map of string pairs to store checksums
		auto referenceChecksums = std::make_shared<std::map<std::string, std::string>>();
		bool sucessful_read = false;

		// get file including checksums
		const auto file_including_checksums = open_file_including_checksums();

		// if file was not successfully opened, don't try to close file_including_checksums 
		// as is the nullptr, so trying to close the file here would result in an read access 
		// violation
		if (file_including_checksums == nullptr)
			return nullptr;

		// read names of executables from file
		std::string line;
		while (std::getline(*file_including_checksums, line))
		{
			strip(line);
			if (line == name)
			{
				sucessful_read = true;
				break;
			}
		}

		// check if names are successfully read
		if (!sucessful_read)
		{
			log_->print(grampcd::DebugType::Error) << "[ChecksumHandler::get_reference_checksums]: File for executable '" + name + "' not found";
			file_including_checksums->close();
			return nullptr;
		}

		while (std::getline(*file_including_checksums, line))
		{
			strip(line);
			if (line != "")
			{
				std::vector<std::string> elements = split(line);
				if (elements.size() < 2)
					break;
				referenceChecksums->insert(std::pair<std::string, std::string>(elements[0], elements[1]));
			}

		}

		// close file
		file_including_checksums->close();
		return referenceChecksums;
	}

	void ChecksumHandler::update_file_including_checksums(const unsigned int startDel, const unsigned int endDel, const std::string& textToInsert) const
	{
		// open file stream for temporary file "referenceChecksums.csv_temp"
		std::ofstream file_stream;
		file_stream.open(PATH_REFERENCE_CHECKSUMS_ + "_temp");

		// open file including checksums
		const auto referenceChecksums_file = open_file_including_checksums();

		// check if file was successfully opened
		if (referenceChecksums_file == nullptr)
			return;

		// iterate over original referenceChecksums.csv file and replace the lines from line startDel to line endDel with the string textToInsert
		std::string line;
		unsigned int row = 0;
		while (std::getline(*referenceChecksums_file, line))
		{
			++row;
			// copy every line from referenceChecksums.csv into referenceChecksums.csv_temp file except the lines between startDel and endDel 
			if (row < startDel || row > endDel)
				file_stream << line << std::endl;
			//insert the string textToString at line startDel
			else if (row == startDel)
				file_stream << textToInsert;
		}

		// close file
		referenceChecksums_file->close();

		// close stream
		file_stream.close();

		//replacing original "referenceChecksums.csv" file with referenceChecksums.csv_temp file
		const auto result_remove = std::remove(PATH_REFERENCE_CHECKSUMS_.c_str());
		const auto result_rename = std::rename((PATH_REFERENCE_CHECKSUMS_ + "_temp").c_str(), PATH_REFERENCE_CHECKSUMS_.c_str());
	}

	void ChecksumHandler::add_unit_test(const std::string& nameOfSimExecutable) const
	{
		//check if there is already an entry for that simulation in referenceChecksums.csv
		const auto linesOfSimEntry = evaluate_row_in_checksum_file(nameOfSimExecutable);

		//check if there went something wrong during reading the referenceChecksums.csv file
		if (linesOfSimEntry->first == referenceChecksum_file_could_not_be_opened)
			return;

		//if there is already an entry for that simulation, ask user what to do
		if (linesOfSimEntry->first > 0) {
			log_->print(grampcd::DebugType::Base) << "There is already a entry for \"" << nameOfSimExecutable
				<< "\"\nWhat do you want to do?\nc = cancel (will stop the whole executable)\ns"
				<< "= skip (will skip this operation, but continue the program)\nu = update entry"
				<< " (will run the simulation and update the entry files with the new data)" << std::endl;
			char input_option;

			while (std::cin >> input_option, input_option != 'c' && input_option != 's' && input_option != 'u')
				log_->print(grampcd::DebugType::Base) << "unknown option. Please choose c (cancel), s (skip) or u (update)" << std::endl;

			if (input_option == 'c')
				std::exit(EXIT_FAILURE);
			else if (input_option == 's')
				return;
			//if input_option is u, resume with program
		}

		log_->print(grampcd::DebugType::Base) << "running simulation " << nameOfSimExecutable << std::endl;

		//remove all .txt files before running simulation to make sure that there are no old simulation files  
		FileHandler::remove_text_files(".");

		//run simulation
		const std::string pathToSimExecutable = nameOfSimExecutable + ".exe";
		system(pathToSimExecutable.c_str());

		//get a list of names of the simulation files
		const auto simFileList = FileHandler::get_list_of_text_files(".");

		std::string newSimEntry = nameOfSimExecutable + "\n";
		//iterate over each simFile, create checksum and add it to newSimEntry
		for (const auto& simFileName : *simFileList)
		{
			const std::string checksum = generate_checksum(simFileName);
			newSimEntry.append(simFileName + ";" + checksum + "\n");
		}

		//check if checksums have to be added or updated
		if (linesOfSimEntry->first > 0)
		{
			//referenceChecksums.csv needs to be updated
			update_file_including_checksums(linesOfSimEntry->first, linesOfSimEntry->second, newSimEntry);
		}
		else
		{
			//open output stream to referenceChecksums.csv file in appendix mode
			std::ofstream referenceChecksums_file;
			referenceChecksums_file.open(PATH_REFERENCE_CHECKSUMS_, std::ios::app);

			//check if opening file was successful
			if (!referenceChecksums_file.is_open())
			{
				log_->print(grampcd::DebugType::Error) << "[ChecksumHandler::addUnitTest]: "
					<< "ReferenceChecksums.csv could not be opened. Path: " << PATH_REFERENCE_CHECKSUMS_ << std::endl;
				return;
			}

			//add entry with the name of the current simulation and close file
			referenceChecksums_file << "\n\n" << newSimEntry;
			referenceChecksums_file.close();
		}

		log_->print(grampcd::DebugType::Base) << std::endl;
	}

	//returns the number of the line containing the string <entry_name>
	const std::shared_ptr<std::pair<int, int>> ChecksumHandler::evaluate_row_in_checksum_file(const std::string& entry_name) const
	{
		//initialize line pair with case that no match was found 
		const auto linesOfSimEntry = std::make_shared<std::pair<int, int>>(no_match_for_requested_simulation_found, 0);

		// open file including checksums
		const auto referenceChecksums_file = open_file_including_checksums();

		// check if file is successfully opened
		if (referenceChecksums_file == nullptr)
		{
			// save error cause in first entry
			linesOfSimEntry->first = referenceChecksum_file_could_not_be_opened;
			return linesOfSimEntry;
		}

		std::string line;
		unsigned int row = 0;
		bool matchFound = false;
		while (std::getline(*referenceChecksums_file, line))
		{
			++row;
			strip(line);
			if (line == entry_name)
			{
				linesOfSimEntry->first = row;
				matchFound = true;
			}
			else if (matchFound && line == "")
			{
				linesOfSimEntry->second = row - 1;
				return linesOfSimEntry;
			}
		}

		// close file
		referenceChecksums_file->close();

		// check if match is found
		if (matchFound)
			linesOfSimEntry->second = row;

		return linesOfSimEntry;
	}
}

