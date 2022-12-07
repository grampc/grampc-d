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

#include FILESYSTEM_INCLUDE_FILE_
#include <vector>
#include <string>

#include "grampcd/util/types.hpp"

namespace grampcd {
	DMPC_CLASS_FORWARD(Logging)
}

namespace unit_test {

	struct FileHandler
	{
		static const int TEXT_FILES_COULD_NOT_BE_REMOVED_;

		static const std::size_t get_number_of_text_files(const FILESYSTEM_::path& path);

		static const std::shared_ptr<std::vector<std::string>> get_list_of_text_files(const FILESYSTEM_::path& path);

		static void remove_text_files(const FILESYSTEM_::path& path, const grampcd::LoggingPtr log);

		static const bool is_txt(const FILESYSTEM_::path path);

		static void run_simulation_executable_in_working_directory(const std::string executable_name, const grampcd::LoggingPtr log);
	};

}
