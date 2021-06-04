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

#include <filesystem>
#include <vector>
#include <string>

namespace unit_test {

	struct FileHandler
	{
		static const std::size_t get_number_of_text_files(const std::filesystem::path& path);

		static const std::shared_ptr<std::vector<std::string>> get_list_of_text_files(const std::filesystem::path& path);

		static void remove_text_files(const std::filesystem::path& path);

		static const bool is_txt(const std::filesystem::path path);
	};

}
