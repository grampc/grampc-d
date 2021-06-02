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
