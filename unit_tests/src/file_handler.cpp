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

#include "unit_tests/file_handler.hpp"
#include "grampcd/util/logging.hpp"

#include <memory>
#include <algorithm>

const int unit_test::FileHandler::TEXT_FILES_COULD_NOT_BE_REMOVED_ = -2;

const bool unit_test::FileHandler::is_txt(const FILESYSTEM_::path path)
{
    return path.extension() == ".txt";
}

const std::size_t unit_test::FileHandler::get_number_of_text_files(const FILESYSTEM_::path& path)
{
    return std::count_if(FILESYSTEM_::directory_iterator(path), FILESYSTEM_::directory_iterator(), is_txt);
}

const std::shared_ptr<std::vector<std::string>> unit_test::FileHandler::get_list_of_text_files(const FILESYSTEM_::path& path)
{
    // prepare string
    const auto listOfTextFiles = std::make_shared<std::vector<std::string>>();

    //iterate over all files
    for (const auto& entry : FILESYSTEM_::directory_iterator(path)) 
    {
        if (is_txt(entry)) 
        {
            std::string filename = entry.path().string();
            filename.erase(filename.begin(), filename.begin() + 2);
            listOfTextFiles->push_back(filename);
        }
    }

    return listOfTextFiles;
}

void unit_test::FileHandler::remove_text_files(const FILESYSTEM_::path& path, const grampcd::LoggingPtr log_)
{
    //system("pause");
    for (const auto& dirEntry : FILESYSTEM_::recursive_directory_iterator(path)) 
    {
        if (!is_txt(dirEntry))
            continue;
        
        const bool removed_file_successfully = FILESYSTEM_::remove(dirEntry.path());
        if (!removed_file_successfully)
        {
            log_->print(grampcd::DebugType::Error) << "[FileHandler::remove_text_files]: "
                << "Fatal Error: text files were not removed sucessfully. Please check if "
                << "they are used by some other process and try again." << std::endl;
            exit(TEXT_FILES_COULD_NOT_BE_REMOVED_);
        }
    }
}

void unit_test::FileHandler::run_simulation_executable_in_working_directory(const std::string executable_name, const grampcd::LoggingPtr log_) {
    if (WINDOWS)
        const auto result = system(executable_name.c_str());
    else if (LINUX_OSX)
        const auto result = system(("./" + executable_name).c_str());
    else
    {
        log_->print(grampcd::DebugType::Error) << "[unit_test::FileHandler::run_simulation_executable_in_working_directory]: "
            << "Error: Operating system is not supported. Check file 'unit_tests/CmakeLists.txt' for more information" << std::endl;
    }
}