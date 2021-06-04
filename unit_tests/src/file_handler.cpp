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

#include <memory>

const bool unit_test::FileHandler::is_txt(const std::filesystem::path path)
{
    return path.extension() == ".txt";
}

const std::size_t unit_test::FileHandler::get_number_of_text_files(const std::filesystem::path& path)
{
    return std::count_if(std::filesystem::directory_iterator(path), std::filesystem::directory_iterator(), is_txt);
}

const std::shared_ptr<std::vector<std::string>> unit_test::FileHandler::get_list_of_text_files(const std::filesystem::path& path)
{
    // prepare string
    const auto listOfTextFiles = std::make_shared<std::vector<std::string>>();

    //iterate over all files
    for (const auto& entry : std::filesystem::directory_iterator(path)) 
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

void unit_test::FileHandler::remove_text_files(const std::filesystem::path& path)
{
    //system("pause");
    for (const auto& dirEntry : std::filesystem::recursive_directory_iterator(path)) 
    {
        if (is_txt(dirEntry)) 
            std::remove(dirEntry.path().string().c_str());
    }
}