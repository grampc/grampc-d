# This file is part of GRAMPC-D - (https://github.com/grampc-d/grampc-d.git)
#
# GRAMPC-D -- A software framework for distributed model predictive control (DMPC)
# 
#
# Copyright 2023 by Daniel Burk, Maximilian Pierer von Esch, Andreas Voelz, Knut Graichen
# All rights reserved.
#
# GRAMPC-D is distributed under the BSD-3-Clause license, see LICENSE.txt

cmake_minimum_required(VERSION 3.0)
project(grampcd VERSION 1.0.0)

set(PYTHON_AVAILABLE FALSE)
set(UNIT_TEST FALSE)
set(AUTO_TUNE FALSE)

# add support for C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# include missing compiler flags
if(MSVC)
    SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -D_WIN32_WINNT=0x0A00 -bigobj")
    SET(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ")
else()
    SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wno-psabi -fpic")
    SET(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -fpic")
endif()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# link_directories should treat paths as relative to the source directory
cmake_policy(SET CMP0015 NEW)

#add libs
add_subdirectory(libs)

# add grampc-d
add_subdirectory(grampc-d)

# add simulation examples
add_subdirectory(examples)

# add model descriptions
add_subdirectory(model_description)

# add unit tests
if(UNIT_TEST)  
    add_subdirectory(unit_tests)
endif()