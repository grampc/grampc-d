# This file is part of GRAMPC-D - (https://github.com/grampc-d/grampc-d.git)
#
# GRAMPC-D -- A software framework for distributed model predictive control (DMPC)
# based on the alternating direction method of multipliers (ADMM).
#
# Copyright 2020 by Daniel Burk, Andreas Voelz, Knut Graichen
# All rights reserved.
#
# GRAMPC-D is distributed under the BSD-3-Clause license, see LICENSE.txt

include_directories(
    ${CMAKE_SOURCE_DIR}/libs/sha/src/sha
)

include(ExternalProject)
find_package(Git REQUIRED)

# include SHA1 as external project
ExternalProject_Add(sha
    GIT_REPOSITORY https://github.com/983/SHA1.git
    GIT_TAG 0ec907360c70b7eac33ec87927d227c602bb1265
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
  	    BUILD_COMMAND ""
  	    INSTALL_COMMAND ""
	    LOG_DOWNLOAD ON
        PREFIX  ${CMAKE_SOURCE_DIR}/libs/sha
)

add_library(unit_test_handler
    src/checksum_handler.cpp
    src/file_handler.cpp
)

target_include_directories(unit_test_handler PUBLIC ${CMAKE_SOURCE_DIR}/grampcd/include)

target_link_libraries(unit_test_handler
    grampcd
)

add_dependencies(unit_test_handler sha)

add_executable(run_unit_test
    src/run_unit_test.cpp
)

target_include_directories(run_unit_test PUBLIC ${CMAKE_SOURCE_DIR}/grampcd/include)

target_link_libraries(run_unit_test
    grampcd
    unit_test_handler
)

add_executable(setup_unit_test
    src/setup_unit_test.cpp
)

target_include_directories(setup_unit_test PUBLIC ${CMAKE_SOURCE_DIR}/grampcd/include)

target_link_libraries( setup_unit_test
    unit_test_handler
)