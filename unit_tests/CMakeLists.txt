# This file is part of GRAMPC-D - (https://github.com/grampc-d/grampc-d.git)
#
# GRAMPC-D -- A software framework for distributed model predictive control (DMPC)
# based on the alternating direction method of multipliers (ADMM).
#
# Copyright 2020 by Daniel Burk, Andreas Voelz, Knut Graichen
# All rights reserved.
#
# GRAMPC-D is distributed under the BSD-3-Clause license, see LICENSE.txt

add_library(unit_test_handler
    src/checksum_handler.cpp
    src/file_handler.cpp
)

target_include_directories(unit_test_handler PUBLIC
    include
    ${SHA_INCLUDE_DIR}
)

target_link_libraries(unit_test_handler
    grampcd
)

add_dependencies(unit_test_handler 
    sha
    cereal
)

add_executable(run_unit_test
    src/run_unit_test.cpp
)

target_link_libraries(run_unit_test
    unit_test_handler
)

add_executable(setup_unit_test
    src/setup_unit_test.cpp
)

target_link_libraries( setup_unit_test
    unit_test_handler
)

if(${CMAKE_SYSTEM_NAME} MATCHES Windows)
    target_compile_definitions(unit_test_handler PRIVATE WINDOWS=1 LINUX_OSX=0)
    add_definitions( -DFILESYSTEM_=std::filesystem )
    add_definitions( -DFILESYSTEM_INCLUDE_FILE_=<filesystem> )
elseif(${CMAKE_SYSTEM_NAME} MATCHES Linux)
    target_compile_definitions(unit_test_handler PRIVATE WINDOWS=0 LINUX_OSX=1)
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
        add_definitions( -DFILESYSTEM_=std::experimental::filesystem )
        add_definitions( -DFILESYSTEM_INCLUDE_FILE_=<experimental/filesystem> )
        target_link_libraries(unit_test_handler stdc++fs)
    else()
        add_definitions( -DFILESYSTEM_=std::filesystem )
        add_definitions( -DFILESYSTEM_INCLUDE_FILE_=<filesystem> )
    endif()
else()
    target_compile_definitions(unit_test_handler PRIVATE WINDOWS=0 LINUX_OSX=0)
endif()