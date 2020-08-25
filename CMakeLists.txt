# This file is part of GRAMPC-D - (https://github.com/DanielBurk/GRAMPC-D.git)
#
# GRAMPC-D -- A software framework for distributed model predictive control (DMPC)
# based on the alternating direction method of multipliers (ADMM).
#
# Copyright 2020 by Daniel Burk, Andreas Voelz, Knut Graichen
# All rights reserved.
#
# GRAMPC-D is distributed under the BSD-3-Clause license, see LICENSE.txt

cmake_minimum_required(VERSION 2.6)
project(grampcd)

set(PYTHON_AVAILABLE FALSE)

# add support for C++17 standard
if(MSVC)
    SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} /std:c++17 /DASIO_STANDALONE -D_WIN32_WINNT=0x0A00")
    SET(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ")
    set(CMAKE_CXX_STANDARD 17)
else()
    SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++17 -DASIO_STANDALONE -Wno-psabi")
    SET(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ")
    set(CMAKE_CXX_STANDARD 17)
endif()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# link_directories should treat paths as relative to the source directory
cmake_policy(SET CMP0015 NEW)


# where to search for header files
include_directories(
    include
    model_description
    libs/grampc/include
    libs/grampc/cpp/include
    libs/asio/asio/include
)

# add simulation examples
add_subdirectory(main)

# add model descriptions
add_subdirectory(model_description)

# compile all source files into library
add_library(dmpc
    src/agent/agent.cpp
    src/agent/neighbor.cpp
    src/coord/coordinator.cpp
    src/simulator/simulator.cpp
      
    src/interface/interface.cpp
    src/interface/dmpc_interface.cpp

    src/comm/communication_interface.cpp
    src/comm/communication_interface_central.cpp  
    src/comm/communication_interface_local.cpp    

    src/model/agent_model.cpp
    src/model/coupling_model.cpp
    src/model/model_factory.cpp

    src/optim/optim_util.cpp
    src/util/data_conversion.cpp
    src/util/protocol_communication.cpp
    src/util/logging.cpp
    src/optim/solution.cpp

    src/optim/solver_central.cpp
    src/optim/solver_local.cpp

    src/optim/problem_description_central.cpp
    src/optim/problem_description_local_default.cpp
    src/optim/problem_description_local_neighbor_approximation.cpp
    src/optim/approximate_neighbor.cpp
)

add_library(grampc
  libs/grampc/src/euler1.c
  libs/grampc/src/eulermod2.c
  libs/grampc/src/grampc_alloc.c
  libs/grampc/src/grampc_fixedsize.c
  libs/grampc/src/grampc_init.c
  libs/grampc/src/grampc_mess.c
  libs/grampc/src/grampc_run.c
  libs/grampc/src/grampc_setopt.c
  libs/grampc/src/grampc_setparam.c
  libs/grampc/src/grampc_util.c
  libs/grampc/src/heun2.c
  libs/grampc/src/rodas.c
  libs/grampc/src/ruku45.c
  libs/grampc/src/simpson.c
  libs/grampc/src/trapezodial.c
  libs/grampc/cpp/src/grampc.cpp
  libs/grampc/cpp/src/problem_description.cpp
)

target_link_libraries(dmpc grampc factory)

if(MSVC)
else()
    target_link_libraries(dmpc pthread)
endif()

if(PYTHON_AVAILABLE)
	add_subdirectory(libs/pybind11)

    pybind11_add_module(grampcd_interface SHARED 
        src/interface/python_interface.cpp 
    ) 

    target_link_libraries(grampcd_interface PRIVATE dmpc )
    target_link_libraries(grampcd_interface PRIVATE factory )
endif()

