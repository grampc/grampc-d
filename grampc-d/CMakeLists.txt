# This file is part of GRAMPC-D - (https://github.com/grampc-d/grampc-d.git)
#
# GRAMPC-D -- A software framework for distributed model predictive control (DMPC)
# 
#
# Copyright 2023 by Daniel Burk, Maximilian Pierer von Esch, Andreas Voelz, Knut Graichen
# All rights reserved.
#
# GRAMPC-D is distributed under the BSD-3-Clause license, see LICENSE.txt

# compile all source files into library
add_library(grampcd
    src/agent/agent.cpp
    src/agent/neighbor.cpp
    src/coord/coordinator.cpp
    src/simulator/simulator.cpp
      
    src/interface/interface.cpp
    src/interface/dmpc_interface.cpp

    src/comm/communication_interface.cpp
    src/comm/communication_interface_central.cpp  
    src/comm/communication_interface_local.cpp    
    src/comm/message_handler.cpp    
    src/comm/message.cpp    
    src/comm/message_definition.cpp

    src/model/agent_model.cpp
    src/model/coupling_model.cpp
    src/model/model_factory.cpp

    src/optim/optim_util.cpp
    src/util/data_conversion.cpp
    
    src/util/logging.cpp
    src/optim/solution.cpp

    src/optim/solver_central.cpp
    src/optim/solver_local.cpp
    src/optim/solver_local_admm.cpp
    src/optim/solver_local_sensi.cpp

    src/optim/problem_description_central.cpp
    src/optim/problem_description_local_default.cpp
    src/optim/problem_description_local_neighbor_approximation.cpp
    src/optim/problem_description_local_sensi.cpp
    src/optim/approximate_neighbor.cpp
    src/util/auto_tune.cpp
    src/agent/step_selector.cpp
    src/agent/sync_step_selector.cpp
    src/agent/async_step_selector.cpp
 )

target_include_directories(grampcd PUBLIC 
    include
    ${ASIO_INCLUDE_DIR}
    ${CEREAL_INCLUDE_DIR}
    ${GALGO_INCLUDE_DIR}
)

add_dependencies(grampcd 
    asio
    cereal
)

# define flag to use asio as standalone library
target_compile_definitions(grampcd PRIVATE -DASIO_STANDALONE )
if(AUTO_TUNE)
    target_compile_definitions(grampcd PRIVATE -DAUTO_TUNE)
endif()

target_link_libraries(grampcd 
    grampc 
    factory
)

if(MSVC)
else()
    target_link_libraries(grampcd pthread)
endif()

if(PYTHON_AVAILABLE)
	add_subdirectory(../libs/pybind11 ${CMAKE_SOURCE_DIR}/bin)

    pybind11_add_module(grampcd_interface SHARED 
        src/interface/python_interface.cpp 
    ) 

    target_link_libraries(grampcd_interface PRIVATE
        grampcd
        factory
    )
endif()