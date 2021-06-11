
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

target_include_directories(grampcd PUBLIC 
    ${CMAKE_SOURCE_DIR}/libs/asio/src/asio/asio/include
    include)

add_dependencies(grampcd asio)

target_link_libraries(grampcd grampc factory)

if(MSVC)
else()
    target_link_libraries(grampcd pthread)
endif()

if(PYTHON_AVAILABLE)
	add_subdirectory(../libs/pybind11 ${CMAKE_SOURCE_DIR}/bin)

    pybind11_add_module(grampcd_interface SHARED 
        src/interface/python_interface.cpp 
    ) 

    target_link_libraries(grampcd_interface PRIVATE grampcd )
    target_link_libraries(grampcd_interface PRIVATE factory )
endif()