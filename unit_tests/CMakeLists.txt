
include_directories(include)

add_library(unit_test_handler
    src/checksum_handler.cpp
    src/file_handler.cpp
)

target_link_libraries(unit_test_handler
    dmpc
)

add_dependencies(unit_test_handler sha)

add_executable(run_unit_test
    src/run_unit_test.cpp
)

target_link_libraries(run_unit_test
    dmpc
    unit_test_handler
)

add_executable(setup_unit_test
    src/setup_unit_test.cpp
)

target_link_libraries( setup_unit_test
    unit_test_handler
)
