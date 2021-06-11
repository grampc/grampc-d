
add_library(grampc 
  grampc/src/euler1.c
  grampc/src/eulermod2.c
  grampc/src/grampc_alloc.c
  grampc/src/grampc_fixedsize.c
  grampc/src/grampc_init.c
  grampc/src/grampc_mess.c
  grampc/src/grampc_run.c
  grampc/src/grampc_setopt.c
  grampc/src/grampc_setparam.c
  grampc/src/grampc_util.c
  grampc/src/heun2.c
  grampc/src/rodas.c
  grampc/src/ruku45.c
  grampc/src/simpson.c
  grampc/src/trapezodial.c
  grampc/cpp/src/grampc.cpp
  grampc/cpp/src/problem_description.cpp
)

target_include_directories(grampc PUBLIC
    ${CMAKE_SOURCE_DIR}/libs/grampc/include
    ${CMAKE_SOURCE_DIR}/libs/grampc/cpp/include)

# include asio
include(ExternalProject)
find_package(Git REQUIRED)

ExternalProject_Add(asio
	GIT_REPOSITORY https://github.com/chriskohlhoff/asio.git
	GIT_TAG 53d83a3321e62e54d352425479e59a91edd1a34b
	UPDATE_COMMAND ""
	CONFIGURE_COMMAND ""
    	BUILD_COMMAND ""
    	INSTALL_COMMAND ""
    	LOG_DOWNLOAD ON
        PREFIX  ${CMAKE_SOURCE_DIR}/libs/asio
)
ExternalProject_Get_Property(asio source_dir)
set(ASIO_INCLUDE_DIR ${source_dir}/asio/include)

# define flag to use asio as standalone library
add_definitions(-DASIO_STANDALONE )