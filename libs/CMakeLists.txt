# This file is part of GRAMPC-D - (https://github.com/grampc-d/grampc-d.git)
#
# GRAMPC-D -- A software framework for distributed model predictive control (DMPC)
# based on the alternating direction method of multipliers (ADMM).
#
# Copyright 2020 by Daniel Burk, Andreas Voelz, Knut Graichen
# All rights reserved.
#
# GRAMPC-D is distributed under the BSD-3-Clause license, see LICENSE.txt

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
    grampc/include
    grampc/cpp/include)
    
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
set(ASIO_INCLUDE_DIR ${source_dir}/asio/include CACHE INTERNAL "ASIO_INCLUDE_DIR")

# include SHA1 as external project
ExternalProject_Add(cereal
    GIT_REPOSITORY https://github.com/USCiLab/cereal.git
    GIT_TAG 46a4a910077bf9e9f8327c8f6ea761c89b06da53
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
  	    BUILD_COMMAND ""
  	    INSTALL_COMMAND ""
	    LOG_DOWNLOAD ON
        PREFIX  ${CMAKE_SOURCE_DIR}/libs/cereal
)
ExternalProject_Get_Property(cereal SOURCE_DIR)
set(CEREAL_INCLUDE_DIR ${SOURCE_DIR}/include CACHE INTERNAL "CEREAL_INCLUDE_DIR")

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

ExternalProject_Get_Property(sha SOURCE_DIR)
set(SHA_INCLUDE_DIR ${SOURCE_DIR} CACHE INTERNAL "SHA_INCLUDE_DIR")

# include galgo as external project
if(AUTO_TUNE)
    ExternalProject_Add(galgo
        GIT_REPOSITORY https://github.com/olmallet81/GALGO-2.0.git
        GIT_TAG 6ee1cd951a2c5abbbcd7cdadf6598eac42e3530b
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND ""
  	        BUILD_COMMAND ""
  	        INSTALL_COMMAND ""
	        LOG_DOWNLOAD ON
            PREFIX  ${CMAKE_SOURCE_DIR}/libs/galgo
    )

    ExternalProject_Get_Property(galgo SOURCE_DIR)
    set(GALGO_INCLUDE_DIR ${SOURCE_DIR}/src CACHE INTERNAL "GALGO_INCLUDE_DIR")
endif()