#######################################################################
#                                                                     #
# CMake script for configuration of micromorphic_fibre                #
#                                                                     #
# Usage: cmake -DCMAKE_BUILD_TYPE=Debug                               #
#                    or                                               #
#        cmake -DCMAKE_BUILD_TYPE=Release                             #
#                                                                     #
#######################################################################

CMAKE_MINIMUM_REQUIRED(VERSION 3.5)

PROJECT("DEAL.II STEP 46")

#name of executable and source files from which executable is compiled
SET(_target "dealii_step46")
SET(_target_src ${_target}.cc)

#Determine whether to build in Release or Debug mode
IF(${CMAKE_BUILD_TYPE} MATCHES "Release")
    MESSAGE("Building target in Release mode")
    SET(_build "Release")
    ADD_DEFINITIONS(-DRELEASE)
ELSE()
    MESSAGE("Building target in Debug mode")
    SET(_build "Debug")
    ADD_DEFINITIONS(-DDEBUG)
ENDIF()
STRING(TOUPPER ${_build} _build_toupper)

#locate deal.II
FIND_PACKAGE(deal.II 9.1.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

#locate galerkin_tools
FIND_PACKAGE(galerkin_tools 1.0.0 QUIET HINTS $ENV{GALERKIN_TOOLS_DIR} ${GALERKIN_TOOLS_DIR})
IF(NOT ${galerkin_tools_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of galerkin_tools. ***\n\n"
    "You may want to either pass a flag -DGALERKIN_TOOLS_DIR=/path/to/galerkin_tools to cmake\n"
    "or set an environment variable \"GALERKIN_TOOLS_DIR\" that contains this path."
    )
ENDIF()

#locate incremental_fe
FIND_PACKAGE(incremental_fe 1.0.0 QUIET HINTS $ENV{INCREMENTAL_FE_DIR} ${INCREMENTAL_FE_DIR})
IF(NOT ${incremental_fe_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of incremental_fe. ***\n\n"
    "You may want to either pass a flag -DINCREMENTAL_FE_DIR=/path/to/incremental_fe to cmake\n"
    "or set an environment variable \"INCREMENTAL_FE_DIR\" that contains this path."
    )
ENDIF()

#set up the target
ADD_EXECUTABLE(${_target} ${_target_src})

#define compile flags and definitions (according to what has been used to compile incrementalFE)
SET_PROPERTY(TARGET ${_target} APPEND_STRING PROPERTY COMPILE_FLAGS " ${INCREMENTAL_FE_COMPILE_FLAGS_${_build_toupper}}")
SET_PROPERTY(TARGET ${_target} APPEND PROPERTY COMPILE_DEFINITIONS "${INCREMENTAL_FE_COMPILE_DEFINITIONS_${_build_toupper}}")

#link against deal.II, include headers of deal.II
TARGET_LINK_LIBRARIES(${_target} ${DEAL_II_LIBRARIES_${_build_toupper}})
TARGET_INCLUDE_DIRECTORIES(${_target} PRIVATE ${DEAL_II_INCLUDE_DIRS})

#link against incremental_fe, include headers of incremental_fe
#work around the cache of cmake
UNSET(INCREMENTAL_FE_LIBRARY CACHE)
FIND_LIBRARY(INCREMENTAL_FE_LIBRARY ${INCREMENTAL_FE_TARGET_${_build_toupper}} HINTS ${INCREMENTAL_FE_INSTALL_PATH_LIBRARIES})
IF(${INCREMENTAL_FE_LIBRARY} MATCHES "INCREMENTAL_FE_LIBRARY-NOTFOUND")
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate incremental_fe ***\n")
ELSE()
  TARGET_LINK_LIBRARIES(${_target} ${INCREMENTAL_FE_LIBRARY})
  TARGET_INCLUDE_DIRECTORIES(${_target} PRIVATE ${INCREMENTAL_FE_INSTALL_PATH_HEADERS})
ENDIF()

#link against galerkin_tools, include headers of galerkin_tools
#work around the cache of cmake
UNSET(GALERKIN_TOOLS_LIBRARY CACHE)
FIND_LIBRARY(GALERKIN_TOOLS_LIBRARY ${GALERKIN_TOOLS_TARGET_${_build_toupper}} HINTS ${GALERKIN_TOOLS_INSTALL_PATH_LIBRARIES})
IF(${GALERKIN_TOOLS_LIBRARY} MATCHES "GALERKIN_TOOLS_LIBRARY-NOTFOUND")
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate galerkin_tools ***\n")
ELSE()
  TARGET_LINK_LIBRARIES(${_target} ${GALERKIN_TOOLS_LIBRARY})
  TARGET_INCLUDE_DIRECTORIES(${_target} PRIVATE ${GALERKIN_TOOLS_INSTALL_PATH_HEADERS})
ENDIF()
