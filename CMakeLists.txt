#######################################################################
#                                                                     #
# CMake script for configuration of GalerkinTools library             #
#                                                                     #
# Usage: cmake -DCMAKE_INSTALL_PREFIX=/installation/path /source/path #
#                                                                     #
#######################################################################

# This may run with older versions, it is not tested however
cmake_minimum_required(VERSION 3.16)

project(galerkin_tools VERSION 0.1.0
                      DESCRIPTION "A set of tools based on deal.II for automatic assembly of finite element systems using the Galerkin method"
                      LANGUAGES CXX)

#locate deal.II
FIND_PACKAGE(deal.II 9.4.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "The version of deal.II required is 9.4.0. Higher versions may work. However, they are not tested yet.\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

#define paths for installation of libraries and headers of library galerkin_tools
SET(GALERKIN_TOOLS_INSTALL_PATH_LIBRARIES ${CMAKE_INSTALL_PREFIX}/lib)
SET(GALERKIN_TOOLS_INSTALL_PATH_HEADERS ${CMAKE_INSTALL_PREFIX}/include)

#set library names for debug and release versions of galerkin_tools library, respectively
SET(GALERKIN_TOOLS_TARGET_DEBUG ${PROJECT_NAME}.debug)
SET(GALERKIN_TOOLS_TARGET_RELEASE ${PROJECT_NAME}.release)

#this does the actual definition of the target, link libraries, etc.
ADD_SUBDIRECTORY(source)

#get the properties for compiling, linking and headers (these are written into the config file)
GET_PROPERTY(GALERKIN_TOOLS_LIBRARIES_DEBUG GLOBAL PROPERTY GALERKIN_TOOLS_LIBRARIES_DEBUG)
GET_PROPERTY(GALERKIN_TOOLS_LIBRARIES_RELEASE GLOBAL PROPERTY GALERKIN_TOOLS_LIBRARIES_RELEASE)
GET_PROPERTY(GALERKIN_TOOLS_INCLUDE_DIRECTORIES GLOBAL PROPERTY GALERKIN_TOOLS_INCLUDE_DIRECTORIES)
GET_PROPERTY(GALERKIN_TOOLS_COMPILE_FLAGS_DEBUG GLOBAL PROPERTY GALERKIN_TOOLS_COMPILE_FLAGS_DEBUG)
GET_PROPERTY(GALERKIN_TOOLS_COMPILE_FLAGS_RELEASE GLOBAL PROPERTY GALERKIN_TOOLS_COMPILE_FLAGS_RELEASE)
GET_PROPERTY(GALERKIN_TOOLS_COMPILE_DEFINITIONS_DEBUG GLOBAL PROPERTY GALERKIN_TOOLS_COMPILE_DEFINITIONS_DEBUG)
GET_PROPERTY(GALERKIN_TOOLS_COMPILE_DEFINITIONS_RELEASE GLOBAL PROPERTY GALERKIN_TOOLS_COMPILE_DEFINITIONS_RELEASE)
GET_PROPERTY(GALERKIN_TOOLS_LINK_FLAGS_DEBUG GLOBAL PROPERTY GALERKIN_TOOLS_LINK_FLAGS_DEBUG)
GET_PROPERTY(GALERKIN_TOOLS_LINK_FLAGS_RELEASE GLOBAL PROPERTY GALERKIN_TOOLS_LINK_FLAGS_RELEASE)

#configure the tests
#make sure that the tests know of the headers of the galerkin_tools library and are linked against the galerkin_tools library
INCLUDE_DIRECTORIES(${GALERKIN_TOOLS_INCLUDE_DIRECTORIES})
SET(TEST_LIBRARIES_DEBUG ${GALERKIN_TOOLS_LIBRARIES_DEBUG})
SET(TEST_LIBRARIES_RELEASE ${GALERKIN_TOOLS_LIBRARIES_RELEASE})
#set compile and link flags (use those of deal.II)
SET(DEAL_II_CXX_FLAGS "")
SET(DEAL_II_USER_DEFINITIONS "")
SET(DEAL_II_LINKER_FLAGS "")
FOREACH(_build "debug" "release")
  STRING(TOUPPER ${_build} _build_toupper)
  SET(DEAL_II_CXX_FLAGS_${_build_toupper} ${GALERKIN_TOOLS_COMPILE_FLAGS_${_build_toupper}})
  SET(DEAL_II_USER_DEFINITIONS_${_build_toupper} ${GALERKIN_TOOLS_COMPILE_DEFINITIONS_${_build_toupper}})
  SET(DEAL_II_LINKER_FLAGS_${_build_toupper} ${GALERKIN_TOOLS_LINK_FLAGS_${_build_toupper}})
ENDFOREACH(_build)
#build the tests for both, debug and release mode
SET(CMAKE_BUILD_TYPE "DebugRelease")
#enable testing
ENABLE_TESTING()
#this does the actual work
ADD_SUBDIRECTORY(tests)

#set up and install the configuration information of the library galerkin_tools
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/cmake/config/Config.cmake.in "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/galerkin_toolsConfig.cmake" @ONLY)
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/cmake/config/ConfigVersion.cmake.in "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/galerkin_toolsConfigVersion.cmake" @ONLY)
INSTALL(FILES ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/galerkin_toolsConfig.cmake DESTINATION ${GALERKIN_TOOLS_INSTALL_PATH_LIBRARIES}/cmake/galerkin_tools/)
INSTALL(FILES ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/galerkin_toolsConfigVersion.cmake DESTINATION ${GALERKIN_TOOLS_INSTALL_PATH_LIBRARIES}/cmake/galerkin_tools/)

