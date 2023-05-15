#source files from which library galerkin_tools is compiled
SET(_target_src
    assembly_helper.cc
    dependent_field.cc
    dirichlet_constraint.cc
    dof_handler_system.cc
    fe_values_interface.cc
    scalar_functional.cc
    independent_field.cc
    ldr.cc
    linear_material.cc
    tools.cc
    total_potential_contribution.cc
    total_potential.cc
    triangulation_system.cc
    solver_wrapper.cc
    dof_renumbering.cc
    two_block_sparsity_pattern.cc
    two_block_matrix.cc
    fe_dgq_anisotropic.cc
)

#link against lapacke and include header lapacke.h
FIND_LIBRARY(LAPACKE_LIBRARY lapacke HINTS ${LAPACKE_DIR} $ENV{LAPACKE_DIR} PATH_SUFFIXES lib)
FIND_PATH(LAPACKE_INSTALL_PATH_HEADERS NAMES lapacke.h HINTS ${LAPACKE_DIR} $ENV{LAPACKE_DIR} PATH_SUFFIXES include include/lapacke)
IF( (${LAPACKE_LIBRARY} MATCHES "LAPACKE_LIBRARY-NOTFOUND") OR (${LAPACKE_INSTALL_PATH_HEADERS} MATCHES "LAPACKE_INSTALL_PATH_HEADERS-NOTFOUND") )
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate lapacke ***\n\n"
    "You may want to either pass a flag -DLAPACKE_DIR=/path/to/lapacke to cmake\n"
    "or set an environment variable \"LAPACKE_DIR\" that contains this path.")
ELSE()
  LINK_LIBRARIES(${LAPACKE_LIBRARY})
  INCLUDE_DIRECTORIES(${LAPACKE_INSTALL_PATH_HEADERS})
ENDIF()

#find PARDISO
IF( DEFINED ENV{PARDISO_LIBRARY} )
  FIND_PACKAGE(OpenMP)
  IF(NOT ${Open_MP_FOUND})
    MESSAGE(FATAL_ERROR "\n"
      "*** You want to use PARDISO, which requires OpenMP. But did not find OpenMP.")
  ENDIF()
  #TODO: deal.II can in principle be configured without LAPACK. In this case, LAPACK would have to be pulled in here ...
  IF(NOT ${DEAL_II_WITH_LAPACK})
    MESSAGE(FATAL_ERROR "\n"
      "*** You want to use PARDISO, which requires LAPACK. But deal.II was not configured with LAPACK.")
  ENDIF()
  LINK_LIBRARIES($ENV{PARDISO_LIBRARY})
  OPTION(GALERKIN_TOOLS_WITH_PARDISO "Enable PARDISO" ON)
ELSE()
  MESSAGE("*** Could not locate PARDISO ***\n\n"
    "If you want to use PARDISO, you must set an environment variable \"PARDISO_LIBRARY\" that contains the path to the library (e.g. /path/to/pardiso.so).")
ENDIF()

#find MA57
IF( DEFINED ENV{MA57_LIBRARIES} )
  #TODO: deal.II can in principle be configured without LAPACK. In this case, LAPACK would have to be pulled in here ...
  IF(NOT ${DEAL_II_WITH_LAPACK})
    MESSAGE(FATAL_ERROR "\n"
      "*** You want to use MA57, which requires LAPACK. But deal.II was not configured with LAPACK.")
  ENDIF()
  LINK_LIBRARIES($ENV{MA57_LIBRARIES})
  OPTION(GALERKIN_TOOLS_WITH_MA57 "Enable MA57" ON)
ELSE()
  MESSAGE("*** Could not locate MA57 ***\n\n"
    "If you want to use MA57, you must set an environment variable \"MA57_LIBRARIES\" that contains the path to the MA57 library and Metis4 (e.g. MA57_LIBRARIES = /path/to/libma57.so;/path/to/libmetis.a).")
ENDIF()

#define the debug and release targets of the library galerkin_tools
FOREACH(_build "debug" "release")

    #set up the target
    STRING(TOUPPER ${_build} _build_toupper)
    SET(_target ${GALERKIN_TOOLS_TARGET_${_build_toupper}})
    ADD_LIBRARY(${_target} SHARED ${_target_src})

    #define compile flags and definitions and link flags (according to what has been used to compile deal.II)
    SET_PROPERTY(TARGET ${_target} APPEND_STRING PROPERTY COMPILE_FLAGS " ${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${_build_toupper}} -Wno-deprecated-declarations -Wno-cast-function-type")
    SET_PROPERTY(TARGET ${_target} APPEND PROPERTY COMPILE_DEFINITIONS "${DEAL_II_USER_DEFINITIONS};${DEAL_II_USER_DEFINITIONS_${_build_toupper}}")
    SET_PROPERTY(TARGET ${_target} APPEND PROPERTY LINK_FLAGS "${DEAL_II_LINKER_FLAGS} ${DEAL_II_LINKER_FLAGS_${_build_toupper}}")
    IF(GALERKIN_TOOLS_WITH_PARDISO)
      SET_PROPERTY(TARGET ${_target} APPEND_STRING PROPERTY COMPILE_FLAGS " -fopenmp")
      SET_PROPERTY(TARGET ${_target} APPEND_STRING PROPERTY LINK_FLAGS " -fopenmp")
      #the pardiso library is not linked against lapack, blas and openmp - so make sure that these libraries are nevertheless listed in the shared library (otherwise they would have
      #to be included into LD_PRELOAD
      SET_PROPERTY(TARGET ${_target} APPEND_STRING PROPERTY LINK_FLAGS " -Wl,--no-as-needed")
    ENDIF()

    #link against deal.II
    TARGET_LINK_LIBRARIES(${_target} ${DEAL_II_LIBRARIES_${_build_toupper}})

    #install path for the library GalerkinTools (debug and release)
    INSTALL(TARGETS ${_target} LIBRARY DESTINATION ${GALERKIN_TOOLS_INSTALL_PATH_LIBRARIES})

    #the libraries to link to when using GalerkinTools (this will end up in the config file)
    GET_TARGET_PROPERTY(LINK_LIBRARIES_ ${_target} LINK_LIBRARIES)
    SET(GALERKIN_TOOLS_LIBRARIES "${GALERKIN_TOOLS_INSTALL_PATH_LIBRARIES}/libgalerkin_tools.${_build}.so")
    FOREACH(LINK_LIBRARY ${LINK_LIBRARIES_})
      STRING(APPEND GALERKIN_TOOLS_LIBRARIES ";${LINK_LIBRARY}")
    ENDFOREACH(LINK_LIBRARY)
    SET_PROPERTY(GLOBAL PROPERTY GALERKIN_TOOLS_LIBRARIES_${_build_toupper} ${GALERKIN_TOOLS_LIBRARIES})

    #flags (this will end up in the config file)
    GET_TARGET_PROPERTY(GALERKIN_TOOLS_COMPILE_FLAGS ${_target} COMPILE_FLAGS)
    SET_PROPERTY(GLOBAL PROPERTY GALERKIN_TOOLS_COMPILE_FLAGS_${_build_toupper} ${GALERKIN_TOOLS_COMPILE_FLAGS})
    GET_TARGET_PROPERTY(GALERKIN_TOOLS_COMPILE_DEFINITIONS ${_target} COMPILE_DEFINITIONS)
    SET_PROPERTY(GLOBAL PROPERTY GALERKIN_TOOLS_COMPILE_DEFINITIONS_${_build_toupper} ${GALERKIN_TOOLS_COMPILE_DEFINITIONS})
    GET_TARGET_PROPERTY(GALERKIN_TOOLS_LINK_FLAGS ${_target} LINK_FLAGS)
    SET_PROPERTY(GLOBAL PROPERTY GALERKIN_TOOLS_LINK_FLAGS_${_build_toupper} ${GALERKIN_TOOLS_LINK_FLAGS})

ENDFOREACH(_build)
INCLUDE_DIRECTORIES(${DEAL_II_INCLUDE_DIRS})

#the include directories needed when using GalerkinTools (this will end up in the config file)
GET_PROPERTY(INCLUDE_DIRECTORIES_ DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
SET(GALERKIN_TOOLS_INCLUDE_DIRECTORIES "${GALERKIN_TOOLS_INSTALL_PATH_HEADERS}")
FOREACH(INCLUDE_DIRECTORY ${INCLUDE_DIRECTORIES_})
  STRING(APPEND GALERKIN_TOOLS_INCLUDE_DIRECTORIES ";${INCLUDE_DIRECTORY}")
ENDFOREACH(INCLUDE_DIRECTORY)
SET_PROPERTY(GLOBAL PROPERTY GALERKIN_TOOLS_INCLUDE_DIRECTORIES ${GALERKIN_TOOLS_INCLUDE_DIRECTORIES})

#headers of library galerkin_tools itself
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../include)

#make config.h
CONFIGURE_FILE(../include/galerkin_tools/config.h.in ../include/galerkin_tools/config.h)

#install the headers of the library galerkin_tools
INSTALL(DIRECTORY ${CMAKE_SOURCE_DIR}/include/galerkin_tools DESTINATION ${GALERKIN_TOOLS_INSTALL_PATH_HEADERS} FILES_MATCHING PATTERN "*.h")
INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../include/galerkin_tools DESTINATION ${GALERKIN_TOOLS_INSTALL_PATH_HEADERS} FILES_MATCHING PATTERN "*.h")

