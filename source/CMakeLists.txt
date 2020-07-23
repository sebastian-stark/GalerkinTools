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
)

#define the debug and release targets of the library galerkin_tools
FOREACH(_build "debug" "release")

    #set up the target
    STRING(TOUPPER ${_build} _build_toupper)
    SET(_target ${GALERKIN_TOOLS_TARGET_${_build_toupper}})
    ADD_LIBRARY(${_target} SHARED ${_target_src})

    #define compile flags and definitions (according to what has been used to compile deal.II)
    SET_PROPERTY(TARGET ${_target} APPEND_STRING PROPERTY COMPILE_FLAGS " ${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${_build_toupper}}")
    SET_PROPERTY(TARGET ${_target} APPEND PROPERTY COMPILE_DEFINITIONS "${DEAL_II_USER_DEFINITIONS};${DEAL_II_USER_DEFINITIONS_${_build_toupper}}")
    SET_PROPERTY(TARGET ${_target} APPEND PROPERTY LINK_FLAGS "${DEAL_II_LINKER_FLAGS} ${DEAL_II_LINKER_FLAGS_${_build_toupper}}")

    #link against lapacke and include header lapacke.h
    FIND_LIBRARY(LAPACKE_LIBRARY lapacke HINTS ${LAPACKE_DIR} $ENV{LAPACKE_DIR} PATH_SUFFIXES lib)
    FIND_PATH(LAPACKE_INSTALL_PATH_HEADERS NAMES lapacke.h HINTS ${LAPACKE_DIR} $ENV{LAPACKE_DIR} PATH_SUFFIXES include include/lapacke)
    IF( (${LAPACKE_LIBRARY} MATCHES "LAPACKE_LIBRARY-NOTFOUND") OR (${LAPACKE_INSTALL_PATH_HEADERS} MATCHES "LAPACKE_INSTALL_PATH_HEADERS-NOTFOUND") )
      MESSAGE(FATAL_ERROR "\n"
        "*** Could not locate lapacke ***\n\n"
        "You may want to either pass a flag -DLAPACKE_DIR=/path/to/lapacke to cmake\n"
        "or set an environment variable \"LAPACKE_DIR\" that contains this path.")
    ELSE()
      TARGET_LINK_LIBRARIES(${_target} ${LAPACKE_LIBRARY})
      TARGET_INCLUDE_DIRECTORIES(${_target} PRIVATE ${LAPACKE_INSTALL_PATH_HEADERS})
      OPTION(GALERKIN_TOOLS_WITH_UMFPACK "Enable UMFPACK" ON)
    ENDIF()

    #link against UMFPACK and include corresponding header
    FIND_LIBRARY(UMFPACK_LIBRARY umfpack HINTS ${SUITE_SPARSE_DIR} $ENV{SUITE_SPARSE_DIR} ${UMFPACK_DIR} $ENV{UMFPACK_DIR} PATH_SUFFIXES lib)
    FIND_PATH(UMFPACK_INSTALL_PATH_HEADERS NAMES umfpack.h HINTS ${SUITE_SPARSE_DIR} $ENV{SUITE_SPARSE_DIR} ${UMFPACK_DIR} $ENV{UMFPACK_DIR} PATH_SUFFIXES include)
    IF( (${UMFPACK_LIBRARY} MATCHES "UMFPACK_LIBRARY-NOTFOUND") OR (${UMFPACK_INSTALL_PATH_HEADERS} MATCHES "UMFPACK_INSTALL_PATH_HEADERS-NOTFOUND") )
      MESSAGE("*** Could not locate UMFPACK ***\n\n"
        "You may want to either pass a flag -DSUITE_SPARSE_DIR=/path/to/suite_sparse to cmake\n"
        "or set an environment variable \"SUITE_SPARSE_DIR\" that contains this path.")
    ELSE()
      TARGET_LINK_LIBRARIES(${_target} ${UMFPACK_LIBRARY})
      TARGET_INCLUDE_DIRECTORIES(${_target} PRIVATE ${UMFPACK_INSTALL_PATH_HEADERS})
    ENDIF()

    #link against deal.II and include headers of deal.II
    TARGET_LINK_LIBRARIES(${_target} ${DEAL_II_LIBRARIES_${_build_toupper}})
    TARGET_INCLUDE_DIRECTORIES(${_target} PRIVATE ${DEAL_II_INCLUDE_DIRS})

    #headers of library galerkin_tools
    TARGET_INCLUDE_DIRECTORIES(${_target} PRIVATE ${CMAKE_SOURCE_DIR}/include)

    #install the library galerkin_tools (debug and release)
    INSTALL(TARGETS ${_target} LIBRARY DESTINATION ${GALERKIN_TOOLS_INSTALL_PATH_LIBRARIES})

ENDFOREACH(_build)

#make config.h
CONFIGURE_FILE(../include/galerkin_tools/config.h.in ../include/galerkin_tools/config.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/../include)

#install the headers of the library galerkin_tools
INSTALL(DIRECTORY ${CMAKE_SOURCE_DIR}/include/galerkin_tools DESTINATION ${GALERKIN_TOOLS_INSTALL_PATH_HEADERS} FILES_MATCHING PATTERN "*.h")
INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../include/galerkin_tools DESTINATION ${GALERKIN_TOOLS_INSTALL_PATH_HEADERS} FILES_MATCHING PATTERN "*.h")

