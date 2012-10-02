#####
##
## Copyright (C) 2012 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/template-arguments.in
  ${CMAKE_CURRENT_BINARY_DIR}/template-arguments
  )

#
# Setup and install the convenience macros:
#

LIST(APPEND DEAL_II_MACROS
  ${CMAKE_INSTALL_PREFIX}/${DEAL_II_CMAKE_MACROS_RELDIR}/macro_deal_ii_initialize_cached_variables.cmake
  ${CMAKE_INSTALL_PREFIX}/${DEAL_II_CMAKE_MACROS_RELDIR}/macro_deal_ii_setup_target.cmake
  )

INSTALL(FILES
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_initialize_cached_variables.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_setup_target.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_setup_compiler.cmake
  DESTINATION ${DEAL_II_CMAKE_MACROS_RELDIR}
  COMPONENT library
  )

#
# Configure and install the cmake project configuration:
#

IF(NOT "${CONFIG_EXTERNAL_LIBRARIES}" STREQUAL "")
  SET(CONFIG_EXTERNAL_LIBRARIES
    general
    ${DEAL_II_EXTERNAL_LIBRARIES}
    )
ENDIF()

IF(CMAKE_BUILD_TYPE MATCHES "Debug")
  SET(CONFIG_EXTERNAL_LIBRARIES_DEBUG
    ${DEAL_II_EXTERNAL_LIBRARIES}
    ${DEAL_II_EXTERNAL_LIBRARIES_DEBUG}
    )
  LIST(APPEND CONFIG_LIBRARIES
    debug
    ${DEAL_II_EXTERNAL_LIBRARIES_DEBUG}
    )
ENDIF()

IF(CMAKE_BUILD_TYPE MATCHES "Release")
  SET(CONFIG_EXTERNAL_LIBRARIES_RELEASE
    ${DEAL_II_EXTERNAL_LIBRARIES}
    ${DEAL_II_EXTERNAL_LIBRARIES_RELEASE}
    )
  LIST(APPEND CONFIG_LIBRARIES
    optimized
    ${DEAL_II_EXTERNAL_LIBRARIES_RELEASE}
    )
ENDIF()

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  )

GET_CMAKE_PROPERTY(res VARIABLES)
FOREACH(var ${res})
  IF(var MATCHES "DEAL_II_WITH")
    FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
      "SET(${var} ${${var}})\n"
      )
  ENDIF()
ENDFOREACH()

#
# Include the convenience macros directly at the end of our project
# config:
#
FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  "\n#\n# Include some convenience macros directly in the project config:\n#\n"
  )
FOREACH(var ${DEAL_II_MACROS})
  FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
    "INCLUDE(${var})\n"
    )
ENDFOREACH()

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/ConfigVersion.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}ConfigVersion.cmake
  @ONLY
  )

INSTALL(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}ConfigVersion.cmake
  DESTINATION ${DEAL_II_PROJECT_CONFIG_RELDIR}
  COMPONENT library
  )


IF(DEAL_II_COMPONENT_COMPAT_FILES)
  #
  # Configure and install the old Makefile configuration:
  #

  #
  # Therefore, transform some cmake lists into a string that the old
  # Makefile mechanism actually understands:
  #

  TO_STRING_AND_ADD_PREFIX(MAKEFILE_INCLUDE_DIRS "${CMAKE_INCLUDE_FLAG_C}"
    ${DEAL_II_USER_INCLUDE_DIRS}
    ${DEAL_II_INCLUDE_DIRS}
    )

  FOREACH(build ${DEAL_II_BUILD_TYPES})
    TO_STRING_AND_ADD_PREFIX(MAKEFILE_DEFINITIONS_${build}
      "-D"
      ${DEAL_II_USER_DEFINITIONS}
      ${DEAL_II_USER_DEFINITIONS_${build}}
      )

    #
    # Add an rpath directive in front of each library, so that libraries
    # outside of the default search directories will be found by the runtime
    # dynamic linker:
    #
    FOREACH(lib
        ${DEAL_II_EXTERNAL_LIBRARIES}
        ${DEAL_II_EXTERNAL_LIBRARIES_${build}}
        ${DEAL_II_LIBRARIES_${build}})
      GET_FILENAME_COMPONENT(path ${lib} PATH)
      LIST(APPEND MAKEFILE_LIBS_${build} "-Wl,-rpath -Wl,${path}")
      LIST(APPEND MAKEFILE_LIBS_${build} ${lib})
    ENDFOREACH()
    TO_STRING(MAKEFILE_LIBS_${build} ${MAKEFILE_LIBS_${build}})
    #
    # Put our linker flags directly in front of this string:
    #
    SET(MAKEFILE_LIBS_${build}
      "${CMAKE_SHARED_LINKER_FLAGS} ${DEAL_II_SHARED_LINKER_FLAGS_${build}} ${MAKEFILE_LIBS_${build}}"
      )
  ENDFOREACH()

  #
  # Boilerplate: The Make.global_options expects variables to be set to
  # yes, as is common for Makefiles.
  #
  COND_SET_TO_YES(DEAL_II_WITH_TBB MAKEFILE_enablethreads)
  COND_SET_TO_YES(DEAL_II_WITH_FUNCTIONPARSER MAKEFILE_enableparser)
  COND_SET_TO_YES(BUILD_SHARED_LIBS MAKEFILE_enableshared)

  COND_SET_TO_YES(DEAL_II_WITH_PETSC MAKEFILE_PETSC)
  COND_SET_TO_YES(DEAL_II_USE_PETSC_DEV MAKEFILE_PETSC_DEV)
  COND_SET_TO_YES(DEAL_II_WITH_TRILINOS MAKEFILE_TRILINOS)
  COND_SET_TO_YES(DEAL_II_WITH_BLAS MAKEFILE_BLAS)
  COND_SET_TO_YES(DEAL_II_WITH_LAPACK MAKEFILE_LAPACK)
  COND_SET_TO_YES(DEAL_II_WITH_ARPACK MAKEFILE_ARPACK)
  COND_SET_TO_YES(DEAL_II_WITH_METIS MAKEFILE_METIS)
  COND_SET_TO_YES(DEAL_II_WITH_UMFPACK MAKEFILE_UMFPACK)
  COND_SET_TO_YES(DEAL_II_WITH_P4EST MAKEFILE_P4EST)
  COND_SET_TO_YES(DEAL_II_WITH_MPI MAKEFILE_MPI)

  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/Make.global_options.in
    ${CMAKE_CURRENT_BINARY_DIR}/Make.global_options
    )

  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/Version.in
    ${CMAKE_CURRENT_BINARY_DIR}/Version
    )

  INSTALL(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/template-arguments
    ${CMAKE_CURRENT_BINARY_DIR}/Make.global_options
    DESTINATION common
    COMPONENT compat_files
    )

  INSTALL(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/Version
    DESTINATION ${CMAKE_INSTALL_PREFIX}
    COMPONENT compat_files
    )

ENDIF()

