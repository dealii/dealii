## ---------------------------------------------------------------------
## $Id$
##
## Copyright (C) 2012 - 2013 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

#
# This file sets up the project configuration consisting of
#
#   Make.global_options (legacy)
#   deal.IIConfig.cmake
#   deal.IIVersionConfig.cmake
#
# and copies it (a) to the build directory and (b) prepares it for later
# installation.
#

#
# Read in auxiliary include directories for the build directory
# configuration:
#
FILE(STRINGS
  ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/deal_ii_source_includes
  deal_ii_source_includes
  )

#
# Configure the template-arguments file
#
CONFIGURE_FILE( # for binary dir:
  ${CMAKE_CURRENT_SOURCE_DIR}/template-arguments.in
  ${CMAKE_BINARY_DIR}/${DEAL_II_COMMON_RELDIR}/template-arguments
  )


########################################################################
#                                                                      #
#              Setup and install the convenience macros:               #
#                                                                      #
########################################################################

FILE(COPY # for binary dir:
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_initialize_cached_variables.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_invoke_autopilot.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_setup_target.cmake
  DESTINATION ${CMAKE_BINARY_DIR}/${DEAL_II_CMAKE_MACROS_RELDIR}
  )
INSTALL(FILES # for installation:
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_initialize_cached_variables.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_invoke_autopilot.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_setup_target.cmake
  DESTINATION ${DEAL_II_CMAKE_MACROS_RELDIR}
  COMPONENT library
  )


########################################################################
#                                                                      #
#        Configure and install the cmake project configuration:        #
#                                                                      #
########################################################################

#
# Strip -Wno-deprecated from CMAKE_CXX_FLAGS so that deprecation warnings
# are actually shown for user code:
#
STRIP_FLAG(CMAKE_CXX_FLAGS "-Wno-deprecated")
STRIP_FLAG(CMAKE_CXX_FLAGS "-Wno-deprecated-declarations")

#
# Populate a bunch of CONFIG_* variables with useful information:
#

FOREACH(_build ${DEAL_II_BUILD_TYPES})

  IF(NOT "${CMAKE_BUILD_TYPE}" STREQUAL "DebugRelease")
    SET(_keyword "general")
  ELSE()
    IF(_build MATCHES DEBUG)
      SET(_keyword "debug")
    ELSE()
      SET(_keyword "optimized")
    ENDIF()
  ENDIF()

  #
  # Get library name directly from the target:
  #
  GET_TARGET_PROPERTY(_lib ${DEAL_II_BASE_NAME}${DEAL_II_${_build}_SUFFIX} LOCATION)
  GET_FILENAME_COMPONENT(CONFIG_LIBRARY_${_build} "${_lib}" NAME)

  SET(CONFIG_LIBRARIES_${_build}
    \${DEAL_II_PATH}/\${DEAL_II_LIBRARY_RELDIR}/${CONFIG_LIBRARY_${_build}}
    ${DEAL_II_EXTERNAL_LIBRARIES_${_build}}
    ${DEAL_II_EXTERNAL_LIBRARIES}
    )
  LIST(APPEND CONFIG_LIBRARIES
    ${_keyword}
    \${DEAL_II_LIBRARIES_${_build}}
    )

  SET(CONFIG_TARGET_${_build}
    ${DEAL_II_BASE_NAME}${DEAL_II_${_build}_SUFFIX}
    )
  LIST(APPEND CONFIG_TARGET
    ${_keyword}
    \${DEAL_II_TARGET_${_build}}
    )
ENDFOREACH()

#
# For binary dir:
#
SET(CONFIG_INCLUDE_DIRS_BINARY
  \${DEAL_II_PATH}/include
  \${DEAL_II_PATH}/include/deal.II
  ${CMAKE_SOURCE_DIR}/include/
  ${CMAKE_SOURCE_DIR}/include/deal.II
  ${deal_ii_source_includes}
  ${DEAL_II_USER_INCLUDE_DIRS}
  )
#
# For installation:
#
SET(CONFIG_INCLUDE_DIRS_INSTALLATION
  \${DEAL_II_PATH}/\${DEAL_II_INCLUDE_RELDIR}
  \${DEAL_II_PATH}/\${DEAL_II_INCLUDE_RELDIR}/deal.II
  \${DEAL_II_PATH}/\${DEAL_II_INCLUDE_RELDIR}/deal.II/bundled
  ${DEAL_II_USER_INCLUDE_DIRS}
  )
CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  @ONLY
  )
#
# Append feature configuration:
#
SET(_file
  ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  )
FILE(APPEND ${_file}
  "\n\n#\n# Feature configuration:\n#\n\n"
  )
GET_CMAKE_PROPERTY(res VARIABLES)
FOREACH(var ${res})
  IF(var MATCHES "DEAL_II_WITH")
    FILE(APPEND ${_file}
      "SET(${var} ${${var}})\n"
      )
  ENDIF()
ENDFOREACH()

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/ConfigVersion.cmake.in
  ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}ConfigVersion.cmake
  @ONLY
  )

INSTALL(FILES
  ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}ConfigVersion.cmake
  DESTINATION ${DEAL_II_PROJECT_CONFIG_RELDIR}
  COMPONENT library
  )


########################################################################
#                                                                      #
#         Configure and install the old Makefile configuration:        #
#                                                                      #
########################################################################


IF(DEAL_II_COMPONENT_COMPAT_FILES)
  #
  # Transform some cmake lists into a string that the old Makefile
  # mechanism actually understands:
  #

  FOREACH(build ${DEAL_II_BUILD_TYPES})
    TO_STRING_AND_ADD_PREFIX(MAKEFILE_DEFINITIONS_${build}
      "-D"
      ${DEAL_II_USER_DEFINITIONS}
      ${DEAL_II_USER_DEFINITIONS_${build}}
      )

    #
    # Ideally, DEAL_II_EXTERNAL_LIBRARIES should just contain a list of
    # full library paths. Nevertheless, filter out all elements that are
    # not a full path starting with "/".
    #
    SET(_libs)
    FOREACH(_lib
        ${DEAL_II_EXTERNAL_LIBRARIES_${build}}
        ${DEAL_II_EXTERNAL_LIBRARIES}
        )
      IF(_lib MATCHES "^/")
        LIST(APPEND _libs "${_lib}")
      ENDIF()
    ENDFOREACH()

    TO_STRING(MAKEFILE_TARGETS_${build}
        $(D)/${DEAL_II_LIBRARY_RELDIR}/${CONFIG_LIBRARY_${build}}
        ${_libs}
        )

    SET(_paths)
    FOREACH(_lib
        $(D)/${DEAL_II_LIBRARY_RELDIR}/${CONFIG_LIBRARY_${build}}
        ${DEAL_II_EXTERNAL_LIBRARIES_${build}}
        ${DEAL_II_EXTERNAL_LIBRARIES}
        )
      LIST(APPEND MAKEFILE_LIBS_${build} ${_lib})
      #
      # Extract library paths where possible:
      #
      IF(_lib MATCHES "^(/|\\$\\(D\\))")
        GET_FILENAME_COMPONENT(_path ${_lib} PATH)
        LIST(APPEND _paths ${_path})
      ENDIF()
    ENDFOREACH()

    #
    # And build up an rpath:
    #
    SET(_rpath "")
    LIST(REMOVE_DUPLICATES _paths)
    FOREACH(_path ${_paths})
      SET(_rpath "${_rpath} -Wl,-rpath,${_path}")
    ENDFOREACH()

    TO_STRING(MAKEFILE_LIBS_${build} ${MAKEFILE_LIBS_${build}} ${_rpath})

    #
    # Set up our linker flags:
    #
    SET(MAKEFILE_LDFLAGS_${build}
      "${CMAKE_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${build}} ${DEAL_II_LINKER_FLAGS} ${DEAL_II_LINKER_FLAGS_${build}}"
      )
  ENDFOREACH()

  #
  # Boilerplate: The Make.global_options expects variables to be set to
  # yes, as is common for Makefiles.
  #
  COND_SET_TO_YES(DEAL_II_WITH_THREADS MAKEFILE_enablethreads)
  COND_SET_TO_YES(DEAL_II_WITH_FUNCTIONPARSER MAKEFILE_enableparser)
  COND_SET_TO_YES(BUILD_SHARED_LIBS MAKEFILE_enableshared)

  COND_SET_TO_YES(DEAL_II_WITH_PETSC MAKEFILE_PETSC)
  COND_SET_TO_YES(DEAL_II_WITH_SLEPC MAKEFILE_SLEPC)
  COND_SET_TO_YES(DEAL_II_WITH_TRILINOS MAKEFILE_TRILINOS)
  COND_SET_TO_YES(DEAL_II_WITH_HDF5 MAKEFILE_HDF5)
  COND_SET_TO_YES(DEAL_II_WITH_BLAS MAKEFILE_BLAS)
  COND_SET_TO_YES(DEAL_II_WITH_MUMPS MAKEFILE_MUMPS)
  COND_SET_TO_YES(DEAL_II_WITH_LAPACK MAKEFILE_LAPACK)
  COND_SET_TO_YES(DEAL_II_WITH_ARPACK MAKEFILE_ARPACK)
  COND_SET_TO_YES(DEAL_II_WITH_METIS MAKEFILE_METIS)
  COND_SET_TO_YES(DEAL_II_WITH_UMFPACK MAKEFILE_UMFPACK)
  COND_SET_TO_YES(DEAL_II_WITH_P4EST MAKEFILE_P4EST)
  COND_SET_TO_YES(DEAL_II_WITH_MPI MAKEFILE_MPI)

  #
  # For binary dir:
  #
  SET(MAKEFILE_D_BINARY ${CMAKE_BINARY_DIR})
  TO_STRING_AND_ADD_PREFIX(MAKEFILE_INCLUDE_DIRS_BINARY "${CMAKE_INCLUDE_FLAG_CXX}"
    $(D)/install
    $(D)/install/deal.II
    ${CMAKE_SOURCE_DIR}/include/
    ${CMAKE_SOURCE_DIR}/include/deal.II
    ${deal_ii_source_includes}
    ${DEAL_II_USER_INCLUDE_DIRS}
    )
  #
  # For installation:
  #
  SET(MAKEFILE_D_INSTALLATION ${CMAKE_INSTALL_PREFIX})
  TO_STRING_AND_ADD_PREFIX(MAKEFILE_INCLUDE_DIRS_INSTALLATION "${CMAKE_INCLUDE_FLAG_CXX}"
    $(D)/${DEAL_II_INCLUDE_RELDIR}
    $(D)/${DEAL_II_INCLUDE_RELDIR}/deal.II
    $(D)/${DEAL_II_INCLUDE_RELDIR}/deal.II/bundled
    ${DEAL_II_USER_INCLUDE_DIRS}
    )
  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/Make.global_options.in
    ${CMAKE_BINARY_DIR}/${DEAL_II_COMMON_RELDIR}/Make.global_options
    )

  INSTALL(FILES
    ${CMAKE_BINARY_DIR}/${DEAL_II_COMMON_RELDIR}/template-arguments
    ${CMAKE_BINARY_DIR}/${DEAL_II_COMMON_RELDIR}/Make.global_options
    DESTINATION ${DEAL_II_COMMON_RELDIR}
    COMPONENT compat_files
    )
ENDIF()

#
# Finally, add a target to create the "binary" file in
# ${DEAL_II_PROJECT_CONFIG_RELDIR} and add it to the "all" target:
#
ADD_CUSTOM_TARGET(setup_build_dir ALL
  COMMAND ${CMAKE_COMMAND} -E touch
    ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/binary
  COMMAND ${CMAKE_COMMAND} -E touch
    ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  COMMENT "Update build directory"
  )
FOREACH(_build ${DEAL_II_BUILD_TYPES})
  ADD_DEPENDENCIES(setup_build_dir ${DEAL_II_BASE_NAME}${DEAL_II_${_build}_SUFFIX})
ENDFOREACH()

#
# And a script to remove it upon installation from the install prefix:
# This is necessary if somebody wants to install into the build directory
# (yes this is a valid use case...).
#
INSTALL(CODE
  "
  FILE(REMOVE ${CMAKE_INSTALL_PREFIX}/${DEAL_II_PROJECT_CONFIG_RELDIR}/binary)
  EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E touch
    ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
    )
  "
  COMPONENT library
  )
