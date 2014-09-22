## ---------------------------------------------------------------------
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
# We support two configurations out of which deal.II can be used - directly
# from the build directory or after installation. So we have to prepare
# two distinct set ups.
#

MESSAGE(STATUS "Setting up project configuration")

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
# Do not force --as-needed for executables on user side:
#

STRIP_FLAG(DEAL_II_LINKER_FLAGS "-Wl,--as-needed")

#
# Strip -Wno-deprecated from DEAL_II_CXX_FLAGS so that deprecation warnings
# are actually shown for user code:
#

STRIP_FLAG(DEAL_II_CXX_FLAGS "-Wno-deprecated")
STRIP_FLAG(DEAL_II_CXX_FLAGS "-Wno-deprecated-declarations")

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

  IF(CMAKE_SYSTEM_NAME MATCHES "CYGWIN" OR CMAKE_SYSTEM_NAME MATCHES "Windows")
    SET(CONFIG_LIBRARIES_${_build}
      \${DEAL_II_PATH}/\${DEAL_II_EXECUTABLE_RELDIR}/${CONFIG_LIBRARY_${_build}}
      )
  ELSE()
    SET(CONFIG_LIBRARIES_${_build}
      \${DEAL_II_PATH}/\${DEAL_II_LIBRARY_RELDIR}/${CONFIG_LIBRARY_${_build}}
      )
  ENDIF()

  LIST(APPEND CONFIG_LIBRARIES_${_build}
    ${DEAL_II_LIBRARIES_${_build}}
    ${DEAL_II_LIBRARIES}
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

SET(CONFIG_BUILD_DIR TRUE)
SET(CONFIG_INCLUDE_DIRS
  \${DEAL_II_PATH}/include
  \${DEAL_II_PATH}/include/deal.II
  ${CMAKE_SOURCE_DIR}/include/
  ${CMAKE_SOURCE_DIR}/include/deal.II
  ${DEAL_II_BUNDLED_INCLUDE_DIRS}
  ${DEAL_II_USER_INCLUDE_DIRS}
  )
CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  @ONLY
  )
CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/ConfigVersion.cmake.in
  ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}ConfigVersion.cmake
  @ONLY
  )

#
# For installation:
#

SET(CONFIG_BUILD_DIR FALSE)
SET(CONFIG_INCLUDE_DIRS
  \${DEAL_II_PATH}/\${DEAL_II_INCLUDE_RELDIR}
  \${DEAL_II_PATH}/\${DEAL_II_INCLUDE_RELDIR}/deal.II
  \${DEAL_II_PATH}/\${DEAL_II_INCLUDE_RELDIR}/deal.II/bundled
  ${DEAL_II_USER_INCLUDE_DIRS}
  )
CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  @ONLY
  )
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

#
# Append feature configuration to both configuration files:
#

SET(_files
  ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  )
FOREACH(_file ${_files})
  FILE(APPEND ${_file}
    "\n\n#\n# Feature configuration:\n#\n\n"
    )
ENDFOREACH()
GET_CMAKE_PROPERTY(res VARIABLES)
FOREACH(var ${res})
  IF(var MATCHES "DEAL_II_WITH")
    FOREACH(_file ${_files})
      FILE(APPEND ${_file}
        "SET(${var} ${${var}})\n"
        )
    ENDFOREACH()
  ENDIF()
ENDFOREACH()


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
    # not a full path starting with "/" or that happen to be a "framework"
    # - we won't deal with that in Make.global_options.
    #

    SET(_libs)
    FOREACH(_lib
        ${DEAL_II_LIBRARIES_${build}}
        ${DEAL_II_LIBRARIES}
        )
      IF(_lib MATCHES "^\\s*\\/" AND NOT _lib MATCHES "\\.framework$")
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
        ${DEAL_II_LIBRARIES_${build}}
        ${DEAL_II_LIBRARIES}
        )

      #
      # Strip leading and trailing whitespace
      #
      STRING(STRIP "${_lib}" _lib)

      IF(_lib MATCHES "^\\/.*\\.framework$")
        #
        # We have a MacOSX framework
        #
        # TODO: What do?
        #

      ELSEIF(_lib MATCHES "^(\\/|\\$\\(D\\))")
        #
        # Full path: /path/to/libfoo.(so|dylib|etc):
        #

        # Extract library paths:
        GET_FILENAME_COMPONENT(_path ${_lib} PATH)
        LIST(APPEND _paths ${_path})

        LIST(APPEND MAKEFILE_LIBS_${build} ${_lib})

      ELSEIF(_lib MATCHES "^-l")
        #
        # Library is of the form "-lfoo":
        #
        LIST(APPEND MAKEFILE_LIBS_${build} ${_lib})

      ELSEIF(NOT _lib MATCHES "^\\s*$" AND NOT _lib MATCHES "\\.so(\\.[0-9]+)*$")
        #
        # Well in this case we just assume that we have to append a "-l"
        #

        LIST(APPEND MAKEFILE_LIBS_${build} "-l${_lib}")

      ELSE()
        #
        # For now, ignore the rest.
        #
        # TODO: What do?
        #
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
      "${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${build}} ${DEAL_II_LINKER_FLAGS} ${DEAL_II_LINKER_FLAGS_${build}}"
      )
  ENDFOREACH()

  #
  # Boilerplate: The Make.global_options expects variables to be set to
  # yes, as is common for Makefiles.
  #
  COND_SET_TO_YES(DEAL_II_WITH_THREADS MAKEFILE_enablethreads)
  COND_SET_TO_YES(DEAL_II_WITH_FUNCTIONPARSER MAKEFILE_enableparser)
  COND_SET_TO_YES(BUILD_SHARED_LIBS MAKEFILE_enableshared)

  COND_SET_TO_YES(DEAL_II_WITH_ARPACK MAKEFILE_ARPACK)
  COND_SET_TO_YES(DEAL_II_WITH_BLAS MAKEFILE_BLAS)
  COND_SET_TO_YES(DEAL_II_WITH_HDF5 MAKEFILE_HDF5)
  COND_SET_TO_YES(DEAL_II_WITH_LAPACK MAKEFILE_LAPACK)
  COND_SET_TO_YES(DEAL_II_WITH_METIS MAKEFILE_METIS)
  COND_SET_TO_YES(DEAL_II_WITH_MPI MAKEFILE_MPI)
  COND_SET_TO_YES(DEAL_II_WITH_MUMPS MAKEFILE_MUMPS)
  COND_SET_TO_YES(DEAL_II_WITH_P4EST MAKEFILE_P4EST)
  COND_SET_TO_YES(DEAL_II_WITH_PETSC MAKEFILE_PETSC)
  COND_SET_TO_YES(DEAL_II_WITH_SLEPC MAKEFILE_SLEPC)
  COND_SET_TO_YES(DEAL_II_WITH_TRILINOS MAKEFILE_TRILINOS)
  COND_SET_TO_YES(DEAL_II_WITH_UMFPACK MAKEFILE_UMFPACK)

  #
  # For binary dir:
  #
  SET(MAKEFILE_D ${CMAKE_BINARY_DIR})
  TO_STRING_AND_ADD_PREFIX(MAKEFILE_INCLUDE_DIRS "${CMAKE_INCLUDE_FLAG_CXX}"
    $(D)/install
    $(D)/install/deal.II
    ${CMAKE_SOURCE_DIR}/include/
    ${CMAKE_SOURCE_DIR}/include/deal.II
    ${DEAL_II_BUNDLED_INCLUDE_DIRS}
    ${DEAL_II_USER_INCLUDE_DIRS}
    )
  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/Make.global_options.in
    ${CMAKE_BINARY_DIR}/${DEAL_II_COMMON_RELDIR}/Make.global_options
    )
  #
  # For installation:
  #
  SET(MAKEFILE_D ${CMAKE_INSTALL_PREFIX})
  TO_STRING_AND_ADD_PREFIX(MAKEFILE_INCLUDE_DIRS "${CMAKE_INCLUDE_FLAG_CXX}"
    $(D)/${DEAL_II_INCLUDE_RELDIR}
    $(D)/${DEAL_II_INCLUDE_RELDIR}/deal.II
    $(D)/${DEAL_II_INCLUDE_RELDIR}/deal.II/bundled
    ${DEAL_II_USER_INCLUDE_DIRS}
    )
  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/Make.global_options.in
    ${CMAKE_CURRENT_BINARY_DIR}/Make.global_options
    )
  INSTALL(FILES
    ${CMAKE_BINARY_DIR}/${DEAL_II_COMMON_RELDIR}/template-arguments
    ${CMAKE_CURRENT_BINARY_DIR}/Make.global_options
    DESTINATION ${DEAL_II_COMMON_RELDIR}
    COMPONENT compat_files
    )
ENDIF()

MESSAGE(STATUS "Setting up project configuration - Done")
