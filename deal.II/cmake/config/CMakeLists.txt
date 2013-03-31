#####
##
## Copyright (C) 2012, 2013 by the deal.II authors
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


#
# Configure the template-arguments file
#
CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/template-arguments.in
  ${CMAKE_CURRENT_BINARY_DIR}/template-arguments
  )


###########################################################################
#                                                                         #
#               Setup and install the convenience macros:                 #
#                                                                         #
###########################################################################

INSTALL(FILES
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_initialize_cached_variables.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_invoke_autopilot.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_setup_target.cmake
  DESTINATION ${DEAL_II_CMAKE_MACROS_RELDIR}
  COMPONENT library
  )


###########################################################################
#                                                                         #
#         Configure and install the cmake project configuration:          #
#                                                                         #
###########################################################################

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
    ${DEAL_II_EXTERNAL_LIBRARIES}
    ${DEAL_II_EXTERNAL_LIBRARIES_${_build}}
    )
  LIST(APPEND CONFIG_LIBRARIES
    ${_keyword}
    ${CONFIG_LIBRARIES_${_build}}
    )

  SET(CONFIG_TARGET_${_build}
    ${DEAL_II_BASE_NAME}${DEAL_II_${_build}_SUFFIX}
    )
  LIST(APPEND CONFIG_TARGET
    ${_keyword}
    ${DEAL_II_BASE_NAME}${DEAL_II_${_build}_SUFFIX}
    )
ENDFOREACH()

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

#
# Append feature configuration:
#

SET(_file ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake)
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
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}ConfigVersion.cmake
  @ONLY
  )

INSTALL(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}ConfigVersion.cmake
  DESTINATION ${DEAL_II_PROJECT_CONFIG_RELDIR}
  COMPONENT library
  )


###########################################################################
#                                                                         #
#          Configure and install the old Makefile configuration:          #
#                                                                         #
###########################################################################


IF(DEAL_II_COMPONENT_COMPAT_FILES)
  #
  # Transform some cmake lists into a string that the old Makefile
  # mechanism actually understands:

  TO_STRING_AND_ADD_PREFIX(MAKEFILE_INCLUDE_DIRS "${CMAKE_INCLUDE_FLAG_CXX}"
    $(D)/${DEAL_II_INCLUDE_RELDIR}
    $(D)/${DEAL_II_INCLUDE_RELDIR}/deal.II
    $(D)/${DEAL_II_INCLUDE_RELDIR}/deal.II/bundled
    ${DEAL_II_USER_INCLUDE_DIRS}
    ${DEAL_II_INCLUDE_DIRS}
    )

  FOREACH(build ${DEAL_II_BUILD_TYPES})
    TO_STRING_AND_ADD_PREFIX(MAKEFILE_DEFINITIONS_${build}
      "-D"
      ${DEAL_II_USER_DEFINITIONS}
      ${DEAL_II_USER_DEFINITIONS_${build}}
      )

    TO_STRING(MAKEFILE_TARGETS_${build}
        $(D)/${DEAL_II_LIBRARY_RELDIR}/${CONFIG_LIBRARY_${build}}
        ${DEAL_II_EXTERNAL_LIBRARIES}
        ${DEAL_II_EXTERNAL_LIBRARIES_${build}}
        )

    #
    # TODO: Come up with an idea how to get the link line directly from
    # CMake...
    #
    FOREACH(_lib
        $(D)/${DEAL_II_LIBRARY_RELDIR}/${CONFIG_LIBRARY_${build}}
        ${DEAL_II_EXTERNAL_LIBRARIES}
        ${DEAL_II_EXTERNAL_LIBRARIES_${build}}
        )
      LIST(APPEND MAKEFILE_LIBS_${build} ${_lib})
      #
      # Add an rpath directive for each library:
      #
      GET_FILENAME_COMPONENT(_path ${_lib} PATH)
      LIST(APPEND MAKEFILE_LIBS_${build} "-Wl,-rpath,${_path}")
    ENDFOREACH()

    TO_STRING(MAKEFILE_LIBS_${build} ${MAKEFILE_LIBS_${build}})

    #
    # Set up our linker flags:
    #
    SET(MAKEFILE_LDFLAGS_${build}
      "${CMAKE_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${build}} ${CMAKE_SHARED_LINKER_FLAGS} ${DEAL_II_SHARED_LINKER_FLAGS_${build}}"
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

  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/Make.global_options.in
    ${CMAKE_CURRENT_BINARY_DIR}/Make.global_options
    )

  INSTALL(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/template-arguments
    ${CMAKE_CURRENT_BINARY_DIR}/Make.global_options
    DESTINATION ${DEAL_II_COMMON_RELDIR}
    COMPONENT compat_files
    )
ENDIF()

