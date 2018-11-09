## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2017 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

#
# This file sets up the project configuration consisting of
#
#   deal.IIConfig.cmake
#   deal.IIVersionConfig.cmake
#
# We support two configurations out of which deal.II can be used - directly
# from the build directory or after installation. So we have to prepare
# two distinct setups.
#

MESSAGE(STATUS "Setting up project configuration")

#
# Configure the template-arguments file
#
CONFIGURE_FILE( # for binary dir:
  ${CMAKE_CURRENT_SOURCE_DIR}/template-arguments.in
  ${CMAKE_BINARY_DIR}/${DEAL_II_SHARE_RELDIR}/template-arguments
  )


########################################################################
#                                                                      #
#                   Setup and install cmake macros:                    #
#                                                                      #
########################################################################

SET(_macros
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_initialize_cached_variables.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_invoke_autopilot.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_setup_target.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_query_git_information.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_add_test.cmake
  ${CMAKE_SOURCE_DIR}/cmake/macros/macro_deal_ii_pickup_tests.cmake
  )
FILE(COPY ${_macros}
  DESTINATION ${CMAKE_BINARY_DIR}/${DEAL_II_SHARE_RELDIR}/macros
  )
INSTALL(FILES ${_macros}
  DESTINATION ${DEAL_II_SHARE_RELDIR}/macros
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
# Strip -Wno-deprecated-declarations from DEAL_II_CXX_FLAGS so that
# deprecation warnings are actually shown for user code:
#

STRIP_FLAG(DEAL_II_CXX_FLAGS "-Wno-deprecated-declarations")

#
# Strip -Wsuggest-override. Enabling this warning is a stylistic choice
# that might create a substantial amount of warnings in user codebases.
# Let's be nice and disable the flag for now.
#

STRIP_FLAG(DEAL_II_CXX_FLAGS "-Wsuggest-override")

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
  # Build up library name depending on link type and platform:
  #

  IF(BUILD_SHARED_LIBS)
    SET(_type "SHARED")
  ELSE()
    SET(_type "STATIC")
  ENDIF()

  SET(_name "${CMAKE_${_type}_LIBRARY_PREFIX}${DEAL_II_BASE_NAME}${DEAL_II_${_build}_SUFFIX}${CMAKE_${_type}_LIBRARY_SUFFIX}")
  SET(CONFIG_LIBRARIES_${_build}
    "\${DEAL_II_PATH}/${DEAL_II_LIBRARY_RELDIR}/${_name}"
    ${DEAL_II_LIBRARIES_${_build}}
    ${DEAL_II_LIBRARIES}
    )
  LIST(APPEND CONFIG_LIBRARIES ${_keyword} \${DEAL_II_LIBRARIES_${_build}})

  SET(CONFIG_TARGET_${_build} ${DEAL_II_BASE_NAME}${DEAL_II_${_build}_SUFFIX})
  LIST(APPEND CONFIG_TARGET ${_keyword} \${DEAL_II_TARGET_${_build}})
ENDFOREACH()

#
# For binary dir:
#

SET(CONFIG_BUILD_DIR TRUE)
SET(CONFIG_INCLUDE_DIRS
  \${DEAL_II_PATH}/include
  ${CMAKE_SOURCE_DIR}/include/
  ${DEAL_II_BUNDLED_INCLUDE_DIRS}
  ${DEAL_II_USER_INCLUDE_DIRS}
  )

SET(CONFIG_INCLUDE_DIRS_BUILD_PC ${CONFIG_INCLUDE_DIRS}) # for pkgconfig

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
  ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}FeatureConfig.cmake
  DESTINATION ${DEAL_II_PROJECT_CONFIG_RELDIR}
  COMPONENT library
  )

#
# Append feature configuration to all configuration files:
#

SET(_files
  ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}Config.cmake
  )
FOREACH(_file ${_files})
  FILE(APPEND ${_file} "\n\n#\n# Feature configuration:\n#\n\n")
  FILE(APPEND ${_file} "SET(DEAL_II_WITH_CXX11 ON)\n")
ENDFOREACH()

GET_CMAKE_PROPERTY(_res VARIABLES)
FOREACH(_var ${_res})
  IF(_var MATCHES "DEAL_II_.*WITH")
    LIST(APPEND _additional_config_variables ${_var})
  ENDIF()
ENDFOREACH()

SET(_deal_ii_features_sorted ${DEAL_II_FEATURES})
LIST(SORT _deal_ii_features_sorted)
FOREACH(_name ${_deal_ii_features_sorted})
  SET(_var DEAL_II_WITH_${_name})

  FOREACH(_file ${_files})
    FILE(APPEND ${_file} "SET(${_var} ${${_var}})\n")
    #
    # Do not pollute deal.IIConfig.cmake with package details of
    # unconfigured features.
    #
    IF(${_var})
      IF(NOT "${${_name}_VERSION}" STREQUAL "")
        FILE(APPEND ${_file}
          "SET(DEAL_II_${_name}_VERSION \"${${_name}_VERSION}\")\n"
          )
      ENDIF()
      FOREACH(_additional ${_additional_config_variables})
        IF(_additional MATCHES "DEAL_II_${_name}_WITH")
          FILE(APPEND ${_file} "SET(${_additional} ${${_additional}})\n")
        ENDIF()
      ENDFOREACH()
    ENDIF()
  ENDFOREACH()
ENDFOREACH()

#####################################################################
# pkgconfig files
#
# Set up the pkgconfig configuration files consisting of
#
#   deal.II_base       # common part to all the rest
#   deal.II_debug      # generated when buildtype is debug
#   deal.II_release    # generated when buildtype is release
#
# We support two configurations out of which deal.II can be used - directly
# from the build directory or after installation. So we have to prepare
# two distinct setups.
#
# pkgconfig looks for *.pc files in an environmental variable called
# PKG_CONFIG_PATH. So, to use the library in the build directory issue
# export PKG_CONFIG_PATH=/path/to/BUILD_DIR/lib/pkgconfig:$PKG_CONFIG_PATH
#
# To use the library in the installed location
# export PKG_CONFIG_PATH=/path/to/INSTALL_DIR/lib/pkgconfig:$PKG_CONFIG_PATH
#

#
# For binary dir (AKA build dir):
#

set(DEAL_II_PATH_PC ${CMAKE_BINARY_DIR})

set (CONFIG_INCLUDE_DIRS_PC "")
foreach(_dir ${CONFIG_INCLUDE_DIRS})
  string(CONCAT CONFIG_INCLUDE_DIRS_PC "${CONFIG_INCLUDE_DIRS_PC} -I${_dir} "  )
endforeach()

# populate the variables
# CONFIG_LIBRARIES_DEBUG_PC 
# CONFIG_LIBRARIES_RELEASE_PC

foreach (_build ${DEAL_II_BUILD_TYPES})
  if(_build MATCHES DEBUG)
    set(_build "debug")
  else()
    set(_build "release")
  endif()

  string(TOUPPER ${_build} _B)
  set(PC_LIBRARIES ${CONFIG_LIBRARIES_${_B}})
  set(CONFIG_LIBRARIES_PC "")
  foreach(_s ${PC_LIBRARIES})
    # check if _s contains lib
    # if _s is somethink like /path/to/trilinos/lib/libteuchos.so
    # we have  to add /path/to/trilinos/lib to the runpath
    if(${_s} MATCHES "lib")
      string(CONCAT CONFIG_LIBRARIES_PC " ${CONFIG_LIBRARIES_PC} ${_s} ")

      # work for the runpath
      string(FIND ${_s} "/" _i REVERSE)
      string(SUBSTRING ${_s} 0 ${_i} _p)
      list(APPEND PC_RPATH ${_p})
    else()
      # here _s is something like "m" for the libm.so so we add the -l flag
      string(CONCAT CONFIG_LIBRARIES_PC " ${CONFIG_LIBRARIES_PC} -l${_s} "  )
    endif()
  endforeach()

  set(CONFIG_LIBRARIES_${_B}_PC ${CONFIG_LIBRARIES_PC})
  
  list(REMOVE_DUPLICATES PC_RPATH)
  foreach(_r ${PC_RPATH})
    string(CONCAT CONFIG_LIBRARIES_${_B}_PC " ${CONFIG_LIBRARIES_${_B}_PC} -Wl,-rpath,${_r} "  )
  endforeach()

endforeach()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/deal.II_base.pc.in
  ${CMAKE_BINARY_DIR}/${DEAL_II_LIBRARY_RELDIR}/pkgconfig/${DEAL_II_PROJECT_CONFIG_NAME}_base.pc
  @ONLY
  )

# generate the requested .pc files
foreach (_build ${DEAL_II_BUILD_TYPES})
  if(_build MATCHES DEBUG)
    set(_build "debug")
  else()
    set(_build "release")
  endif()

  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/deal.II_${_build}.pc.in
    ${CMAKE_BINARY_DIR}/${DEAL_II_LIBRARY_RELDIR}/pkgconfig/${DEAL_II_PROJECT_CONFIG_NAME}_${_build}.pc
    @ONLY
    )
endforeach()

#############################################################################
#
# For installation:
#

set(DEAL_II_PATH_PC ${CMAKE_INSTALL_PREFIX})

set (CONFIG_INCLUDE_DIRS_PC "")
foreach(_dir ${CONFIG_INCLUDE_DIRS})
  string(CONCAT CONFIG_INCLUDE_DIRS_PC "${CONFIG_INCLUDE_DIRS_PC} -I${_dir} "  )
endforeach()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/deal.II_base.pc.in
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}_base.pc
  @ONLY
  )

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}_base.pc
  DESTINATION ${DEAL_II_LIBRARY_RELDIR}/pkgconfig
  COMPONENT library
  )


# generate the requested .pc files
foreach (_build ${DEAL_II_BUILD_TYPES})
  if(NOT "${CMAKE_BUILD_TYPE}" STREQUAL "DebugRelease")
    set(_build "general")
  else()
    if(_build MATCHES DEBUG)
      set(_build "debug")
    else()
      set(_build "release")
    endif()
  endif()

  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/deal.II_${_build}.pc.in
    ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}_${_build}.pc
    @ONLY
    )
  install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_NAME}_${_build}.pc
    DESTINATION ${DEAL_II_LIBRARY_RELDIR}/pkgconfig
    COMPONENT library
    )
endforeach()








# end of pkgconfig files
####################################################################

MESSAGE(STATUS "Setting up project configuration - Done")
