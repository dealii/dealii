## ---------------------------------------------------------------------
##
## Copyright (C) 2014 - 2022 by the deal.II authors
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


########################################################################
#                                                                      #
#                Query for git repository information:                 #
#                                                                      #
########################################################################

deal_ii_query_git_information("DEAL_II")

file(WRITE ${CMAKE_BINARY_DIR}/revision.log
"###
#
#  Git information:
#        Branch:    ${DEAL_II_GIT_BRANCH}
#        Revision:  ${DEAL_II_GIT_REVISION}
#        Timestamp: ${DEAL_II_GIT_TIMESTAMP}
#
###"
  )


########################################################################
#                                                                      #
#              Write a nice configuration summary to file:             #
#                                                                      #
########################################################################

set(_log_detailed "${CMAKE_BINARY_DIR}/detailed.log")
set(_log_summary  "${CMAKE_BINARY_DIR}/summary.log")
file(REMOVE ${_log_detailed} ${_log_summary})

macro(_both)
  # Write to both log files:
  file(APPEND ${_log_detailed} "${ARGN}")
  file(APPEND ${_log_summary} "${ARGN}")
endmacro()

macro(_detailed)
  # Only write to detailed.log:
  file(APPEND ${_log_detailed} "${ARGN}")
endmacro()

macro(_summary)
  # Only write to summary.log:
  file(APPEND ${_log_summary} "${ARGN}")
endmacro()

_both(
"###
#
#  ${DEAL_II_PACKAGE_NAME} configuration:
#        CMAKE_BUILD_TYPE:       ${CMAKE_BUILD_TYPE}
#        BUILD_SHARED_LIBS:      ${BUILD_SHARED_LIBS}
#        CMAKE_INSTALL_PREFIX:   ${CMAKE_INSTALL_PREFIX}
#        CMAKE_SOURCE_DIR:       ${CMAKE_SOURCE_DIR}
"
  )
if("${DEAL_II_GIT_SHORTREV}" STREQUAL "")
  _both("#                                (version ${DEAL_II_PACKAGE_VERSION})\n")
else()
  _both("#                                (version ${DEAL_II_PACKAGE_VERSION}, shortrev ${DEAL_II_GIT_SHORTREV})\n")
endif()
_both(
"#        CMAKE_BINARY_DIR:       ${CMAKE_BINARY_DIR}
#        CMAKE_CXX_COMPILER:     ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} on platform ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}
#                                ${CMAKE_CXX_COMPILER}
"
  )
if(DEAL_II_HAVE_CXX20)
  _both("#        C++ language standard:  C++20\n")
elseif(DEAL_II_HAVE_CXX17)
  _both("#        C++ language standard:  C++17\n")
elseif(DEAL_II_HAVE_CXX14)
  _both("#        C++ language standard:  C++14\n")
endif()

if(CMAKE_C_COMPILER_WORKS)
  _detailed("#        CMAKE_C_COMPILER:       ${CMAKE_C_COMPILER}\n")
endif()
if(CMAKE_Fortran_COMPILER_WORKS)
  _detailed("#        CMAKE_Fortran_COMPILER: ${CMAKE_Fortran_COMPILER}\n")
endif()
_detailed("#        CMAKE_GENERATOR:        ${CMAKE_GENERATOR}\n")

if(CMAKE_CROSSCOMPILING)
  _both(
    "#\n#        CROSSCOMPILING!\n"
    )
endif()

_both("#\n")

_detailed(
"#  Base configuration (prior to feature configuration):
#        DEAL_II_CXX_FLAGS:            ${BASE_CXX_FLAGS}
"
  )
if(CMAKE_BUILD_TYPE MATCHES "Release")
  _detailed("#        DEAL_II_CXX_FLAGS_RELEASE:    ${BASE_CXX_FLAGS_RELEASE}\n")
endif()
if(CMAKE_BUILD_TYPE MATCHES "Debug")
  _detailed("#        DEAL_II_CXX_FLAGS_DEBUG:      ${BASE_CXX_FLAGS_DEBUG}\n")
endif()

_detailed("#        DEAL_II_LINKER_FLAGS:         ${BASE_LINKER_FLAGS}\n")
if(CMAKE_BUILD_TYPE MATCHES "Release")
  _detailed("#        DEAL_II_LINKER_FLAGS_RELEASE: ${BASE_LINKER_FLAGS_RELEASE}\n")
endif()
if(CMAKE_BUILD_TYPE MATCHES "Debug")
  _detailed("#        DEAL_II_LINKER_FLAGS_DEBUG:   ${BASE_LINKER_FLAGS_DEBUG}\n")
endif()

_detailed("#        DEAL_II_DEFINITIONS:          ${BASE_DEFINITIONS}\n")
if(CMAKE_BUILD_TYPE MATCHES "Release")
  _detailed("#        DEAL_II_DEFINITIONS_RELEASE:  ${BASE_DEFINITIONS_RELEASE}\n")
endif()
if(CMAKE_BUILD_TYPE MATCHES "Debug")
  _detailed("#        DEAL_II_DEFINITIONS_DEBUG:    ${BASE_DEFINITIONS_DEBUG}\n")
endif()

_detailed("#        DEAL_II_INCLUDE_DIRS          ${BASE_INCLUDE_DIRS}\n")
_detailed("#        DEAL_II_BUNDLED_INCLUDE_DIRS: ${BASE_BUNDLED_INCLUDE_DIRS}\n")

_detailed("#        DEAL_II_LIBRARIES:            ${BASE_LIBRARIES}\n")
if(CMAKE_BUILD_TYPE MATCHES "Release")
  _detailed("#        DEAL_II_LIBRARIES_RELEASE:    ${BASE_LIBRARIES_RELEASE}\n")
endif()
if(CMAKE_BUILD_TYPE MATCHES "Debug")
  _detailed("#        DEAL_II_LIBRARIES_DEBUG:      ${BASE_LIBRARIES_DEBUG}\n")
endif()
_detailed("#        DEAL_II_VECTORIZATION_WIDTH_IN_BITS: ${DEAL_II_VECTORIZATION_WIDTH_IN_BITS}\n")

if(DEAL_II_HAVE_CXX20)
  _detailed("#        DEAL_II_HAVE_CXX20\n")
elseif(DEAL_II_HAVE_CXX17)
  _detailed("#        DEAL_II_HAVE_CXX17\n")
elseif(DEAL_II_HAVE_CXX14)
  _detailed("#        DEAL_II_HAVE_CXX14\n")
endif()


_detailed("#\n")

_both("#  Configured Features (")
if(DEFINED DEAL_II_ALLOW_BUNDLED)
  _both("DEAL_II_ALLOW_BUNDLED = ${DEAL_II_ALLOW_BUNDLED}, ")
endif()
if(DEAL_II_FORCE_AUTODETECTION)
  _both("!!! DEAL_II_FORCE_AUTODETECTION=ON !!!, ")
endif()
_both("DEAL_II_ALLOW_AUTODETECTION = ${DEAL_II_ALLOW_AUTODETECTION}):\n")


set(_deal_ii_features_sorted ${DEAL_II_FEATURES})
list(SORT _deal_ii_features_sorted)
foreach(_feature ${_deal_ii_features_sorted})
  set(_var DEAL_II_WITH_${_feature})

  if(${${_var}})
    #
    # The feature is enabled:
    #
    if(FEATURE_${_feature}_EXTERNAL_CONFIGURED)
      _both("#        ${_var} set up with external dependencies\n")
    elseif(DEAL_II_FEATURE_${_feature}_BUNDLED_CONFIGURED)
      if(DEAL_II_FORCE_BUNDLED_${_feature})
        _both("#        ${_var} set up with bundled packages (forced)\n")
      else()
        _both("#        ${_var} set up with bundled packages\n")
      endif()
    else()
     _both("#        ${_var} = ${${_var}}\n")
    endif()

    #
    # Print out version number:
    #
    if(DEFINED ${_feature}_VERSION)
      _detailed("#            ${_feature}_VERSION = ${${_feature}_VERSION}\n")
    endif()

    #
    # Special version numbers:
    #
    if(_feature MATCHES "THREADS" AND DEFINED TBB_VERSION)
      _detailed("#            TBB_VERSION = ${TBB_VERSION}\n")
    endif()

    if(_feature MATCHES "KOKKOS")
      _detailed("#            KOKKOS_BACKENDS = ${Kokkos_DEVICES}\n")
      _detailed("#            KOKKOS_ARCHITECTURES = ${Kokkos_ARCH}\n")
    endif()

    #
    # Print out ${_feature}_DIR:
    #
    if(NOT "${${_feature}_DIR}" STREQUAL "")
      _detailed("#            ${_feature}_DIR = ${${_feature}_DIR}\n")
    endif()

    if(NOT "${${_feature}_SPLIT_CONFIGURATION}" STREQUAL "")
      _detailed("#            ${_feature}_SPLIT_CONFIGURATION = ${${_feature}_SPLIT_CONFIGURATION}\n")
    endif()

    #
    # Print the feature configuration:
    #
    foreach(_var2
      C_COMPILER CXX_COMPILER Fortran_COMPILER WITH_64BIT_BLAS_INDICES
      ${DEAL_II_STRING_SUFFIXES} ${DEAL_II_LIST_SUFFIXES}
      )
      if(DEFINED ${_feature}_${_var2})
        _detailed("#            ${_feature}_${_var2} = ${${_feature}_${_var2}}\n")
      endif()
    endforeach()
  else()
    # FEATURE is disabled
    _both("#      ( ${_var} = ${${_var}} )\n")
  endif()
endforeach()

_both(
  "#\n#  Component configuration:\n"
  )

foreach(_component ${DEAL_II_COMPONENTS})
  set(_var DEAL_II_COMPONENT_${_component})
  if(${${_var}})
    _both("#        ${_var}\n")
  else()
    _both("#      ( ${_var} = ${${_var}} )\n")
  endif()
endforeach()

_summary(
"#\n#  Detailed information (compiler flags, feature configuration) can be found in detailed.log
#\n#  Run  $ "
  )
if(CMAKE_GENERATOR MATCHES "Ninja")
  _summary("ninja ")
else()
_summary("make ")
endif()
_summary("info  to print a help message with a list of top level targets\n")

_both("#\n###")
