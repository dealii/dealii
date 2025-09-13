## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2014 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------


########################################################################
#                                                                      #
#              Write a nice configuration summary to file:             #
#                                                                      #
########################################################################

set(_log_detailed "${CMAKE_BINARY_DIR}/detailed.log")
set(_log_summary  "${CMAKE_BINARY_DIR}/summary.log")
file(REMOVE ${_log_detailed} ${_log_summary})

function(_both)
  # Write to both log files:
  file(APPEND ${_log_detailed} "${ARGN}")
  file(APPEND ${_log_summary} "${ARGN}")
endfunction()

function(_detailed)
  # Only write to detailed.log:
  file(APPEND ${_log_detailed} "${ARGN}")
endfunction()

function(_summary)
  # Only write to summary.log:
  file(APPEND ${_log_summary} "${ARGN}")
endfunction()

function(_print_target _target)
  print_target_properties(${_target} _messages)
  foreach(_message ${_messages})
    _detailed("#    ${_message}\n")
  endforeach()
endfunction()

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
if(CMAKE_C_COMPILER_WORKS)
  _detailed("#        CMAKE_C_COMPILER:       ${CMAKE_C_COMPILER}\n")
endif()
if(CMAKE_Fortran_COMPILER_WORKS)
  _detailed("#        CMAKE_Fortran_COMPILER: ${CMAKE_Fortran_COMPILER}\n")
endif()
_detailed("#        CMAKE_GENERATOR:        ${CMAKE_GENERATOR}\n")

if(DEAL_II_HAVE_CXX23)
  _both("#        C++ language standard:  C++23\n")
elseif(DEAL_II_HAVE_CXX20)
  _both("#        C++ language standard:  C++20\n")
elseif(DEAL_II_HAVE_CXX17)
  _both("#        C++ language standard:  C++17\n")
endif()

if(DEAL_II_HAVE_CXX23 OR DEAL_II_HAVE_CXX20)
  if(DEAL_II_WITH_CXX20_MODULE)
    _both("#        Building C++20 module:  ON\n")
  else()
    _both("#        Building C++20 module:  OFF\n")
  endif()
endif()

_both("#        Vectorization level:    ${DEAL_II_VECTORIZATION_WIDTH_IN_BITS} bit")
set(_instructions)
if(DEAL_II_HAVE_SSE2)
  list(APPEND _instructions "sse2")
endif()
if(DEAL_II_HAVE_AVX)
  list(APPEND _instructions "avx2")
endif()
if(DEAL_II_HAVE_AVX512)
  list(APPEND _instructions "avx512*")
endif()
if(DEAL_II_HAVE_ALTIVEC)
  list(APPEND _instructions "altivec")
endif()
if(DEAL_II_HAVE_ARM_NEON)
  list(APPEND _instructions "arm_neon")
endif()
if(NOT "${_instructions}" STREQUAL "")
  to_string(_string ${_instructions})
  _both(" (${_string})\n")
else()
  _both("\n")
endif()

if(CMAKE_CROSSCOMPILING)
  _both(
    "#\n#        CROSSCOMPILING!\n"
    )
endif()

_both("#\n")

_detailed("#  Exported compiler and linker flags:\n")
_detailed("#        DEAL_II_CXX_FLAGS:            ${DEAL_II_CXX_FLAGS}\n")
if(CMAKE_BUILD_TYPE MATCHES "Release")
  _detailed("#        DEAL_II_CXX_FLAGS_RELEASE:    ${DEAL_II_CXX_FLAGS_RELEASE}\n")
endif()
if(CMAKE_BUILD_TYPE MATCHES "Debug")
  _detailed("#        DEAL_II_CXX_FLAGS_DEBUG:      ${DEAL_II_CXX_FLAGS_DEBUG}\n")
endif()
_detailed("#        DEAL_II_WARNING_FLAGS:        ${DEAL_II_WARNING_FLAGS}\n")

_detailed("#        DEAL_II_LINKER_FLAGS:         ${DEAL_II_LINKER_FLAGS}\n")
if(CMAKE_BUILD_TYPE MATCHES "Release")
  _detailed("#        DEAL_II_LINKER_FLAGS_RELEASE: ${DEAL_II_LINKER_FLAGS_RELEASE}\n")
endif()
if(CMAKE_BUILD_TYPE MATCHES "Debug")
  _detailed("#        DEAL_II_LINKER_FLAGS_DEBUG:   ${DEAL_II_LINKER_FLAGS_DEBUG}\n")
endif()

_detailed("#  Library targets:\n")

_detailed("#\n")
_print_target(${DEAL_II_TARGET_NAME})

foreach(_build ${DEAL_II_BUILD_TYPES})
  string(TOLOWER "${_build}" _build_lowercase)
  _detailed("#\n")
  _print_target(${DEAL_II_TARGET_NAME}_${_build_lowercase})
endforeach()
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
    # Print some non interface target related configuration values.
    # Specifically, go through all combinations of
    #   ${FEATURE}_EXECUTABLE
    #   ${FEATURE}_VERSION
    #   ...
    # and if a variable of that name is found, prints its value.
    #

    foreach(_conf_var
      EXECUTABLE VERSION DIR C_COMPILER CXX_COMPILER Fortran_COMPILER
      WITH_64BIT_BLAS_INDICES
      )
      if(NOT "${${_feature}_${_conf_var}}" STREQUAL "")
        _detailed("#            ${_feature}_${_conf_var} = ${${_feature}_${_conf_var}}\n")
      endif()
    endforeach()

    #
    # Special information:
    #

    if(_feature MATCHES "MPI")
      _detailed("#            MPI_CXX_VERSION = ${MPI_CXX_VERSION}\n")
    endif()

    if(_feature MATCHES "KOKKOS")
      _detailed("#            KOKKOS_BACKENDS = ${Kokkos_DEVICES}\n")
      _detailed("#            KOKKOS_ARCHITECTURES = ${Kokkos_ARCH}\n")
    endif()
  else()
    # FEATURE is disabled
    _both("#      ( ${_var} = ${${_var}} )\n")
  endif()
endforeach()

_detailed("#\n#  Interface targets:\n")

foreach(_feature ${_deal_ii_features_sorted})
  if(${DEAL_II_WITH_${_feature}})
    string(TOLOWER ${_feature} _feature_lowercase)
    foreach(_target ${DEAL_II_BUILD_TYPES}
      interface_${_feature_lowercase} interface_${_feature_lowercase}_debug interface_${_feature_lowercase}_release
      bundled_${_feature_lowercase} bundled_${_feature_lowercase}_debug bundled_${_feature_lowercase}_release
      )
      if(TARGET ${_target})
        _detailed("#\n")
        _print_target(${_target})
      endif()
    endforeach()
  endif()
endforeach()

_both("#\n#  Component configuration:\n")

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
