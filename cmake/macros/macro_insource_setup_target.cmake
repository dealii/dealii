## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2013 - 2024 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

#
# This file provides an "insource" version of the DEAL_II_SETUP_TARGET macro.
#
# Usage:
#       insource_setup_target(target build)
#
# This appends necessary include directories, linker flags, compile
# definitions and the deal.II library link interface to the given target.
#

function(insource_setup_target _target _build)
  string(TOLOWER ${_build} _build_lowercase)

  set_target_properties(${_target} PROPERTIES
    LINKER_LANGUAGE "CXX"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    )

  target_compile_flags(${_target} PRIVATE
    "${DEAL_II_WARNING_FLAGS} ${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${_build}}"
    )

  get_property(_type TARGET ${_target} PROPERTY TYPE)
  if(NOT "${_type}" STREQUAL "OBJECT_LIBRARY")
    target_link_flags(${_target} PRIVATE
      "${DEAL_II_LINKER_FLAGS} ${DEAL_II_LINKER_FLAGS_${_build}}"
      )
  endif()

  target_include_directories(${_target}
    PRIVATE
      "${CMAKE_BINARY_DIR}/include"
      "${CMAKE_SOURCE_DIR}/include"
    SYSTEM PRIVATE
      ${DEAL_II_BUNDLED_INCLUDE_DIRS}
      ${DEAL_II_INCLUDE_DIRS}
    )

  target_link_libraries(${_target} ${DEAL_II_TARGET_NAME}_${_build_lowercase})
endfunction()
