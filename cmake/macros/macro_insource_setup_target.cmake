## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2018 by the deal.II authors
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
    LINK_FLAGS "${DEAL_II_LINKER_FLAGS} ${DEAL_II_LINKER_FLAGS_${_build}}"
    LINKER_LANGUAGE "CXX"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    )

  target_include_directories(${_target}
    PRIVATE
      "${CMAKE_BINARY_DIR}/include"
      "${CMAKE_SOURCE_DIR}/include"
      ${DEAL_II_BUNDLED_INCLUDE_DIRS}
    )
  target_include_directories(${_target} SYSTEM PRIVATE ${DEAL_II_INCLUDE_DIRS})

  if(CMAKE_VERSION VERSION_LESS 3.9 OR CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set_property(TARGET ${_target} APPEND_STRING PROPERTY
      COMPILE_FLAGS " ${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${_build}}"
      )
    set_property(TARGET ${_target} APPEND PROPERTY
      COMPILE_DEFINITIONS "${DEAL_II_DEFINITIONS};${DEAL_II_DEFINITIONS_${_build}}"
      )

  else()

    set(_flags "${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${_build}}")
    separate_arguments(_flags)
    target_compile_options(${_target} PUBLIC
      $<$<COMPILE_LANGUAGE:CXX>:${_flags}>
      )

    target_compile_definitions(${_target}
      PUBLIC ${DEAL_II_DEFINITIONS} ${DEAL_II_DEFINITIONS_${_build}}
      )
  endif()

  get_property(_type TARGET ${_target} PROPERTY TYPE)
  if(NOT "${_type}" STREQUAL "OBJECT_LIBRARY")
    target_link_libraries(${_target} ${DEAL_II_NAMESPACE}_${_build_lowercase})
  endif()

endfunction()
