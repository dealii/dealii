## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2023 - 2024 by the deal.II authors
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
# print_target_properties(<target> [variable])
#
# Prints all relevant information of the specified target.
#

function(print_target_properties _target)

  if(NOT TARGET ${_target})
    return()
  endif()

  set(_messages)
  list(APPEND _messages "Target: ${_target}")

  set(_properties
    INTERFACE_LINK_LIBRARIES INTERFACE_INCLUDE_DIRECTORIES
    INTERFACE_SYSTEM_INCLUDE_DIRECTORIES INTERFACE_COMPILE_DEFINITIONS
    INTERFACE_COMPILE_OPTIONS INTERFACE_LINK_OPTIONS
    )
  if(NOT CMAKE_VERSION VERSION_LESS 3.19)
    set(_properties
      TYPE VERSION SOVERSION LINK_LIBRARIES INCLUDE_DIRECTORIES
      COMPILE_DEFINITIONS COMPILE_FEATURES COMPILE_OPTIONS LINK_OPTIONS
      ${_properties}
      )
  endif()

  foreach(_property ${_properties})
    get_target_property(_value ${_target} ${_property})
    if(NOT "${_value}" MATCHES "-NOTFOUND" AND NOT "${_value}" STREQUAL "")
      string(REPLACE ";" "," _value "${_value}") # workaround: we cannot use ";"
      list(APPEND _messages "    ${_property}: ${_value}")
    endif()
  endforeach()

  if("${ARGN}" STREQUAL "")
    foreach(_message ${_messages})
      message(STATUS "${_message}")
    endforeach()
  else()
    set(${ARGN} "${_messages}" PARENT_SCOPE)
  endif()
endfunction()
