## ---------------------------------------------------------------------
##
## Copyright (C) 2022 - 2022 by the deal.II authors
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

  foreach(_property
      TYPE
      LINK_LIBRARIES INTERFACE_LINK_LIBRARIES
      INCLUDE_DIRECTORIES INTERFACE_INCLUDE_DIRECTORIES
      COMPILE_DEFINITIONS INTERFACE_COMPILE_DEFINITIONS
      COMPILE_OPTIONS INTERFACE_COMPILE_OPTIONS
      LINK_OPTIONS INTERFACE_LINK_OPTIONS
      )
    get_target_property(_value ${_target} ${_property})
    if(NOT "${_value}" MATCHES "-NOTFOUND")
      string(REPLACE ";" " " _value "${_value}")
      list(APPEND _messages "  ${_property}: ${_value}")
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

