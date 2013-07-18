## ---------------------------------------------------------------------
## $Id$
##
## Copyright (C) 2013 by the deal.II authors
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
# A small macro to split a list of libraries with debug, optimized and
# general keywords into two lists consisting of all libraries necessary for
# the debug and release target only. If no keyword is given "optimized" is
# assumed.
#
# Usage:
#     SPLIT_DEBUG_RELEASE(list_debug list_release <...list of libraries...>)
#
#

MACRO(SPLIT_DEBUG_RELEASE _list_debug _list_release)

  SET(_toggle "optimized")
  FOREACH(_tmp ${ARGN})
    IF("${_tmp}" STREQUAL "debug" OR
       "${_tmp}" STREQUAL "optimized" OR
       "${_tmp}" STREQUAL "general")
      SET(_toggle "${_tmp}")
    ELSE()
      IF("${_toggle}" STREQUAL "general")
        LIST(APPEND ${_list_debug} "${_tmp}")
        LIST(APPEND ${_list_release} "${_tmp}")
      ELSEIF("${_toggle}" STREQUAL "debug")
        LIST(APPEND ${_list_debug} "${_tmp}")
      ELSEIF("${_toggle}" STREQUAL "optimized")
        LIST(APPEND ${_list_release} "${_tmp}")
      ENDIF()
    ENDIF()
  ENDFOREACH()

  IF("${${_list_debug}}" STREQUAL "")
    SET(${_list_debug} ${${_list_release}})
  ELSEIF("${${_list_release}}" STREQUAL "")
    SET(${_list_release} ${${_list_debug}})
  ENDIF()

ENDMACRO()
