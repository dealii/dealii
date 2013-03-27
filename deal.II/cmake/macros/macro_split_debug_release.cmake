#####
##
## Copyright (C) 2013 by the deal.II authors
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
