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
# A small macro to add libraries to
#   DEAL_II_EXTERNAL_LIBRARIES
#   DEAL_II_EXTERNAL_LIBRARIES_DEBUG
#   DEAL_II_EXTERNAL_LIBRARIES_RELEASE
# depending on the "optmized", "debug" or "general" keyword
#
# Usage:
#   DEAL_II_APPEND_LIBRARIES(<list of libraries>)
#

MACRO(DEAL_II_APPEND_LIBRARIES)

  SET(_toggle "general")
  FOREACH(_tmp ${ARGN})
    IF( "${_tmp}" STREQUAL "debug" OR
        "${_tmp}" STREQUAL "optimized" OR
        "${_tmp}" STREQUAL "general" )
      SET(_toggle "${_tmp}")
    ELSE()
      IF("${_toggle}" STREQUAL "general")
        LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES ${_tmp})
      ELSEIF("${_toggle}" STREQUAL "debug")
        LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES_DEBUG ${_tmp})
      ELSEIF("${_toggle}" STREQUAL "optimized")
        LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES_RELEASE ${_tmp})
      ENDIF()
    ENDIF()
  ENDFOREACH()

ENDMACRO()
