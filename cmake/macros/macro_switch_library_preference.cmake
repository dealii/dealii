## ---------------------------------------------------------------------
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
# This macro toggles the preference for static/shared libraries if
# DEAL_II_PREFER_STATIC_LIBS=TRUE but the final executable will still be
# dynamically linked, i.e. DEAL_II_STATIC_EXECUTABLE=OFF
#
# Usage:
#     SWITCH_LIBRARY_PREFERENCE()
#

MACRO(SWITCH_LIBRARY_PREFERENCE)
  IF(DEAL_II_PREFER_STATIC_LIBS AND NOT DEAL_II_STATIC_EXECUTABLE)
    #
    # Invert the search order for libraries when DEAL_II_PREFER_STATIC_LIBS
    # is set. This will prefer static archives instead of shared libraries:
    LIST(REVERSE CMAKE_FIND_LIBRARY_SUFFIXES)
  ENDIF()
ENDMACRO()
