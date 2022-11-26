## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2017 by the deal.II authors
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
# A small wrapper around
# SET_TARget_property(... PROPERTIES COMPILE_DEFINITIONS ...)
# to _add_ compile definitions to every target we have specified.
#

macro(deal_ii_add_definitions _name)

  foreach(_build ${DEAL_II_BUILD_TYPES})
    string(TOLOWER ${_build} _build_lowercase)

    set_property(TARGET ${_name}_${_build_lowercase}
      APPEND PROPERTY COMPILE_DEFINITIONS "${ARGN}"
      )
  endforeach()

endmacro()
