## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2013 by the deal.II authors
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
# A small wrapper around ADD_DEPENDENCIES to add the specified dependencies
# to every ${target}_${build} target, where build runs through all build
# types specified in DEAL_II_BUILD_TYPES
#

MACRO(DEAL_II_ADD_DEPENDENCIES _name _target)

  FOREACH(_build ${DEAL_II_BUILD_TYPES})
    STRING(TOLOWER ${_build} _build_lowercase)
    ADD_DEPENDENCIES(${_name}.${_build_lowercase}
      ${_target}.${_build_lowercase}
      )
  ENDFOREACH()

ENDMACRO()
