## ---------------------------------------------------------------------
##
## Copyright (C) 2014 by the deal.II authors
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
# Remove all cached and non cached variables associated with a feature.
#
# Usage:
#     purge_feature(feature)
#

macro(purge_feature _feature)
  #
  # uncached:
  #
  clear_feature(${_feature})
  unset(${_feature}_FOUND)
  unset(${_feature}_VERSION)

  #
  # cached:
  #
  foreach(_var ${${_feature}_CLEAR_VARIABLES})
    set(${_var})
    unset(${_var} CACHE)
  endforeach()

  unset(${_feature}_CLEAR_VARIABLES CACHE)

  mark_as_advanced(CLEAR ${_feature}_DIR ${_feature}_ARCH)
endmacro()
