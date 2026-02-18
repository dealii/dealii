## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2014 - 2022 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

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
