## ---------------------------------------------------------------------
##
## Copyright (C) 2020 by the deal.II authors
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
# This macro is used for the feature configuration in deal.II. It clears
# all individual FEATURE_* configuration variables for a given feature.
#
# Usage:
#     clear_feature(feature)
#
# and all other suffixes defined in DEAL_II_LIST_SUFFIXES and
# DEAL_II_STRING_SUFFIXES to the corresponding DEAL_II_* variables
#

macro(clear_feature _feature)
  foreach(_var ${DEAL_II_LIST_SUFFIXES} ${DEAL_II_STRING_SUFFIXES})
    unset(${_feature}_${_var})
    #
    # If an unset variable is still defined then it is a cached variable.
    # We want to retain this cached variable because it has been set by the
    # user (and might have been used for the feature configuration). But it
    # is important that we clear the variable so that the process_feature()
    # macro can function as intended. As a workaround let's simply set it
    # to the empty string.
    #
    if(DEFINED ${_feature}_${_var})
      set(${_feature}_${_var} "")
    endif()
  endforeach()
  unset(${_FEATURE}_FOUND)
  unset(${_FEATURE}_SPLIT_CONFIGURATION)
endmacro()
