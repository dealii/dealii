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
#     CLEAR_FEATURE(feature)
#
# and all other suffixes defined in DEAL_II_LIST_SUFFIXES and
# DEAL_II_STRING_SUFFIXES to the corresponding DEAL_II_* variables
#

MACRO(CLEAR_FEATURE _feature)
  FOREACH(_var ${DEAL_II_LIST_SUFFIXES})
    unset(${_feature}_${_var})
  ENDFOREACH()
  FOREACH(_var ${DEAL_II_STRING_SUFFIXES})
    unset(${_feature}_${_var})
  ENDFOREACH()
ENDMACRO()
