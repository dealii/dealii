## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2020 - 2022 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

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
