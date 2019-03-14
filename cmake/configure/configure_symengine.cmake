## ---------------------------------------------------------------------
##
## Copyright (C) 2019 by the deal.II authors
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
# Configuration for the SymEngine library:
#

#
# We require at least version 0.3 of the symengine library:
#
SET(SYMENGINE_MINIMUM_REQUIRED_VERSION "0.3")


MACRO(FEATURE_SYMENGINE_CONFIGURE_EXTERNAL)
  SET(DEAL_II_SYMENGINE_WITH_LLVM ${SYMENGINE_WITH_LLVM})

  IF(DEAL_II_SYMENGINE_WITH_LLVM)
    MESSAGE(STATUS "Configured with SymEngine LLVM capabilities.")
  ENDIF()

  #
  # Overwrite the compiler flags imported from SymEngine
  #
  SET(SYMENGINE_CXX_FLAGS_DEBUG)
  SET(SYMENGINE_CXX_FLAGS_RELEASE)
ENDMACRO()


CONFIGURE_FEATURE(SYMENGINE)
