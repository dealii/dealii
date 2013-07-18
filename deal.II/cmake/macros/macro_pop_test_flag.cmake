## ---------------------------------------------------------------------
## $Id$
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
# A small macro used in the platform checks to remove the right most flag in
# CMAKE_REQUIRED_FLAGS
#
# We assume that the flags in CMAKE_REQUIRED_FLAGS are space separated
#
# Usage:
#     POP_TEST_FLAG()
#

MACRO(POP_TEST_FLAG)
  SET(CMAKE_REQUIRED_FLAGS " ${CMAKE_REQUIRED_FLAGS}")
  STRING(REGEX REPLACE " [^ ]+$" ""
    CMAKE_REQUIRED_FLAGS
    "${CMAKE_REQUIRED_FLAGS}"
    )
  STRING(STRIP "${CMAKE_REQUIRED_FLAGS}" CMAKE_REQUIRED_FLAGS)
ENDMACRO()

