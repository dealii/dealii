## ---------------------------------------------------------------------
##
## Copyright (C) 2017 by the deal.II authors
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
# Try to find the ASSIMP libraries
#
# This module exports
#
#   ASSIMP_LIBRARIES
#   ASSIMP_INCLUDE_DIRS
#

SET(ASSIMP_DIR "" CACHE PATH "An optional hint to a Assimp installation")
SET_IF_EMPTY(ASSIMP_DIR "$ENV{ASSIMP_DIR}")

DEAL_II_FIND_LIBRARY(ASSIMP_LIB NAMES assimp
  HINTS ${ASSIMP_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

DEAL_II_FIND_PATH(ASSIMP_INC assimp/defs.h
  HINTS ${ASSIMP_DIR}
  PATH_SUFFIXES include
  )

DEAL_II_PACKAGE_HANDLE(ASSIMP
  LIBRARIES REQUIRED ASSIMP_LIB
  INCLUDE_DIRS REQUIRED ASSIMP_INC
  USER_INCLUDE_DIRS REQUIRED ASSIMP_INC
  CLEAR ASSIMP_LIB ASSIMP_INC
  )
