## ---------------------------------------------------------------------
##
## Copyright (C) 2018 - 2020 by the deal.II authors
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
# Try to find the GINKGO library
#
# This module exports
#
#   GINKGO_INCLUDE_DIRS
#

SET(GINKGO_DIR "" CACHE PATH "An optional hint to a GINKGO installation")
SET_IF_EMPTY(GINKGO_DIR "$ENV{GINKGO_DIR}")

DEAL_II_FIND_LIBRARY(GINKGO_LIBRARY
  NAMES ginkgo
  HINTS ${GINKGO_DIR}
  PATH_SUFFIXES
    lib${LIB_SUFFIX} lib64 lib
    # This is a hint, isn't it?
    build/${CMAKE_CXX_PLATFORM_ID}-${CMAKE_SYSTEM_PROCESSOR}/libginkgo
  )
DEAL_II_FIND_LIBRARY(GINKGO_REFERENCE_LIBRARY
  NAMES ginkgo_reference
  HINTS ${GINKGO_DIR}
  PATH_SUFFIXES
    lib${LIB_SUFFIX} lib64 lib
    # This is a hint, isn't it?
    build/${CMAKE_CXX_PLATFORM_ID}-${CMAKE_SYSTEM_PROCESSOR}/libginkgo_reference
  )
DEAL_II_FIND_LIBRARY(GINKGO_OMP_LIBRARY
  NAMES ginkgo_omp
  HINTS ${GINKGO_DIR}
  PATH_SUFFIXES
    lib${LIB_SUFFIX} lib64 lib
    # This is a hint, isn't it?
    build/${CMAKE_CXX_PLATFORM_ID}-${CMAKE_SYSTEM_PROCESSOR}/libginkgo_omp
  )
DEAL_II_FIND_LIBRARY(GINKGO_CUDA_LIBRARY
  NAMES ginkgo_cuda
  HINTS ${GINKGO_DIR}
  PATH_SUFFIXES
    lib${LIB_SUFFIX} lib64 lib
    # This is a hint, isn't it?
    build/${CMAKE_CXX_PLATFORM_ID}-${CMAKE_SYSTEM_PROCESSOR}/libginkgo_cuda
  )

DEAL_II_FIND_PATH(GINKGO_INCLUDE_DIR ginkgo/ginkgo.hpp
  HINTS ${GINKGO_DIR}
  PATH_SUFFIXES include
  )

DEAL_II_PACKAGE_HANDLE(GINKGO
  LIBRARIES
    REQUIRED GINKGO_LIBRARY GINKGO_REFERENCE_LIBRARY GINKGO_OMP_LIBRARY GINKGO_CUDA_LIBRARY
  INCLUDE_DIRS REQUIRED GINKGO_INCLUDE_DIR
  USER_INCLUDE_DIRS REQUIRED GINKGO_INCLUDE_DIR
  CLEAR
    GINKGO_LIBRARY GINKGO_REFERENCE_LIBRARY GINKGO_OMP_LIBRARY GINKGO_CUDA_LIBRARY GINKGO_INCLUDE_DIR
  )
