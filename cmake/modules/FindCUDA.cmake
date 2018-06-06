## ---------------------------------------------------------------------
##
## Copyright (C) 2014 - 2018 by the deal.II authors
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
# Try to find cuda
#
# This module exports:
#
#     CUDA_FOUND
#     CUDA_INCLUDE_DIRS
#     CUDA_LIBRARIES
#     CUDA_NVCC_EXECUTABLE
#     CUDA_VERSION
#     CUDA_VERSION_MAJOR
#     CUDA_VERSION_MINOR
#

SET(CUDA_DIR "" CACHE PATH "An optional hint to a CUDA installation")
SET_IF_EMPTY(CUDA_DIR "$ENV{CUDA_DIR}")

IF(NOT "${CUDA_DIR}" STREQUAL "")
  SET(CUDA_TOOLKIT_ROOT_DIR "${CUDA_DIR}")
ENDIF()

# temporarily disable ${CMAKE_SOURCE_DIR}/cmake/modules for module lookup
LIST(REMOVE_ITEM CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)
FIND_PACKAGE(CUDA)
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)

IF(CUDA_FOUND)
  MESSAGE(STATUS "Configured to use CUDA installation at ${CUDA_TOOLKIT_ROOT_DIR}")
ENDIF()

# cuSOLVER requires OpenMP
FIND_PACKAGE(OpenMP)
SET(_cuda_libraries ${CUDA_LIBRARIES} ${CUDA_cusparse_LIBRARY}
  ${CUDA_cusolver_LIBRARY} ${OpenMP_CXX_FLAGS})
SET(_cuda_include_dirs ${CUDA_INCLUDE_DIRS})
DEAL_II_PACKAGE_HANDLE(CUDA
  LIBRARIES REQUIRED _cuda_libraries
  INCLUDE_DIRS REQUIRED  _cuda_include_dirs
  USER_INCLUDE_DIRS REQUIRED _cuda_include_dirs
  CLEAR
    CUDA_cublas_device_LIBRARY
    CUDA_cublas_LIBRARY
    CUDA_cudadevrt_LIBRARY
    CUDA_cudart_static_LIBRARY
    CUDA_cufft_LIBRARY
    CUDA_cupti_LIBRARY
    CUDA_curand_LIBRARY
    CUDA_HOST_COMPILER
    CUDA_nppc_LIBRARY
    CUDA_nppi_LIBRARY
    CUDA_npps_LIBRARY
    CUDA_rt_LIBRARY
    CUDA_SDK_ROOT_DIR
    CUDA_TOOLKIT_ROOT_DIR
    CUDA_USE_STATIC_CUDA_RUNTIME
  )
