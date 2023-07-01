## ---------------------------------------------------------------------
##
## Copyright (C) 2014 - 2023 by the deal.II authors
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

set(CUDA_DIR "" CACHE PATH "An optional hint to a CUDA installation")
set_if_empty(CUDA_DIR "$ENV{CUDA_DIR}")

if(NOT "${CUDA_DIR}" STREQUAL "")
  set(CUDA_TOOLKIT_ROOT_DIR "${CUDA_DIR}")
endif()

find_package(CUDA)

if(CUDA_FOUND)
  message(STATUS "Configured to use CUDA installation at ${CUDA_TOOLKIT_ROOT_DIR}")
endif()

enable_if_supported(_cuda_flags "-Xcudafe --diag_suppress=unsigned_compare_with_zero")
enable_if_supported(_cuda_flags "-Xcudafe --diag_suppress=integer_sign_change")
enable_if_supported(_cuda_flags "-Xcudafe --diag_suppress=20208") # long double treated as double
enable_if_supported(_cuda_flags "-Xcudafe --diag_suppress=1301")  # non-template friend
enable_if_supported(_cuda_flags "-Wno-non-template-friend")
enable_if_supported(_cuda_flags "-Xcudafe --diag_suppress=loop_not_reachable")

process_feature(CUDA
  LIBRARIES REQUIRED CUDA_LIBRARIES CUDA_cusparse_LIBRARY CUDA_cusolver_LIBRARY
  INCLUDE_DIRS REQUIRED CUDA_INCLUDE_DIRS
  CXX_FLAGS OPTIONAL _cuda_flags
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
