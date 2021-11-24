## ---------------------------------------------------------------------
##
## Copyright (C) 2021 by the deal.II authors
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
# Try to find the Kokkos library
#
# This module exports
#
#   KOKKOS_INCLUDE_DIRS
#   KOKKOS_INTERFACE_LINK_FLAGS
#

SET(KOKKOS_DIR "" CACHE PATH "An optional hint to a Kokkos installation")
SET_IF_EMPTY(KOKKOS_DIR "$ENV{KOKKOS_DIR}")


IF(DEAL_II_TRILINOS_WITH_KOKKOS)
  # Let ArborX know that we have found Kokkos
  SET(Kokkos_FOUND ON)
  # Let deal.II know that we have found Kokkos
  SET(KOKKOS_FOUND ON)
ELSE()
  FIND_PACKAGE(Kokkos
    HINTS ${KOKKOS_DIR} ${Kokkos_DIR} $ENV{Kokkos_DIR}
    )

  IF(Kokkos_FOUND)
    # We are only interested in Kokkos if it is not part of Trilinos
    IF(TARGET Kokkos::kokkos AND TARGET Kokkos::kokkoscore)
      GET_PROPERTY(KOKKOS_INSTALL_INCLUDE_DIR TARGET Kokkos::kokkos PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
      GET_PROPERTY(KOKKOS_EXTRA_LD_FLAGS_FULL TARGET Kokkos::kokkoscore PROPERTY INTERFACE_LINK_OPTIONS)
      STRING(REGEX REPLACE "\\$<\\$<LINK_LANGUAGE:CXX>:([^>]*)>" "\\1" KOKKOS_EXTRA_LD_FLAGS "${KOKKOS_EXTRA_LD_FLAGS_FULL}")
      GET_PROPERTY(KOKKOS_COMPILE_FLAGS_FULL TARGET Kokkos::kokkoscore PROPERTY INTERFACE_COMPILE_OPTIONS)
      STRING(REGEX REPLACE "\\$<\\$<COMPILE_LANGUAGE:CXX>:([^>]*)>" "\\1" KOKKOS_COMPILE_FLAGS "${KOKKOS_COMPILE_FLAGS_FULL}")

      # We need to disable SIMD vectorization for CUDA device code.
      # Otherwise, nvcc compilers from version 9 on will emit an error message like:
      # "[...] contains a vector, which is not supported in device code". We
      # would like to set the variable in check_01_cpu_feature but at that point
      # we don't know if CUDA support is enabled in Kokkos
      IF(Kokkos_ENABLE_CUDA)
        SET(DEAL_II_VECTORIZATION_WIDTH_IN_BITS 0)
        IF(DEAL_II_WITH_KOKKOS_BACKEND)
          # Require lambda support to use Kokkos as a backend
          KOKKOS_CHECK(OPTIONS CUDA_LAMBDA)
        ENDIF()
      ENDIF()

      DEAL_II_FIND_LIBRARY(KOKKOS_CORE_LIBRARY
        NAMES kokkoscore
        HINTS ${KOKKOS_DIR}/lib ${Kokkos_DIR}/lib $ENV{Kokkos_DIR}/lib
        PATH_SUFFIXES
          lib${LIB_SUFFIX} lib64 lib
        )

      DEAL_II_FIND_LIBRARY(KOKKOS_CONTAINERS_LIBRARY
        NAMES kokkoscontainers
        HINTS ${KOKKOS_DIR}/lib ${Kokkos_DIR}/lib $ENV{Kokkos_DIR}/lib
        PATH_SUFFIXES
          lib${LIB_SUFFIX} lib64 lib
        )
    ELSE()
      SET(Kokkos_FOUND FALSE)
    ENDIF()
  ENDIF()

  DEAL_II_PACKAGE_HANDLE(KOKKOS
    LIBRARIES REQUIRED KOKKOS_CORE_LIBRARY KOKKOS_CONTAINERS_LIBRARY
    INCLUDE_DIRS REQUIRED KOKKOS_INSTALL_INCLUDE_DIR
    USER_INCLUDE_DIRS REQUIRED KOKKOS_INSTALL_INCLUDE_DIR
    CXX_FLAGS OPTIONAL KOKKOS_COMPILE_FLAGS
    LINKER_FLAGS OPTIONAL KOKKOS_EXTRA_LD_FLAGS
    CLEAR KOKKOS_DIR KOKKOS_CORE_LIBRARY KOKKOS_CONTAINERS_LIBRARY
    )
ENDIF()
