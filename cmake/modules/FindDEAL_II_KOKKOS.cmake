## ---------------------------------------------------------------------
##
## Copyright (C) 2021 - 2022 by the deal.II authors
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

set(KOKKOS_DIR "" CACHE PATH "An optional hint to a Kokkos installation")
set_if_empty(KOKKOS_DIR "$ENV{KOKKOS_DIR}")


if(DEAL_II_TRILINOS_WITH_KOKKOS OR DEAL_II_PETSC_WITH_KOKKOS)
  # Let ArborX know that we have found Kokkos
  set(Kokkos_FOUND ON)
  # Let deal.II know that we have found Kokkos
  set(KOKKOS_FOUND ON)
else()

  # silence a warning when including FindKOKKOS.cmake
  set(CMAKE_CXX_EXTENSIONS OFF)
  find_package(Kokkos 3.7.0 QUIET
    HINTS ${KOKKOS_DIR} ${Kokkos_DIR} $ENV{Kokkos_DIR}
    )

  if(Kokkos_FOUND)
    # We are only interested in Kokkos if it is not part of Trilinos
    if(TARGET Kokkos::kokkos AND TARGET Kokkos::kokkoscore)
      get_property(KOKKOS_INSTALL_INCLUDE_DIR TARGET Kokkos::kokkos PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
      get_property(KOKKOS_EXTRA_LD_FLAGS_FULL TARGET Kokkos::kokkoscore PROPERTY INTERFACE_LINK_OPTIONS)
      string(REGEX REPLACE "\\$<\\$<LINK_LANGUAGE:CXX>:([^>]*)>" "\\1" KOKKOS_EXTRA_LD_FLAGS "${KOKKOS_EXTRA_LD_FLAGS_FULL}")
      string(REPLACE ";" " " KOKKOS_EXTRA_LD_FLAGS "${KOKKOS_EXTRA_LD_FLAGS}")
      get_property(KOKKOS_COMPILE_FLAGS_FULL TARGET Kokkos::kokkoscore PROPERTY INTERFACE_COMPILE_OPTIONS)
      string(REGEX REPLACE "\\$<\\$<COMPILE_LANGUAGE:CXX>:([^>]*)>" "\\1" KOKKOS_COMPILE_FLAGS "${KOKKOS_COMPILE_FLAGS_FULL}")
      string(REPLACE ";" " " KOKKOS_COMPILE_FLAGS "${KOKKOS_COMPILE_FLAGS}")

      # Kokkos links transitively with OpenMP so that we need to find OpenMP again.
      # Since we are not using target_link_libraries, we have to extract the compile flag manually.
      if(Kokkos_VERSION VERSION_GREATER_EQUAL 4.0.00 AND Kokkos_ENABLE_OPENMP)
        string(APPEND KOKKOS_COMPILE_FLAGS " ${OpenMP_CXX_FLAGS}")
      endif()

      # We need to disable SIMD vectorization for CUDA device code.
      # Otherwise, nvcc compilers from version 9 on will emit an error message like:
      # "[...] contains a vector, which is not supported in device code". We
      # would like to set the variable in check_01_cpu_feature but at that point
      # we don't know if CUDA support is enabled in Kokkos
      if(Kokkos_ENABLE_CUDA)
        set(DEAL_II_VECTORIZATION_WIDTH_IN_BITS 0)
        # Require lambda support to use Kokkos as a backend
        KOKKOS_CHECK(OPTIONS CUDA_LAMBDA)
      endif()

      deal_ii_find_library(KOKKOS_CORE_LIBRARY
        NAMES kokkoscore
        HINTS ${KOKKOS_DIR}/lib ${Kokkos_DIR}/lib $ENV{Kokkos_DIR}/lib
        PATH_SUFFIXES
          lib${LIB_SUFFIX} lib64 lib
        )

      deal_ii_find_library(KOKKOS_CONTAINERS_LIBRARY
        NAMES kokkoscontainers
        HINTS ${KOKKOS_DIR}/lib ${Kokkos_DIR}/lib $ENV{Kokkos_DIR}/lib
        PATH_SUFFIXES
          lib${LIB_SUFFIX} lib64 lib
        )
    else()
      set(Kokkos_FOUND FALSE)
    endif()
  endif()

  process_feature(KOKKOS
    LIBRARIES REQUIRED KOKKOS_CORE_LIBRARY KOKKOS_CONTAINERS_LIBRARY
    INCLUDE_DIRS REQUIRED KOKKOS_INSTALL_INCLUDE_DIR
    CXX_FLAGS OPTIONAL KOKKOS_COMPILE_FLAGS
    LINKER_FLAGS OPTIONAL KOKKOS_EXTRA_LD_FLAGS
    CLEAR KOKKOS_DIR KOKKOS_CORE_LIBRARY KOKKOS_CONTAINERS_LIBRARY
    )
endif()
