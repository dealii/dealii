## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2022 - 2023 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

# Collect files that need to be compiled. This list is created via the following
# command:
#   find . -name *.cpp | egrep -v '(OpenMP|OpenACC|CUDA|Cuda|SYCL|HIP|HPX)'
set(src_kokkos
  simd/src/Kokkos_SIMD_dummy.cpp
  containers/src/impl/Kokkos_UnorderedMap_impl.cpp
  algorithms/src/KokkosAlgorithms_dummy.cpp
  core/src/Serial/Kokkos_Serial.cpp
  core/src/Serial/Kokkos_Serial_Task.cpp
  core/src/impl/Kokkos_SharedAlloc.cpp
  core/src/impl/Kokkos_HostThreadTeam.cpp
  core/src/impl/Kokkos_ExecPolicy.cpp
  core/src/impl/Kokkos_HostSpace.cpp
  core/src/impl/Kokkos_Core.cpp
  core/src/impl/Kokkos_MemoryPool.cpp
  core/src/impl/Kokkos_Profiling.cpp
  core/src/impl/Kokkos_hwloc.cpp
  core/src/impl/Kokkos_Command_Line_Parsing.cpp
  core/src/impl/Kokkos_HostBarrier.cpp
  core/src/impl/Kokkos_CPUDiscovery.cpp
  core/src/impl/Kokkos_Abort.cpp
  core/src/impl/Kokkos_HostSpace_deepcopy.cpp
  core/src/impl/Kokkos_Error.cpp
  core/src/impl/Kokkos_Stacktrace.cpp
)

enable_if_supported(DEAL_II_WARNING_FLAGS -Wno-float-conversion)
enable_if_supported(DEAL_II_WARNING_FLAGS -Wno-missing-field-initializers)
enable_if_supported(DEAL_II_WARNING_FLAGS -Wno-suggest-override)
enable_if_supported(DEAL_II_WARNING_FLAGS -Wno-unused-but-set-parameter)


include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/algorithms/src
  ${CMAKE_CURRENT_SOURCE_DIR}/containers/src
  ${CMAKE_CURRENT_SOURCE_DIR}/core/src
  ${CMAKE_CURRENT_SOURCE_DIR}/simd/src
  ${CMAKE_CURRENT_SOURCE_DIR}/tpls/mdspan/include
  ${CMAKE_CURRENT_SOURCE_DIR}/tpls/desul/include
  )

define_object_library(bundled_kokkos OBJECT ${src_kokkos})
