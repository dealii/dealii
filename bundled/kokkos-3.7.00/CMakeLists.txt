## ---------------------------------------------------------------------
##
## Copyright (C) 2022 by the deal.II authors
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

set(src_kokkos
  #tpls/desul/src/Lock_Array_CUDA.cpp
  #tpls/desul/src/Lock_Array_HIP.cpp
  core/src/impl/Kokkos_Profiling.cpp
  core/src/impl/Kokkos_NumericTraits.cpp
  core/src/impl/Kokkos_HostSpace_deepcopy.cpp
  core/src/impl/Kokkos_Error.cpp
  core/src/impl/Kokkos_Spinwait.cpp
  core/src/impl/Kokkos_MemoryPool.cpp
  core/src/impl/Kokkos_HostThreadTeam.cpp
  core/src/impl/Kokkos_hwloc.cpp
  core/src/impl/Kokkos_MemorySpace.cpp
  core/src/impl/Kokkos_HBWSpace.cpp
  core/src/impl/Kokkos_HostBarrier.cpp
  core/src/impl/Kokkos_Core.cpp
  core/src/impl/Kokkos_HostSpace.cpp
  core/src/impl/Kokkos_CPUDiscovery.cpp
  core/src/impl/Kokkos_Stacktrace.cpp
  core/src/impl/Kokkos_SharedAlloc.cpp
  core/src/impl/Kokkos_Command_Line_Parsing.cpp
  core/src/impl/Kokkos_ExecPolicy.cpp
  #core/src/OpenMPTarget/Kokkos_OpenMPTargetSpace.cpp
  #core/src/OpenMPTarget/Kokkos_OpenMPTarget_Task.cpp
  #core/src/OpenMPTarget/Kokkos_OpenMPTarget_Instance.cpp
  #core/src/OpenMPTarget/Kokkos_OpenMPTarget_Exec.cpp
  #core/src/OpenACC/Kokkos_OpenACCSpace.cpp
  #core/src/OpenACC/Kokkos_OpenACC.cpp
  #core/src/OpenACC/Kokkos_OpenACC_Instance.cpp
  #core/src/Cuda/Kokkos_Cuda_Instance.cpp
  #core/src/Cuda/Kokkos_Cuda_Task.cpp
  #core/src/Cuda/Kokkos_CudaSpace.cpp
  #core/src/Cuda/Kokkos_Cuda_Locks.cpp
  #core/src/HPX/Kokkos_HPX_Task.cpp
  #core/src/HPX/Kokkos_HPX.cpp
  core/src/dummy.cpp
  #core/src/Threads/Kokkos_ThreadsExec.cpp
  #core/src/HIP/Kokkos_HIP_Instance.cpp
  #core/src/HIP/Kokkos_HIP_Space.cpp
  #core/src/HIP/Kokkos_HIP_Locks.cpp
  #core/src/SYCL/Kokkos_SYCL_Space.cpp
  #core/src/SYCL/Kokkos_SYCL.cpp
  #core/src/SYCL/Kokkos_SYCL_Instance.cpp
  core/src/Serial/Kokkos_Serial.cpp
  core/src/Serial/Kokkos_Serial_Task.cpp
  #core/src/OpenMP/Kokkos_OpenMP_Instance.cpp
  #core/src/OpenMP/Kokkos_OpenMP_Task.cpp
  algorithms/src/KokkosAlgorithms_dummy.cpp
  simd/src/Kokkos_SIMD_dummy.cpp
  containers/src/impl/Kokkos_UnorderedMap_impl.cpp
)

enable_if_supported(DEAL_II_CXX_FLAGS -Wno-float-conversion)
enable_if_supported(DEAL_II_CXX_FLAGS -Wno-missing-field-initializers)
enable_if_supported(DEAL_II_CXX_FLAGS -Wno-suggest-override)
enable_if_supported(DEAL_II_CXX_FLAGS -Wno-unused-but-set-parameter)

deal_ii_add_library(obj_kokkos OBJECT ${src_kokkos})
