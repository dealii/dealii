// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


#ifndef dealii_memory_space_data_h
#define dealii_memory_space_data_h

#include <deal.II/base/config.h>

#include <deal.II/base/cuda.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/kokkos.h>

#include <Kokkos_Core.hpp>

#include <functional>
#include <memory>

DEAL_II_NAMESPACE_OPEN

/**
 */
namespace MemorySpace
{
  /**
   * Structure which stores data on the host or the @ref GlossDevice "device" depending on the
   * template parameter @p MemorySpace. Valid choices are MemorySpace::Host,
   * MemorySpace::Default, and MemorySpace::CUDA (if CUDA was enabled in
   * deal.II). The data is copied into the structure which then owns the data
   * and will release the memory when the destructor is called.
   */
  template <typename T, typename MemorySpace>
  struct MemorySpaceData
  {
    MemorySpaceData();

    /**
     * Copy the class member values to @p begin.
     * If the data is on the @ref GlossDevice "device" it is moved to the host.
     */
    void
    copy_to(T *begin, const std::size_t n_elements);

    /**
     * Copy the data in @p begin to the class member values.
     * The pointer @p begin must be on the host.
     */
    void
    copy_from(const T *begin, const std::size_t n_elements);

    /**
     * Kokkos View owning a host buffer used for MPI communication.
     */
    // FIXME Should we move this somewhere else?
#if KOKKOS_VERSION < 40000
    Kokkos::View<T *, Kokkos::HostSpace> values_host_buffer;
#else
    Kokkos::View<T *, Kokkos::SharedHostPinnedSpace> values_host_buffer;
#endif

    /**
     * Kokkos View owning the data on the @ref GlossDevice "device" (unless @p values_sm_ptr is used).
     */
    Kokkos::View<T *, typename MemorySpace::kokkos_space> values;

    /**
     * Pointer to data on the host. The pointer points to the same data as
     * @p values when using shared memory and the memory space is
     * MemorySpace::Host. Otherwise it is not set.
     */
    // This a shared pointer so that MemorySpaceData can be copied and
    // MemorySpaceData::values can be used in Kokkos::parallel_for. This
    // pointer owns the data when using shared memory with MPI. In this case,
    // the Kokkos::View @p values is non-owning. When shared memory with MPI is
    // not used, the @p values_sm_ptr is unused.
    std::shared_ptr<T> values_sm_ptr;

    /**
     * Pointers to the data of the processes sharing the same memory.
     */
    std::vector<ArrayView<const T>> values_sm;
  };



  /**
   * Swap function similar to std::swap.
   */
  template <typename T, typename MemorySpace>
  inline void
  swap(MemorySpaceData<T, MemorySpace> &u, MemorySpaceData<T, MemorySpace> &v);


#ifndef DOXYGEN

  template <typename T, typename MemorySpace>
  MemorySpaceData<T, MemorySpace>::MemorySpaceData()
    : values_host_buffer(
        (dealii::internal::ensure_kokkos_initialized(),
#  if KOKKOS_VERSION < 40000
         Kokkos::View<T *, Kokkos::HostSpace>("host buffer", 0)))
#  else
         Kokkos::View<T *, Kokkos::SharedHostPinnedSpace>("host pinned buffer",
                                                          0)))
#  endif
    , values(Kokkos::View<T *, typename MemorySpace::kokkos_space>(
        "memoryspace data",
        0))
  {}



  template <typename T, typename MemorySpace>
  void
  MemorySpaceData<T, MemorySpace>::copy_to(T                *begin,
                                           const std::size_t n_elements)
  {
    Assert(n_elements <= values.extent(0),
           ExcMessage(
             "n_elements is greater than the size of MemorySpaceData."));
    using ExecutionSpace = typename MemorySpace::kokkos_space::execution_space;
    Kokkos::
      View<T *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        begin_view(begin, n_elements);
    Kokkos::deep_copy(
      ExecutionSpace{},
      begin_view,
      Kokkos::subview(values, Kokkos::make_pair(std::size_t(0), n_elements)));
    ExecutionSpace{}.fence();
  }



  template <typename T, typename MemorySpace>
  void
  MemorySpaceData<T, MemorySpace>::copy_from(const T          *begin,
                                             const std::size_t n_elements)
  {
    Assert(n_elements <= values.extent(0),
           ExcMessage(
             "n_elements is greater than the size of MemorySpaceData."));
    using ExecutionSpace = typename MemorySpace::kokkos_space::execution_space;
    Kokkos::View<const T *,
                 Kokkos::HostSpace,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      begin_view(begin, n_elements);
    Kokkos::deep_copy(
      ExecutionSpace{},
      Kokkos::subview(values, Kokkos::make_pair(std::size_t(0), n_elements)),
      begin_view);
    ExecutionSpace{}.fence();
  }



  /**
   * Swap function similar to std::swap.
   */
  template <typename T, typename MemorySpace>
  inline void
  swap(MemorySpaceData<T, MemorySpace> &u, MemorySpaceData<T, MemorySpace> &v)
  {
    std::swap(u.values_host_buffer, v.values_host_buffer);
    std::swap(u.values, v.values);
    std::swap(u.values_sm_ptr, v.values_sm_ptr);
  }

#endif

} // namespace MemorySpace

DEAL_II_NAMESPACE_CLOSE

#endif
