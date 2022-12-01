// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <Kokkos_Core.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#include <functional>
#include <memory>

DEAL_II_NAMESPACE_OPEN

/**
 */
namespace MemorySpace
{
  /**
   * Data structure
   */
  template <typename T, typename MemorySpace>
  struct MemorySpaceData
  {
    MemorySpaceData()
    {
      static_assert(std::is_same<MemorySpace, Host>::value ||
                      std::is_same<MemorySpace, Device>::value,
                    "MemorySpace should be Host or Device");
    }

    /**
     * Copy the active data (values for Host and values_dev for Device) to @p begin.
     * If the data is on the device it is moved to the host.
     */
    void
    copy_to(T *begin, std::size_t n_elements)
    {
      (void)begin;
      (void)n_elements;
    }

    /**
     * Copy the data in @p begin to the active data of the structure (values for
     * Host and values_dev for Device). The pointer @p begin must be on the host.
     */
    void
    copy_from(T *begin, std::size_t n_elements)
    {
      (void)begin;
      (void)n_elements;
    }

    /**
     * Pointer to data on the host.
     */
    Kokkos::View<T *, Kokkos::HostSpace> values;

    /**
     * Pointer to data on the device.
     */
    Kokkos::View<T *, typename MemorySpace::kokkos_space> values_dev;

    /**
     * Pointer to data on the host. The pointer points to the same data as
     * values when using shared memory. Otherwise it is not set.
     */
    // The pointer is shared so that MemorySpaceData can be copied and
    // MemorySpaceData::values can be used in Kokkos::parallel_for. This
    // pointer owns the data when using shared memory with MPI. In this case,
    // the (host) Kokkos::View is non-owning. When shared memory with MPI is not
    // used, the pointer is not used.
    std::shared_ptr<T> values_sm_ptr;

    /**
     * Pointers to the data of the processes sharing the same memory. Not used for MemorySpace::Device.
     */
    std::vector<ArrayView<const T>> values_sm;
  };



  /**
   * Swap function similar to std::swap.
   */
  template <typename T, typename MemorySpace>
  inline void
  swap(MemorySpaceData<T, MemorySpace> &,
       MemorySpaceData<T, MemorySpace> &)
  {
    static_assert(std::is_same<MemorySpace, Host>::value ||
                    std::is_same<MemorySpace, Device>::value,
                  "MemorySpace should be Host or Device");
  }

#ifndef DOXYGEN

  template <typename T>
  struct MemorySpaceData<T, Host>
  {
    using MemorySpace = Host;

    MemorySpaceData()
    : values((dealii::Impl::ensure_kokkos_initialized(),
              Kokkos::View<T *, Kokkos::HostSpace>("host data", 0)))
    {}

    void
    copy_to(T *begin, std::size_t n_elements)
    {
      Assert(n_elements <= values.extent(0),
             ExcMessage("n_elements greater than the size of values."));
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

    void
    copy_from(T *begin, std::size_t n_elements)
    {
      Assert(n_elements <= values.extent(0),
             ExcMessage("n_elements greater than the size of values."));
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

    Kokkos::View<T *, Kokkos::HostSpace> values;

    // unused
    Kokkos::View<T *, typename MemorySpace::kokkos_space> values_dev;

    std::shared_ptr<T> values_sm_ptr;

    std::vector<ArrayView<const T>> values_sm;
  };



  template <typename T>
  inline void
  swap(MemorySpaceData<T, Host> &u, MemorySpaceData<T, Host> &v)
  {
    std::swap(u.values, v.values);
    std::swap(u.values_sm_ptr, v.values_sm_ptr);
  }



  template <typename T>
  struct MemorySpaceData<T, Device>
  {
    using MemorySpace = Device;

    MemorySpaceData()
    : values((dealii::Impl::ensure_kokkos_initialized(),
              Kokkos::View<T *, Kokkos::HostSpace>("host data", 0)))
    , values_dev(Kokkos::View<T *, typename MemorySpace::kokkos_space>(
        "memoryspace data",
        0))
    {}

    void
    copy_to(T *begin, std::size_t n_elements)
    {
      Assert(n_elements <= values_dev.extent(0),
             ExcMessage("n_elements greater than the size of values."));
      using ExecutionSpace = typename MemorySpace::kokkos_space::execution_space;
      Kokkos::
        View<T *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
          begin_view(begin, n_elements);
      Kokkos::deep_copy(
        ExecutionSpace{},
        begin_view,
        Kokkos::subview(values_dev, Kokkos::make_pair(std::size_t(0), n_elements)));
      ExecutionSpace{}.fence();
    }

    void
    copy_from(T *begin, std::size_t n_elements)
    {
      Assert(n_elements <= values_dev.extent(0),
             ExcMessage("n_elements greater than the size of values."));
      using ExecutionSpace = typename MemorySpace::kokkos_space::execution_space;
      Kokkos::View<const T *,
                   Kokkos::HostSpace,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        begin_view(begin, n_elements);
      Kokkos::deep_copy(
        ExecutionSpace{},
        Kokkos::subview(values_dev, Kokkos::make_pair(std::size_t(0), n_elements)),
        begin_view);
      ExecutionSpace{}.fence();
    }

    Kokkos::View<T *, Kokkos::HostSpace> values;
    
    Kokkos::View<T *, typename MemorySpace::kokkos_space> values_dev;
    
    // unused
    std::shared_ptr<T> values_sm_ptr;
  
    // unused  
    std::vector<ArrayView<const T>> values_sm;
  };



  template <typename T>
  inline void
  swap(MemorySpaceData<T, Device> &u, MemorySpaceData<T, Device> &v)
  {
    std::swap(u.values, v.values);
    std::swap(u.values_dev, v.values_dev);
  }

#endif

} // namespace MemorySpace

DEAL_II_NAMESPACE_CLOSE

#endif
