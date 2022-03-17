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

#include <functional>
#include <memory>

#ifdef DEAL_II_USE_KOKKOS_BACKEND
#  include <Kokkos_Core.hpp>
#endif

DEAL_II_NAMESPACE_OPEN

/**
 */
namespace MemorySpace
{
  /**
   * Structure which stores data on the host or the device depending on the
   * template parameter @p MemorySpace. The data is copied into the structure
   * which then owns the data and will release the memory when the destructor is
   * called.
   */
  template <typename T, typename MemorySpace>
  struct MemorySpaceData
  {
    MemorySpaceData();

    /**
     * Return true if there is data stored on the host
     */
    bool
    has_data_on_host();

    /**
     * Return the pointer to the underlying data on the host
     */
    T *
    data();

    /**
     * Return the pointer to the underlying data on the host
     */
    const T *
    data() const;

    /**
     * Return element `i` of the underlying data.
     */
    T &
    operator()(const unsigned int i);

    /**
     * Return element `i` of the underlying data.
     */
    const T &
    operator()(const unsigned int i) const;

    /**
     * Copy the active data (values for Host and values_dev for Device) to @p begin.
     * If the data is on the device it is moved to the host.
     */
    void
    copy_to(T *begin, const std::size_t n_elements);

    /**
     * Copy the data in @p begin to the active data of the structure (values for
     * Host and values_dev for Device). The pointer @p begin must be on the host.
     */
    void
    copy_from(const T *begin, const std::size_t n_elements);

#ifdef DEAL_II_USE_KOKKOS_BACKEND
    /**
     * Kokkos View to the data
     */
    Kokkos::View<T *, MemorySpace> values;

    /**
     * Pointer to data on the host. The pointer points to the same data as
     * values when using shared memory. Otherwise it is not set.
     */
    // The pointer is shared so that MemorySpaceData can be copied and
    // MemorySpaceData::values can be used in Kokkos::parallel_for. This
    // pointer owns the data when using shared memory with MPI. In this case,
    // the Kokkos::View in non-owning. When shared memory with MPI is not used,
    // the pointer is not used.
    std::shared_ptr<T[]> values_sm_ptr;
#else
    /**
     * Pointer to data on the host.
     */
    std::unique_ptr<T[], std::function<void(T *)>> values;

    /**
     * Pointer to data on the device.
     */
    std::unique_ptr<T[]> values_dev;
#endif

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


#  ifdef DEAL_II_USE_KOKKOS_BACKEND
  template <typename T, typename MemorySpace>
  MemorySpaceData<T, MemorySpace>::MemorySpaceData()
    : values((dealii::Impl::ensure_kokkos_initialized(), Kokkos::View<T *, MemorySpace>("memoryspace data", 0)))
  {}



  template <typename T, typename MemorySpace>
  bool
  MemorySpaceData<T, MemorySpace>::has_data_on_host()
  {
    return values.extent(0) > 0;
  }



  template <typename T, typename MemorySpace>
  T *
  MemorySpaceData<T, MemorySpace>::data()
  {
    return values.data();
  }



  template <typename T, typename MemorySpace>
  const T *
  MemorySpaceData<T, MemorySpace>::data() const
  {
    return values.data();
  }


  template <typename T, typename MemorySpace>
  T &
  MemorySpaceData<T, MemorySpace>::operator()(const unsigned int i)
  {
    return values(i);
  }



  template <typename T, typename MemorySpace>
  const T &
  MemorySpaceData<T, MemorySpace>::operator()(const unsigned int i) const
  {
    return values(i);
  }



  template <typename T, typename MemorySpace>
  void
  MemorySpaceData<T, MemorySpace>::copy_to(T *               begin,
                                           const std::size_t n_elements)
  {
    Assert(n_elements <= values.extent(0),
           ExcMessage("n_elements greater than the size of values."));
    Kokkos::parallel_for(
      "MemorySpaceData::copy_to",
      Kokkos::RangePolicy<typename MemorySpace::execution_space>(0, n_elements),
      KOKKOS_LAMBDA(const int i) { *(begin + i) = values(i); });
  }



  template <typename T, typename MemorySpace>
  void
  MemorySpaceData<T, MemorySpace>::copy_from(const T *         begin,
                                             const std::size_t n_elements)
  {
    Assert(n_elements <= values.extent(0),
           ExcMessage("n_elements greater than the size of values."));
    Kokkos::parallel_for(
      "MemorySpaceData::copy_from",
      Kokkos::RangePolicy<typename MemorySpace::execution_space>(0, n_elements),
      KOKKOS_LAMBDA(const int i) { values(i) = *(begin + i); });
  }



  /**
   * Swap function similar to std::swap.
   */
  template <typename T, typename MemorySpace>
  inline void
  swap(MemorySpaceData<T, MemorySpace> &u, MemorySpaceData<T, MemorySpace> &v)
  {
    auto u_copy = Kokkos::create_mirror(Kokkos::WithoutInitializing, u);
    typename MemorySpace::execution_space exec_space;
    // The first two calls to Kokkos::deep_copy are asynchronous. The last call
    // will wait for the three deep_copy to be done before returning.
    Kokkos::deep_copy(exec_space, u_copy, u);
    Kokkos::deep_copy(exec_space, u, v);
    Kokkos::deep_copy(v, u_copy);
  }

#  else

  // Specialization on the Host
  template <typename T>
  struct MemorySpaceData<T, Host>
  {
    MemorySpaceData()
      : values(nullptr, &std::free)
    {}

    bool
    has_data_on_host()
    {
      return values.get() != nullptr;
    }

    T *
    data()
    {
      return values.get();
    }

    const T *
    data() const
    {
      return values.get();
    }

    T &
    operator()(const unsigned int i)
    {
      return values[i];
    }

    const T &
    operator()(const unsigned int i) const
    {
      return values[i];
    }


    void
    copy_to(T *begin, const std::size_t n_elements)
    {
      std::copy(values.get(), values.get() + n_elements, begin);
    }

    void
    copy_from(const T *begin, const std::size_t n_elements)
    {
      std::copy(begin, begin + n_elements, values.get());
    }

    std::unique_ptr<T[], std::function<void(T *)>> values;

    // This is not used but it allows to simplify the code until we start using
    // CUDA-aware MPI.
    std::unique_ptr<T[]> values_dev;

    std::vector<ArrayView<const T>> values_sm;
  };



  template <typename T>
  inline void
  swap(MemorySpaceData<T, Host> &u, MemorySpaceData<T, Host> &v)
  {
    std::swap(u.values, v.values);
  }



#    ifdef DEAL_II_COMPILER_CUDA_AWARE

  // Specialization when using CUDA
  template <typename T>
  struct MemorySpaceData<T, CUDA>
  {
    MemorySpaceData()
      : values(nullptr, &std::free)
      , values_dev(nullptr, Utilities::CUDA::delete_device_data<T>)
    {}

    bool
    has_data_on_host()
    {
      return values.get() != nullptr;
    }

    T *
    data()
    {
      return values.get();
    }

    const T *
    data() const
    {
      return values.get();
    }

    T &
    operator()(const unsigned int i)
    {
      return values[i];
    }

    const T &
    operator()(const unsigned int i) const
    {
      return values[i];
    }

    void
    copy_to(T *begin, const std::size_t n_elements)
    {
      const cudaError_t cuda_error_code = cudaMemcpy(begin,
                                                     values_dev.get(),
                                                     n_elements * sizeof(T),
                                                     cudaMemcpyDeviceToHost);
      AssertCuda(cuda_error_code);
    }

    void
    copy_from(const T *begin, const std::size_t n_elements)
    {
      const cudaError_t cuda_error_code = cudaMemcpy(values_dev.get(),
                                                     begin,
                                                     n_elements * sizeof(T),
                                                     cudaMemcpyHostToDevice);
      AssertCuda(cuda_error_code);
    }

    std::unique_ptr<T[], std::function<void(T *)>> values;
    std::unique_ptr<T[], void (*)(T *)>            values_dev;

    /**
     * This is currently not used.
     */
    std::vector<ArrayView<const T>> values_sm;
  };



  template <typename T>
  inline void
  swap(MemorySpaceData<T, CUDA> &u, MemorySpaceData<T, CUDA> &v)
  {
    std::swap(u.values, v.values);
    std::swap(u.values_dev, v.values_dev);
  }

#    endif

#  endif
#endif

} // namespace MemorySpace

DEAL_II_NAMESPACE_CLOSE

#endif
