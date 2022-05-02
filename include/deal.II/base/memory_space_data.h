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
  template <typename Number, typename MemorySpace>
  struct MemorySpaceData
  {
    MemorySpaceData()
    {
      static_assert(std::is_same<MemorySpace, Host>::value ||
                      std::is_same<MemorySpace, CUDA>::value,
                    "MemorySpace should be Host or CUDA");
    }

    /**
     * Copy the active data (values for Host and values_dev for CUDA) to @p begin.
     * If the data is on the device it is moved to the host.
     */
    void
    copy_to(Number *begin, std::size_t n_elements)
    {
      (void)begin;
      (void)n_elements;
    }

    /**
     * Copy the data in @p begin to the active data of the structure (values for
     * Host and values_dev for CUDA). The pointer @p begin must be on the host.
     */
    void
    copy_from(Number *begin, std::size_t n_elements)
    {
      (void)begin;
      (void)n_elements;
    }

    /**
     * Pointer to data on the host.
     */
    std::unique_ptr<Number[], std::function<void(Number *)>> values;

    /**
     * Pointer to data on the device.
     */
    std::unique_ptr<Number[]> values_dev;

    /**
     * Pointers to the data of the processes sharing the same memory.
     */
    std::vector<ArrayView<const Number>> values_sm;
  };



  /**
   * Swap function similar to std::swap.
   */
  template <typename Number, typename MemorySpace>
  inline void
  swap(MemorySpaceData<Number, MemorySpace> &,
       MemorySpaceData<Number, MemorySpace> &)
  {
    static_assert(std::is_same<MemorySpace, Host>::value ||
                    std::is_same<MemorySpace, CUDA>::value,
                  "MemorySpace should be Host or CUDA");
  }

#ifndef DOXYGEN

  template <typename Number>
  struct MemorySpaceData<Number, Host>
  {
    MemorySpaceData()
      : values(nullptr, &std::free)
    {}

    void
    copy_to(Number *begin, std::size_t n_elements)
    {
      std::copy(values.get(), values.get() + n_elements, begin);
    }

    void
    copy_from(Number *begin, std::size_t n_elements)
    {
      std::copy(begin, begin + n_elements, values.get());
    }

    std::unique_ptr<Number[], std::function<void(Number *)>> values;

    // This is not used but it allows to simplify the code until we start using
    // CUDA-aware MPI.
    std::unique_ptr<Number[]> values_dev;

    std::vector<ArrayView<const Number>> values_sm;
  };



  template <typename Number>
  inline void
  swap(MemorySpaceData<Number, Host> &u, MemorySpaceData<Number, Host> &v)
  {
    std::swap(u.values, v.values);
  }



#  ifdef DEAL_II_COMPILER_CUDA_AWARE

  template <typename Number>
  struct MemorySpaceData<Number, CUDA>
  {
    MemorySpaceData()
      : values(nullptr, &std::free)
      , values_dev(nullptr, Utilities::CUDA::delete_device_data<Number>)
    {}

    void
    copy_to(Number *begin, std::size_t n_elements)
    {
      const cudaError_t cuda_error_code =
        cudaMemcpy(begin,
                   values_dev.get(),
                   n_elements * sizeof(Number),
                   cudaMemcpyDeviceToHost);
      AssertCuda(cuda_error_code);
    }

    void
    copy_from(Number *begin, std::size_t n_elements)
    {
      const cudaError_t cuda_error_code =
        cudaMemcpy(values_dev.get(),
                   begin,
                   n_elements * sizeof(Number),
                   cudaMemcpyHostToDevice);
      AssertCuda(cuda_error_code);
    }

    std::unique_ptr<Number[], std::function<void(Number *)>> values;
    std::unique_ptr<Number[], void (*)(Number *)>            values_dev;

    /**
     * This is currently not used.
     */
    std::vector<ArrayView<const Number>> values_sm;
  };



  template <typename Number>
  inline void
  swap(MemorySpaceData<Number, CUDA> &u, MemorySpaceData<Number, CUDA> &v)
  {
    std::swap(u.values, v.values);
    std::swap(u.values_dev, v.values_dev);
  }

#  endif

#endif

} // namespace MemorySpace

DEAL_II_NAMESPACE_CLOSE

#endif
