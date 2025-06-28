// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_vector_memory_templates_h
#define dealii_vector_memory_templates_h


#include <deal.II/base/config.h>

#include <deal.II/base/kokkos.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/memory_consumption.h>

#include <deal.II/lac/vector_memory.h>

#include <Kokkos_Core.hpp>

#include <memory>

DEAL_II_NAMESPACE_OPEN


template <typename VectorType>
typename GrowingVectorMemory<VectorType>::Pool &
GrowingVectorMemory<VectorType>::get_pool()
{
  // Kokkos needs to be initialized before constructing the static Pool for
  // vector types that use Kokkos.
  // If Kokkos is initialized by deal.II, this make sure that it is finalized
  // after the Pool has been destroyed.
  // If Kokkos is not initialized by deal.II, we assume that Kokkos is not
  // finalized past program end together with static variables and we need to
  // make sure to empty the Pool when finalizing Kokkos so that the destruction
  // of the Pool doesn't call Kokkos functions.
  static auto pool = []() {
    internal::ensure_kokkos_initialized();
    if (!internal::dealii_initialized_kokkos)
      Kokkos::push_finalize_hook(
        GrowingVectorMemory<VectorType>::release_unused_memory);
    return GrowingVectorMemory<VectorType>::Pool{};
  }();
  return pool;
}



template <typename VectorType>
Threads::Mutex GrowingVectorMemory<VectorType>::mutex;



template <typename VectorType>
inline GrowingVectorMemory<VectorType>::Pool::Pool()
  : data(nullptr)
{}



template <typename VectorType>
inline GrowingVectorMemory<VectorType>::Pool::~Pool()
{
  // Nothing to do if memory was unused.
  if (data == nullptr)
    return;

  // delete the 'data' object. this also releases all vectors
  // that are pointed to by the std::unique_ptrs
  data->clear();
  delete data;
}


template <typename VectorType>
inline void
GrowingVectorMemory<VectorType>::Pool::initialize(const size_type size)
{
  if (data == nullptr)
    {
      data = new std::vector<entry_type>(size);

      for (typename std::vector<entry_type>::iterator i = data->begin();
           i != data->end();
           ++i)
        {
          i->first  = false;
          i->second = std::make_unique<VectorType>();
        }
    }
}


template <typename VectorType>
inline GrowingVectorMemory<VectorType>::GrowingVectorMemory(
  const size_type initial_size,
  const bool      log_statistics)

  : total_alloc(0)
  , current_alloc(0)
  , log_statistics(log_statistics)
{
  // Constructors generally cannot be called in parallel and so it is not
  // necessary to guard access to member variables with a mutex -- unless
  // the member variable being accessed is 'static', which is the case
  // for the things get_pool() returns, and in that case the mutex itself
  // must also be 'static' as it is here.
  std::lock_guard<std::mutex> lock(mutex);
  get_pool().initialize(initial_size);
}


template <typename VectorType>
inline GrowingVectorMemory<VectorType>::~GrowingVectorMemory()
{
  AssertNothrow(current_alloc == 0,
                StandardExceptions::ExcMemoryLeak(current_alloc));
  if (log_statistics)
    {
      deallog << "GrowingVectorMemory:Overall allocated vectors: "
              << total_alloc << std::endl;
      deallog << "GrowingVectorMemory:Maximum allocated vectors: "
              << get_pool().data->size() << std::endl;
    }
}



template <typename VectorType>
inline VectorType *
GrowingVectorMemory<VectorType>::alloc()
{
  std::lock_guard<std::mutex> lock(mutex);

  ++total_alloc;
  ++current_alloc;

  // See if there is a currently unused vector available in our list
  for (entry_type &i : *get_pool().data)
    {
      if (i.first == false)
        {
          i.first = true;
          return i.second.get();
        }
    }

  // No currently unused vector found, so let's just allocate a new one
  // and return it:
  const auto &new_entry =
    get_pool().data->emplace_back(true, std::make_unique<VectorType>());

  return new_entry.second.get();
}



template <typename VectorType>
inline void
GrowingVectorMemory<VectorType>::free(const VectorType *const v)
{
  std::lock_guard<std::mutex> lock(mutex);

  // Find the vector to be de-allocated and mark it as now unused:
  for (entry_type &i : *get_pool().data)
    {
      if (v == i.second.get())
        {
          i.first = false;
          --current_alloc;
          return;
        }
    }

  // If we got here, someone is trying to free a vector that has not
  // been allocated!
  Assert(false, typename VectorMemory<VectorType>::ExcNotAllocatedHere());
}



template <typename VectorType>
inline void
GrowingVectorMemory<VectorType>::release_unused_memory()
{
  std::lock_guard<std::mutex> lock(mutex);

  if (get_pool().data != nullptr)
    get_pool().data->clear();
}



template <typename VectorType>
inline std::size_t
GrowingVectorMemory<VectorType>::memory_consumption() const
{
  std::lock_guard<std::mutex> lock(mutex);

  std::size_t result = sizeof(*this);
  for (const auto &[_, ptr] : *get_pool().data)
    result += sizeof(ptr) + (ptr ? MemoryConsumption::memory_consumption(*ptr) :
                                   MemoryConsumption::memory_consumption(ptr));

  return result;
}


DEAL_II_NAMESPACE_CLOSE

#endif
