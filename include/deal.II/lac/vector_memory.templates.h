// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2018 by the deal.II authors
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

#ifndef dealii_vector_memory_templates_h
#define dealii_vector_memory_templates_h


#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/std_cxx14/memory.h>

#include <deal.II/lac/vector_memory.h>

DEAL_II_NAMESPACE_OPEN


template <typename VectorType>
typename GrowingVectorMemory<VectorType>::Pool &
GrowingVectorMemory<VectorType>::get_pool()
{
  static GrowingVectorMemory<VectorType>::Pool pool;
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
          i->second = std_cxx14::make_unique<VectorType>();
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
  // see if there is a free vector
  // available in our list
  for (typename std::vector<entry_type>::iterator i = get_pool().data->begin();
       i != get_pool().data->end();
       ++i)
    {
      if (i->first == false)
        {
          i->first = true;
          return i->second.get();
        }
    }

  // no free vector found, so let's just allocate a new one
  get_pool().data->emplace_back(true, std_cxx14::make_unique<VectorType>());

  return get_pool().data->back().second.get();
}



template <typename VectorType>
inline void
GrowingVectorMemory<VectorType>::free(const VectorType *const v)
{
  std::lock_guard<std::mutex> lock(mutex);

  for (typename std::vector<entry_type>::iterator i = get_pool().data->begin();
       i != get_pool().data->end();
       ++i)
    {
      if (v == i->second.get())
        {
          i->first = false;
          --current_alloc;
          return;
        }
    }
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

  std::size_t                                            result = sizeof(*this);
  const typename std::vector<entry_type>::const_iterator end =
    get_pool().data->end();
  for (typename std::vector<entry_type>::const_iterator i =
         get_pool().data->begin();
       i != end;
       ++i)
    result += sizeof(*i) + MemoryConsumption::memory_consumption(i->second);

  return result;
}


DEAL_II_NAMESPACE_CLOSE

#endif
