// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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


#include <deal.II/particles/property_pool.h>

DEAL_II_NAMESPACE_OPEN

namespace Particles
{
  const PropertyPool::Handle PropertyPool::invalid_handle =
    static_cast<Handle>(-1);


  PropertyPool::PropertyPool(const unsigned int n_properties_per_slot)
    : n_properties(n_properties_per_slot)
  {}


  PropertyPool::~PropertyPool()
  {
    clear();
  }



  void
  PropertyPool::clear()
  {
    if (n_properties > 0)
      {
        const unsigned int n_open_handles =
          properties.size() / n_properties - currently_available_handles.size();
        (void)n_open_handles;
        AssertThrow(n_open_handles == 0,
                    ExcMessage("This property pool currently still holds " +
                               std::to_string(n_open_handles) +
                               " open handles to memory that was allocated "
                               "via allocate_properties_array() but that has "
                               "not been returned via "
                               "deallocate_properties_array()."));
      }

    // Clear vectors and ensure deallocation of memory
    properties.clear();
    properties.shrink_to_fit();

    currently_available_handles.clear();
    currently_available_handles.shrink_to_fit();
  }



  PropertyPool::Handle
  PropertyPool::allocate_properties_array()
  {
    Handle handle = invalid_handle;
    if (n_properties > 0)
      {
        if (currently_available_handles.size() > 0)
          {
            handle = currently_available_handles.back();
            currently_available_handles.pop_back();
          }
        else
          {
            handle = properties.size() / n_properties;
            properties.resize(properties.size() + n_properties, 0.0);
          }
      }

    return handle;
  }



  void
  PropertyPool::deallocate_properties_array(Handle &handle)
  {
    Assert(
      handle != invalid_handle,
      ExcMessage(
        "This handle is invalid and cannot be deallocated. This can happen if the "
        "handle was deallocated already before calling this function."));
    currently_available_handles.push_back(handle);
    handle = invalid_handle;

    // If this was the last handle, resize containers to 0.
    // This guarantees that all future properties
    // are allocated in a sorted way without requiring reallocation.
    if (currently_available_handles.size() * n_properties == properties.size())
      {
        currently_available_handles.clear();
        properties.clear();
      }
  }



  void
  PropertyPool::reserve(const std::size_t size)
  {
    properties.reserve(size * n_properties);
  }



  unsigned int
  PropertyPool::n_properties_per_slot() const
  {
    return n_properties;
  }
} // namespace Particles
DEAL_II_NAMESPACE_CLOSE
