// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#include <deal.II/particles/property_pool.h>

DEAL_II_NAMESPACE_OPEN

namespace Particles
{
  const PropertyPool::Handle PropertyPool::invalid_handle = nullptr;


  PropertyPool::PropertyPool (const unsigned int n_properties_per_slot)
    :
    n_properties (n_properties_per_slot)
  {}



  PropertyPool::Handle
  PropertyPool::allocate_properties_array ()
  {
    return new double[n_properties];
  }



  void
  PropertyPool::deallocate_properties_array (Handle handle)
  {
    delete[] handle;
  }



  ArrayView<double>
  PropertyPool::get_properties (const Handle handle)
  {
    return ArrayView<double>(handle, n_properties);
  }



  void
  PropertyPool::reserve(const std::size_t size)
  {
    (void)size;
  }



  unsigned int
  PropertyPool::n_properties_per_slot() const
  {
    return n_properties;
  }
}
DEAL_II_NAMESPACE_CLOSE
