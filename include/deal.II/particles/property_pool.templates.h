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

#ifndef dealii__particles_property_pool_templates_h
#define dealii__particles_property_pool_templates_h

#include <deal.II/particles/property_pool.h>
#include <deal.II/base/exceptions.h>

DEAL_II_NAMESPACE_OPEN

namespace Particles
{
  template <typename PropertyType>
  const typename PropertyPool<PropertyType>::Handle PropertyPool<PropertyType>::invalid_handle = NULL;



  template <typename PropertyType>
  PropertyPool<PropertyType>::PropertyPool (const unsigned int n_properties_per_slot)
    :
    n_properties (n_properties_per_slot)
  {}



  template <typename PropertyType>
  typename PropertyPool<PropertyType>::Handle
  PropertyPool<PropertyType>::allocate_properties_array ()
  {
    return new PropertyType[n_properties];
  }


  template <typename PropertyType>
  void
  PropertyPool<PropertyType>::deallocate_properties_array (Handle &handle)
  {
    delete[] handle;
    handle = invalid_handle;
  }



  template <typename PropertyType>
  const ArrayView<PropertyType>
  PropertyPool<PropertyType>::get_properties (const Handle handle)
  {
    Assert(handle != PropertyPool<PropertyType>::invalid_handle,
           ExcMessage("PropertyPool was asked for the properties that belong to an invalid handle. "
                      "Either this handle was never validly created, or it was deallocated before use."));

    return ArrayView<PropertyType>(handle, n_properties);
  }



  template <typename PropertyType>
  const ArrayView<const PropertyType>
  PropertyPool<PropertyType>::get_properties (const Handle handle) const
  {
    Assert(handle != PropertyPool<PropertyType>::invalid_handle,
           ExcMessage("PropertyPool was asked for the properties that belong to an invalid handle. "
                      "Either this handle was never validly created, or it was deallocated before use."));

    return ArrayView<const PropertyType>(handle, n_properties);
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
