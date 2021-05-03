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


#include <deal.II/base/signaling_nan.h>

#include <deal.II/particles/property_pool.h>

DEAL_II_NAMESPACE_OPEN

namespace Particles
{
  template <int dim, int spacedim>
  const typename PropertyPool<dim, spacedim>::Handle
    PropertyPool<dim, spacedim>::invalid_handle = static_cast<Handle>(-1);



  template <int dim, int spacedim>
  PropertyPool<dim, spacedim>::PropertyPool(
    const unsigned int n_properties_per_slot)
    : n_properties(n_properties_per_slot)
  {}



  template <int dim, int spacedim>
  PropertyPool<dim, spacedim>::~PropertyPool()
  {
    clear();
  }



  template <int dim, int spacedim>
  void
  PropertyPool<dim, spacedim>::clear()
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
    locations.clear();
    locations.shrink_to_fit();

    reference_locations.clear();
    reference_locations.shrink_to_fit();

    ids.clear();
    ids.shrink_to_fit();

    properties.clear();
    properties.shrink_to_fit();

    currently_available_handles.clear();
    currently_available_handles.shrink_to_fit();
  }



  template <int dim, int spacedim>
  typename PropertyPool<dim, spacedim>::Handle
  PropertyPool<dim, spacedim>::register_particle()
  {
    Handle handle = invalid_handle;
    if (currently_available_handles.size() > 0)
      {
        handle = currently_available_handles.back();
        currently_available_handles.pop_back();
      }
    else
      {
        handle = locations.size();

        locations.resize(locations.size() + 1);
        reference_locations.resize(reference_locations.size() + 1);
        ids.resize(ids.size() + 1);
        properties.resize(properties.size() + n_properties);
      }

    // Then initialize whatever slot we have taken with invalid locations,
    // but initialize properties with zero.
    set_location(handle, numbers::signaling_nan<Point<spacedim>>());
    set_reference_location(handle, numbers::signaling_nan<Point<dim>>());
    set_id(handle, numbers::invalid_unsigned_int);
    for (double &x : get_properties(handle))
      x = 0;

    return handle;
  }



  template <int dim, int spacedim>
  void
  PropertyPool<dim, spacedim>::deregister_particle(Handle &handle)
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
    if (currently_available_handles.size() == locations.size())
      {
        currently_available_handles.clear();
        properties.clear();
        locations.clear();
        reference_locations.clear();
        ids.clear();
      }
  }



  template <int dim, int spacedim>
  void
  PropertyPool<dim, spacedim>::reserve(const std::size_t size)
  {
    locations.reserve(size);
    reference_locations.reserve(size);
    properties.reserve(size * n_properties);
    ids.reserve(size);
  }



  template <int dim, int spacedim>
  unsigned int
  PropertyPool<dim, spacedim>::n_properties_per_slot() const
  {
    return n_properties;
  }


  // Instantiate the class for all reasonable template arguments
  template class PropertyPool<1, 1>;
  template class PropertyPool<1, 2>;
  template class PropertyPool<1, 3>;
  template class PropertyPool<2, 2>;
  template class PropertyPool<2, 3>;
  template class PropertyPool<3, 3>;
} // namespace Particles
DEAL_II_NAMESPACE_CLOSE
