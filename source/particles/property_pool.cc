// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/array_view.h>
#include <deal.II/base/memory_space.h>
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
    const unsigned int n_open_handles =
      locations.size() - currently_available_handles.size();
    (void)n_open_handles;
    Assert(n_open_handles == 0,
           ExcMessage("This property pool currently still holds " +
                      std::to_string(n_open_handles) +
                      " open handles to memory that was allocated "
                      "via allocate_properties_array() but that has "
                      "not been returned via "
                      "deregister_particle()."));

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

    Assert(
      currently_available_handles.size() < locations.size(),
      ExcMessage(
        "Trying to deallocate a particle when none are allocated. This can happen if all "
        "handles were deallocated already before calling this function."));

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



  template <int dim, int spacedim>
  unsigned int
  PropertyPool<dim, spacedim>::n_registered_slots() const
  {
    Assert(locations.size() == reference_locations.size(),
           ExcMessage("Number of registered locations is not equal to number "
                      "of registered reference locations."));

    Assert(locations.size() == ids.size(),
           ExcMessage("Number of registered locations is not equal to number "
                      "of registered ids."));

    Assert(locations.size() * n_properties == properties.size(),
           ExcMessage("Number of registered locations is not equal to number "
                      "of registered property slots."));

    return locations.size() - currently_available_handles.size();
  }



  template <int dim, int spacedim>
  void
  PropertyPool<dim, spacedim>::sort_memory_slots(
    const std::vector<Handle> &handles_to_sort)
  {
    std::vector<Point<spacedim>>       sorted_locations;
    std::vector<Point<dim>>            sorted_reference_locations;
    std::vector<types::particle_index> sorted_ids;
    std::vector<double>                sorted_properties;

    sorted_locations.reserve(locations.size());
    sorted_reference_locations.reserve(reference_locations.size());
    sorted_ids.reserve(ids.size());
    sorted_properties.reserve(properties.size());

    for (const auto &handle : handles_to_sort)
      {
        Assert(handle != invalid_handle,
               ExcMessage(
                 "Invalid handle detected during sorting particle memory."));

        sorted_locations.push_back(locations[handle]);
        sorted_reference_locations.push_back(reference_locations[handle]);
        sorted_ids.push_back(ids[handle]);

        for (unsigned int j = 0; j < n_properties; ++j)
          sorted_properties.push_back(properties[handle * n_properties + j]);
      }

    Assert(sorted_locations.size() ==
             locations.size() - currently_available_handles.size(),
           ExcMessage("Number of sorted property handles is not equal to "
                      "number of currently registered handles: " +
                      std::to_string(sorted_locations.size()) + " vs " +
                      std::to_string(locations.size()) + " - " +
                      std::to_string(currently_available_handles.size())));

    locations           = std::move(sorted_locations);
    reference_locations = std::move(sorted_reference_locations);
    ids                 = std::move(sorted_ids);
    properties          = std::move(sorted_properties);

    currently_available_handles.clear();
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
