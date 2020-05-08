// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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

#include <deal.II/particles/particle.h>

DEAL_II_NAMESPACE_OPEN

namespace Particles
{
  template <int dim, int spacedim>
  Particle<dim, spacedim>::Particle()
    : location(numbers::signaling_nan<Point<spacedim>>())
    , reference_location(numbers::signaling_nan<Point<dim>>())
    , id(0)
    , property_pool(nullptr)
    , properties(PropertyPool::invalid_handle)
  {}



  template <int dim, int spacedim>
  Particle<dim, spacedim>::Particle(const Point<spacedim> &location,
                                    const Point<dim> &     reference_location,
                                    const types::particle_index id)
    : location(location)
    , reference_location(reference_location)
    , id(id)
    , property_pool(nullptr)
    , properties(PropertyPool::invalid_handle)
  {}



  template <int dim, int spacedim>
  Particle<dim, spacedim>::Particle(const Particle<dim, spacedim> &particle)
    : location(particle.get_location())
    , reference_location(particle.get_reference_location())
    , id(particle.get_id())
    , property_pool(particle.property_pool)
    , properties((particle.has_properties()) ?
                   property_pool->allocate_properties_array() :
                   PropertyPool::invalid_handle)
  {
    if (particle.has_properties())
      {
        const ArrayView<double> my_properties =
          property_pool->get_properties(properties);
        const ArrayView<const double> their_properties =
          particle.get_properties();

        std::copy(their_properties.begin(),
                  their_properties.end(),
                  my_properties.begin());
      }
  }



  template <int dim, int spacedim>
  Particle<dim, spacedim>::Particle(const void *&       data,
                                    PropertyPool *const new_property_pool)
  {
    const types::particle_index *id_data =
      static_cast<const types::particle_index *>(data);
    id                  = *id_data++;
    const double *pdata = reinterpret_cast<const double *>(id_data);

    for (unsigned int i = 0; i < spacedim; ++i)
      location(i) = *pdata++;

    for (unsigned int i = 0; i < dim; ++i)
      reference_location(i) = *pdata++;

    property_pool = new_property_pool;
    if (property_pool != nullptr)
      properties = property_pool->allocate_properties_array();
    else
      properties = PropertyPool::invalid_handle;

    // See if there are properties to load
    if (has_properties())
      {
        const ArrayView<double> particle_properties =
          property_pool->get_properties(properties);
        const unsigned int size = particle_properties.size();
        for (unsigned int i = 0; i < size; ++i)
          particle_properties[i] = *pdata++;
      }

    data = static_cast<const void *>(pdata);
  }



  template <int dim, int spacedim>
  Particle<dim, spacedim>::Particle(Particle<dim, spacedim> &&particle) noexcept
    : location(std::move(particle.location))
    , reference_location(std::move(particle.reference_location))
    , id(std::move(particle.id))
    , property_pool(std::move(particle.property_pool))
    , properties(std::move(particle.properties))
  {
    particle.properties = PropertyPool::invalid_handle;
  }



  template <int dim, int spacedim>
  Particle<dim, spacedim> &
  Particle<dim, spacedim>::operator=(const Particle<dim, spacedim> &particle)
  {
    if (this != &particle)
      {
        location           = particle.location;
        reference_location = particle.reference_location;
        id                 = particle.id;
        property_pool      = particle.property_pool;

        if (particle.has_properties())
          {
            properties = property_pool->allocate_properties_array();
            const ArrayView<const double> their_properties =
              particle.get_properties();
            const ArrayView<double> my_properties =
              property_pool->get_properties(properties);

            std::copy(their_properties.begin(),
                      their_properties.end(),
                      my_properties.begin());
          }
        else
          properties = PropertyPool::invalid_handle;
      }
    return *this;
  }



  template <int dim, int spacedim>
  Particle<dim, spacedim> &
  Particle<dim, spacedim>::
  operator=(Particle<dim, spacedim> &&particle) noexcept
  {
    if (this != &particle)
      {
        location            = particle.location;
        reference_location  = particle.reference_location;
        id                  = particle.id;
        property_pool       = particle.property_pool;
        properties          = particle.properties;
        particle.properties = PropertyPool::invalid_handle;
      }
    return *this;
  }



  template <int dim, int spacedim>
  Particle<dim, spacedim>::~Particle()
  {
    if (property_pool != nullptr && properties != PropertyPool::invalid_handle)
      property_pool->deallocate_properties_array(properties);
  }



  template <int dim, int spacedim>
  void
  Particle<dim, spacedim>::write_data(void *&data) const
  {
    types::particle_index *id_data = static_cast<types::particle_index *>(data);
    *id_data                       = id;
    ++id_data;
    double *pdata = reinterpret_cast<double *>(id_data);

    // Write location data
    for (unsigned int i = 0; i < spacedim; ++i, ++pdata)
      *pdata = location(i);

    // Write reference location data
    for (unsigned int i = 0; i < dim; ++i, ++pdata)
      *pdata = reference_location(i);

    // Write property data
    if (has_properties())
      {
        const ArrayView<double> particle_properties =
          property_pool->get_properties(properties);
        for (unsigned int i = 0; i < particle_properties.size(); ++i, ++pdata)
          *pdata = particle_properties[i];
      }

    data = static_cast<void *>(pdata);
  }



  template <int dim, int spacedim>
  std::size_t
  Particle<dim, spacedim>::serialized_size_in_bytes() const
  {
    std::size_t size = sizeof(types::particle_index) + sizeof(location) +
                       sizeof(reference_location);

    if (has_properties())
      {
        const ArrayView<double> particle_properties =
          property_pool->get_properties(properties);
        size += sizeof(double) * particle_properties.size();
      }
    return size;
  }



  template <int dim, int spacedim>
  void
  Particle<dim, spacedim>::set_location(const Point<spacedim> &new_loc)
  {
    location = new_loc;
  }



  template <int dim, int spacedim>
  const Point<spacedim> &
  Particle<dim, spacedim>::get_location() const
  {
    return location;
  }



  template <int dim, int spacedim>
  void
  Particle<dim, spacedim>::set_reference_location(const Point<dim> &new_loc)
  {
    reference_location = new_loc;
  }



  template <int dim, int spacedim>
  const Point<dim> &
  Particle<dim, spacedim>::get_reference_location() const
  {
    return reference_location;
  }



  template <int dim, int spacedim>
  types::particle_index
  Particle<dim, spacedim>::get_id() const
  {
    return id;
  }



  template <int dim, int spacedim>
  void
  Particle<dim, spacedim>::set_id(const types::particle_index &new_id)
  {
    id = new_id;
  }



  template <int dim, int spacedim>
  void
  Particle<dim, spacedim>::set_property_pool(PropertyPool &new_property_pool)
  {
    property_pool = &new_property_pool;
  }



  template <int dim, int spacedim>
  void
  Particle<dim, spacedim>::set_properties(
    const ArrayView<const double> &new_properties)
  {
    Assert(property_pool != nullptr, ExcInternalError());

    if (properties == PropertyPool::invalid_handle)
      properties = property_pool->allocate_properties_array();

    const ArrayView<double> old_properties =
      property_pool->get_properties(properties);

    Assert(
      new_properties.size() == old_properties.size(),
      ExcMessage(
        std::string(
          "You are trying to assign properties with an incompatible length. ") +
        "The particle has space to store " +
        std::to_string(old_properties.size()) + " properties, " +
        "and this function tries to assign" +
        std::to_string(new_properties.size()) + " properties. " +
        "This is not allowed."));

    if (old_properties.size() > 0)
      std::copy(new_properties.begin(),
                new_properties.end(),
                old_properties.begin());
  }



  template <int dim, int spacedim>
  const ArrayView<const double>
  Particle<dim, spacedim>::get_properties() const
  {
    Assert(has_properties(), ExcInternalError());

    return property_pool->get_properties(properties);
  }



  template <int dim, int spacedim>
  const ArrayView<double>
  Particle<dim, spacedim>::get_properties()
  {
    Assert(property_pool != nullptr, ExcInternalError());

    // If this particle has no properties yet, allocate and initialize them.
    if (properties == PropertyPool::invalid_handle)
      {
        properties = property_pool->allocate_properties_array();

        ArrayView<double> my_properties =
          property_pool->get_properties(properties);

        std::fill(my_properties.begin(), my_properties.end(), 0.0);
      }

    return property_pool->get_properties(properties);
  }



  template <int dim, int spacedim>
  bool
  Particle<dim, spacedim>::has_properties() const
  {
    return (property_pool != nullptr) &&
           (properties != PropertyPool::invalid_handle);
  }
} // namespace Particles

DEAL_II_NAMESPACE_CLOSE

DEAL_II_NAMESPACE_OPEN

#include "particle.inst"

DEAL_II_NAMESPACE_CLOSE
