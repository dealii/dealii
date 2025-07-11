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

#include <deal.II/base/signaling_nan.h>

#include <deal.II/particles/particle.h>

DEAL_II_NAMESPACE_OPEN

namespace Particles
{
  template <int dim, int spacedim>
  PropertyPool<dim, spacedim> Particle<dim, spacedim>::global_property_pool = {
    /* no properties: */ 0};


  template <int dim, int spacedim>
  Particle<dim, spacedim>::Particle()
    : property_pool(&global_property_pool)
    , property_pool_handle(property_pool->register_particle())
  {}



  template <int dim, int spacedim>
  Particle<dim, spacedim>::Particle(const Point<spacedim> &location,
                                    const Point<dim>      &reference_location,
                                    const types::particle_index id)
    : property_pool(&global_property_pool)
    , property_pool_handle(property_pool->register_particle())
  {
    set_location(location);
    set_reference_location(reference_location);
    set_id(id);
  }



  template <int dim, int spacedim>
  Particle<dim, spacedim>::Particle(const Particle<dim, spacedim> &particle)
    : property_pool(particle.property_pool)
    , property_pool_handle(property_pool->register_particle())
  {
    set_location(particle.get_location());
    set_reference_location(particle.get_reference_location());
    set_id(particle.get_id());

    if (particle.has_properties())
      {
        const ArrayView<const double> their_properties =
          particle.get_properties();
        const ArrayView<double> my_properties =
          property_pool->get_properties(property_pool_handle);

        std::copy(their_properties.begin(),
                  their_properties.end(),
                  my_properties.begin());
      }
  }



  template <int dim, int spacedim>
  Particle<dim, spacedim>::Particle(
    const void                       *&data,
    PropertyPool<dim, spacedim> *const new_property_pool)
    : property_pool(new_property_pool != nullptr ? new_property_pool :
                                                   &global_property_pool)
    , property_pool_handle(property_pool->register_particle())
  {
    const types::particle_index *id_data =
      static_cast<const types::particle_index *>(data);
    set_id(*id_data++);
    const double *pdata = reinterpret_cast<const double *>(id_data);

    Point<spacedim> location;
    for (unsigned int i = 0; i < spacedim; ++i)
      location[i] = *pdata++;
    set_location(location);

    Point<dim> reference_location;
    for (unsigned int i = 0; i < dim; ++i)
      reference_location[i] = *pdata++;
    set_reference_location(reference_location);

    // See if there are properties to load
    if (has_properties())
      {
        const ArrayView<double> particle_properties = this->get_properties();
        const unsigned int      size = particle_properties.size();
        for (unsigned int i = 0; i < size; ++i)
          particle_properties[i] = *pdata++;
      }

    data = static_cast<const void *>(pdata);
  }



  template <int dim, int spacedim>
  Particle<dim, spacedim>::Particle(Particle<dim, spacedim> &&particle) noexcept
    : property_pool(std::move(particle.property_pool))
    , property_pool_handle(std::move(particle.property_pool_handle))
  {
    // There is no need to copy locations, properties, and id -- we simply
    // inherit them from the object we move from.

    // We stole the rhs's properties, so we need to invalidate
    // the handle the rhs holds lest it releases the memory that
    // we still reference here.
    particle.property_pool_handle = PropertyPool<dim, spacedim>::invalid_handle;
  }



  template <int dim, int spacedim>
  Particle<dim, spacedim> &
  Particle<dim, spacedim>::operator=(const Particle<dim, spacedim> &particle)
  {
    if (this != &particle)
      {
        Assert(this->property_pool->n_properties_per_slot() ==
                 particle.property_pool->n_properties_per_slot(),
               ExcInternalError());

        set_location(particle.get_location());
        set_reference_location(particle.get_reference_location());
        set_id(particle.get_id());

        if (particle.has_properties())
          {
            const ArrayView<const double> their_properties =
              particle.get_properties();
            const ArrayView<double> my_properties =
              property_pool->get_properties(property_pool_handle);

            std::copy(their_properties.begin(),
                      their_properties.end(),
                      my_properties.begin());
          }
      }

    return *this;
  }



  template <int dim, int spacedim>
  Particle<dim, spacedim> &
  Particle<dim, spacedim>::operator=(
    Particle<dim, spacedim> &&particle) noexcept
  {
    if (this != &particle)
      {
        // If we currently hold a handle, release the memory. (The only way a
        // particle can end up not holding a valid handle is if it has been
        // moved from.)
        if (property_pool_handle != PropertyPool<dim, spacedim>::invalid_handle)
          property_pool->deregister_particle(property_pool_handle);

        property_pool        = std::move(particle.property_pool);
        property_pool_handle = std::move(particle.property_pool_handle);

        // No need to copy locations, properties, and id -- we just get them
        // by taking over the property pool handle.

        // We stole the rhs's properties, so we need to invalidate
        // the handle the rhs holds lest it releases the memory that
        // we still reference here.
        particle.property_pool_handle =
          PropertyPool<dim, spacedim>::invalid_handle;
      }
    return *this;
  }



  template <int dim, int spacedim>
  Particle<dim, spacedim>::~Particle()
  {
    // If we still hold a handle, release the memory. The only way a
    // particle can end up not holding a valid handle is if it has been
    // moved from.
    if (property_pool_handle != PropertyPool<dim, spacedim>::invalid_handle)
      property_pool->deregister_particle(property_pool_handle);
  }



  template <int dim, int spacedim>
  void
  Particle<dim, spacedim>::free_properties()
  {
    // If we still hold a handle, release the memory. The only way a
    // particle can end up not holding a valid handle is if it has been
    // moved from.
    //
    // deregister_particle() automatically invalidates its argument.
    if (property_pool_handle != PropertyPool<dim, spacedim>::invalid_handle)
      property_pool->deregister_particle(property_pool_handle);
  }



  template <int dim, int spacedim>
  void *
  Particle<dim, spacedim>::write_particle_data_to_memory(
    void *data_pointer) const
  {
    types::particle_index *id_data =
      static_cast<types::particle_index *>(data_pointer);
    *id_data = get_id();
    ++id_data;
    double *pdata = reinterpret_cast<double *>(id_data);

    // Write location
    for (unsigned int i = 0; i < spacedim; ++i, ++pdata)
      *pdata = get_location()[i];

    // Write reference location
    for (unsigned int i = 0; i < dim; ++i, ++pdata)
      *pdata = get_reference_location()[i];

    // Write properties
    if (has_properties())
      {
        const ArrayView<double> particle_properties =
          property_pool->get_properties(property_pool_handle);
        for (unsigned int i = 0; i < particle_properties.size(); ++i, ++pdata)
          *pdata = particle_properties[i];
      }

    return static_cast<void *>(pdata);
  }



  template <int dim, int spacedim>
  const void *
  Particle<dim, spacedim>::read_particle_data_from_memory(
    const void *data_pointer)
  {
    const types::particle_index *id_data =
      static_cast<const types::particle_index *>(data_pointer);
    set_id(*id_data++);
    const double *pdata = reinterpret_cast<const double *>(id_data);

    Point<spacedim> location;
    for (unsigned int i = 0; i < spacedim; ++i)
      location[i] = *pdata++;
    set_location(location);

    Point<dim> reference_location;
    for (unsigned int i = 0; i < dim; ++i)
      reference_location[i] = *pdata++;
    set_reference_location(reference_location);

    // See if there are properties to load
    if (has_properties())
      {
        const ArrayView<double> particle_properties =
          property_pool->get_properties(property_pool_handle);
        const unsigned int size = particle_properties.size();
        for (unsigned int i = 0; i < size; ++i)
          particle_properties[i] = *pdata++;
      }

    return static_cast<const void *>(pdata);
  }



  template <int dim, int spacedim>
  std::size_t
  Particle<dim, spacedim>::serialized_size_in_bytes() const
  {
    std::size_t size = sizeof(get_id()) +
                       sizeof(double) * spacedim + // get_location()
                       sizeof(double) * dim;       // get_reference_location()

    if (has_properties())
      {
        const ArrayView<double> particle_properties =
          property_pool->get_properties(property_pool_handle);
        size += sizeof(double) * particle_properties.size();
      }
    return size;
  }



  template <int dim, int spacedim>
  void
  Particle<dim, spacedim>::set_properties(
    const ArrayView<const double> &new_properties)
  {
    const ArrayView<double> property_values =
      property_pool->get_properties(property_pool_handle);

    Assert(new_properties.size() == property_values.size(),
           ExcMessage(
             "You are trying to assign properties with an incompatible length. "
             "The particle has space to store " +
             std::to_string(property_values.size()) +
             " properties, but you are trying to assign " +
             std::to_string(new_properties.size()) +
             " properties. This is not allowed."));

    if (property_values.size() > 0)
      std::copy(new_properties.begin(),
                new_properties.end(),
                property_values.begin());
  }



  template <int dim, int spacedim>
  ArrayView<double>
  Particle<dim, spacedim>::get_properties()
  {
    return property_pool->get_properties(property_pool_handle);
  }
} // namespace Particles

#include "particles/particle.inst"

DEAL_II_NAMESPACE_CLOSE
