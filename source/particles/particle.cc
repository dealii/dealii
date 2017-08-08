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

#include <deal.II/particles/particle.h>

DEAL_II_NAMESPACE_OPEN

namespace Particles
{
  template <int dim>
  Particle<dim>::Particle ()
    :
    location (),
    reference_location(),
    id (0),
    property_pool(NULL),
    properties(PropertyPool::invalid_handle)
  {
  }


  template <int dim>
  Particle<dim>::Particle (const Point<dim> &new_location,
                           const Point<dim> &new_reference_location,
                           const types::particle_index new_id)
    :
    location (new_location),
    reference_location (new_reference_location),
    id (new_id),
    property_pool(NULL),
    properties (PropertyPool::invalid_handle)
  {
  }

  template <int dim>
  Particle<dim>::Particle (const Particle<dim> &particle)
    :
    location (particle.get_location()),
    reference_location (particle.get_reference_location()),
    id (particle.get_id()),
    property_pool(particle.property_pool),
    properties ((property_pool != NULL) ? property_pool->allocate_properties_array() : PropertyPool::invalid_handle)
  {
    if (property_pool != NULL)
      {
        const ArrayView<double> my_properties = property_pool->get_properties(properties);

        if (my_properties.size() != 0)
          {
            const ArrayView<const double> their_properties = particle.get_properties();

            std::copy(&their_properties[0],&their_properties[0]+their_properties.size(),&my_properties[0]);
          }
      }
  }


  template <int dim>
  Particle<dim>::Particle (const void *&data,
                           PropertyPool &new_property_pool)
  {
    const types::particle_index *id_data = static_cast<const types::particle_index *> (data);
    id = *id_data++;
    const double *pdata = reinterpret_cast<const double *> (id_data);

    for (unsigned int i = 0; i < dim; ++i)
      location(i) = *pdata++;

    for (unsigned int i = 0; i < dim; ++i)
      reference_location(i) = *pdata++;

    property_pool = &new_property_pool;
    properties = property_pool->allocate_properties_array();

    // See if there are properties to load
    const ArrayView<double> particle_properties = property_pool->get_properties(properties);
    for (unsigned int i = 0; i < particle_properties.size(); ++i)
      particle_properties[i] = *pdata++;

    data = static_cast<const void *> (pdata);
  }

  template <int dim>
  Particle<dim>::Particle (Particle<dim> &&particle)
    :
    location (particle.location),
    reference_location(particle.reference_location),
    id (particle.id),
    property_pool(particle.property_pool),
    properties(particle.properties)
  {
    particle.properties = PropertyPool::invalid_handle;
  }

  template <int dim>
  Particle<dim> &
  Particle<dim>::operator=(const Particle<dim> &particle)
  {
    if (this != &particle)
      {
        location = particle.location;
        reference_location = particle.reference_location;
        id = particle.id;
        property_pool = particle.property_pool;

        if (property_pool != NULL)
          {
            properties = property_pool->allocate_properties_array();
            const ArrayView<const double> their_properties = particle.get_properties();

            if (their_properties.size() != 0)
              {
                const ArrayView<double> my_properties = property_pool->get_properties(properties);
                std::copy(&their_properties[0],&their_properties[0]+their_properties.size(),&my_properties[0]);
              }
          }
        else
          properties = PropertyPool::invalid_handle;
      }
    return *this;
  }

  template <int dim>
  Particle<dim> &
  Particle<dim>::operator=(Particle<dim> &&particle)
  {
    if (this != &particle)
      {
        location = particle.location;
        reference_location = particle.reference_location;
        id = particle.id;
        property_pool = particle.property_pool;
        properties = particle.properties;
        particle.properties = PropertyPool::invalid_handle;
      }
    return *this;
  }

  template <int dim>
  Particle<dim>::~Particle ()
  {
    if (properties != PropertyPool::invalid_handle)
      property_pool->deallocate_properties_array(properties);
  }

  template <int dim>
  void
  Particle<dim>::write_data (void *&data) const
  {
    types::particle_index *id_data  = static_cast<types::particle_index *> (data);
    *id_data = id;
    ++id_data;
    double *pdata = reinterpret_cast<double *> (id_data);

    // Write location data
    for (unsigned int i = 0; i < dim; ++i,++pdata)
      *pdata = location(i);

    // Write reference location data
    for (unsigned int i = 0; i < dim; ++i,++pdata)
      *pdata = reference_location(i);

    // Write property data
    const ArrayView<double> particle_properties = property_pool->get_properties(properties);
    for (unsigned int i = 0; i < particle_properties.size(); ++i,++pdata)
      *pdata = particle_properties[i];

    data = static_cast<void *> (pdata);
  }

  template <int dim>
  void
  Particle<dim>::set_location (const Point<dim> &new_loc)
  {
    location = new_loc;
  }

  template <int dim>
  const Point<dim> &
  Particle<dim>::get_location () const
  {
    return location;
  }

  template <int dim>
  void
  Particle<dim>::set_reference_location (const Point<dim> &new_loc)
  {
    reference_location = new_loc;
  }

  template <int dim>
  const Point<dim> &
  Particle<dim>::get_reference_location () const
  {
    return reference_location;
  }

  template <int dim>
  types::particle_index
  Particle<dim>::get_id () const
  {
    return id;
  }

  template <int dim>
  void
  Particle<dim>::set_property_pool (PropertyPool &new_property_pool)
  {
    property_pool = &new_property_pool;
  }

  template <int dim>
  void
  Particle<dim>::set_properties (const std::vector<double> &new_properties)
  {
    if (properties == PropertyPool::invalid_handle)
      properties = property_pool->allocate_properties_array();

    const ArrayView<double> old_properties = property_pool->get_properties(properties);

    std::copy(new_properties.begin(),new_properties.end(),&old_properties[0]);
  }

  template <int dim>
  const ArrayView<const double>
  Particle<dim>::get_properties () const
  {
    Assert(property_pool != NULL,
           ExcInternalError());

    return property_pool->get_properties(properties);
  }


  template <int dim>
  const ArrayView<double>
  Particle<dim>::get_properties ()
  {
    Assert(property_pool != NULL,
           ExcInternalError());

    return property_pool->get_properties(properties);
  }
}

DEAL_II_NAMESPACE_CLOSE

DEAL_II_NAMESPACE_OPEN

#include "particle.inst"

DEAL_II_NAMESPACE_CLOSE
