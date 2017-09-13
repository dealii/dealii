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
  template <int dim, int spacedim, typename PropertyType>
  Particle<dim,spacedim,PropertyType>::Particle ()
    :
    location (),
    reference_location(),
    id (0),
    property_pool(NULL),
    properties(PropertyPool<PropertyType>::invalid_handle)
  {
  }



  template <int dim, int spacedim, typename PropertyType>
  Particle<dim,spacedim,PropertyType>::Particle (const Point<spacedim> &location,
                                                 const Point<dim> &reference_location,
                                                 const types::particle_index id)
    :
    location (location),
    reference_location (reference_location),
    id (id),
    property_pool(NULL),
    properties (PropertyPool<PropertyType>::invalid_handle)
  {
  }



  template <int dim, int spacedim, typename PropertyType>
  Particle<dim,spacedim,PropertyType>::Particle (const Particle<dim,spacedim,PropertyType> &particle)
    :
    location (particle.get_location()),
    reference_location (particle.get_reference_location()),
    id (particle.get_id()),
    property_pool(particle.property_pool),
    properties ((property_pool != NULL) ? property_pool->allocate_properties_array() : PropertyPool<PropertyType>::invalid_handle)
  {
    if (property_pool != NULL)
      {
        const ArrayView<PropertyType> my_properties = property_pool->get_properties(properties);

        if (my_properties.size() != 0)
          {
            const ArrayView<const PropertyType> their_properties = particle.get_properties();

            std::copy(&their_properties[0],&their_properties[0]+their_properties.size(),&my_properties[0]);
          }
      }
  }



  template <int dim, int spacedim, typename PropertyType>
  Particle<dim,spacedim,PropertyType>::Particle (const void *&data,
                                                 PropertyPool<PropertyType> *new_property_pool)
  {
    const types::particle_index *id_data = static_cast<const types::particle_index *> (data);
    id = *id_data++;
    const double *pdata = reinterpret_cast<const double *> (id_data);

    for (unsigned int i = 0; i < dim; ++i)
      location(i) = *pdata++;

    for (unsigned int i = 0; i < dim; ++i)
      reference_location(i) = *pdata++;

    data = static_cast<const void *> (pdata);

    // See if there are properties to load
    if (new_property_pool != NULL)
      {
        property_pool = new_property_pool;
        properties = property_pool->allocate_properties_array();

        const PropertyType *property_data = static_cast<const PropertyType *> (data);

        const ArrayView<PropertyType> particle_properties = property_pool->get_properties(properties);
        for (unsigned int i = 0; i < particle_properties.size(); ++i)
          particle_properties[i] = *property_data++;

        data = static_cast<const void *> (property_data);
      }
  }



  template <int dim, int spacedim, typename PropertyType>
  Particle<dim,spacedim,PropertyType>::Particle (Particle<dim,spacedim,PropertyType> &&particle)
    :
    location (particle.location),
    reference_location(particle.reference_location),
    id (particle.id),
    property_pool(particle.property_pool),
    properties(particle.properties)
  {
    particle.properties = PropertyPool<PropertyType>::invalid_handle;
  }



  template <int dim, int spacedim, typename PropertyType>
  Particle<dim,spacedim,PropertyType> &
  Particle<dim,spacedim,PropertyType>::operator=(const Particle<dim,spacedim,PropertyType> &particle)
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
            const ArrayView<const PropertyType> their_properties = particle.get_properties();

            if (their_properties.size() != 0)
              {
                const ArrayView<PropertyType> my_properties = property_pool->get_properties(properties);
                std::copy(&their_properties[0],&their_properties[0]+their_properties.size(),&my_properties[0]);
              }
          }
        else
          properties = PropertyPool<PropertyType>::invalid_handle;
      }
    return *this;
  }



  template <int dim, int spacedim, typename PropertyType>
  Particle<dim,spacedim,PropertyType> &
  Particle<dim,spacedim,PropertyType>::operator=(Particle<dim,spacedim,PropertyType> &&particle)
  {
    if (this != &particle)
      {
        location = particle.location;
        reference_location = particle.reference_location;
        id = particle.id;
        property_pool = particle.property_pool;
        properties = particle.properties;
        particle.properties = PropertyPool<PropertyType>::invalid_handle;
      }
    return *this;
  }



  template <int dim, int spacedim, typename PropertyType>
  Particle<dim,spacedim,PropertyType>::~Particle ()
  {
    if (properties != PropertyPool<PropertyType>::invalid_handle)
      property_pool->deallocate_properties_array(properties);
  }



  template <int dim, int spacedim, typename PropertyType>
  void
  Particle<dim,spacedim,PropertyType>::write_data (void *&data) const
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

    data = static_cast<void *> (pdata);

    // Write property data
    if (has_properties())
      {
        PropertyType *property_data = reinterpret_cast<PropertyType *> (data);

        const ArrayView<PropertyType> particle_properties = property_pool->get_properties(properties);
        for (unsigned int i = 0; i < particle_properties.size(); ++i,++property_data)
          *property_data = particle_properties[i];

        data = static_cast<void *> (pdata);
      }
  }



  template <int dim, int spacedim, typename PropertyType>
  std::size_t
  Particle<dim,spacedim,PropertyType>::size () const
  {
    std::size_t size = sizeof(types::particle_index)
                       + sizeof(location)
                       + sizeof(reference_location);

    if (has_properties())
      {
        const ArrayView<PropertyType> particle_properties = property_pool->get_properties(properties);
        size += sizeof(PropertyType) * particle_properties.size();
      }
    return size;
  }



  template <int dim, int spacedim, typename PropertyType>
  void
  Particle<dim,spacedim,PropertyType>::set_location (const Point<spacedim> &new_loc)
  {
    location = new_loc;
  }



  template <int dim, int spacedim, typename PropertyType>
  const Point<spacedim> &
  Particle<dim,spacedim,PropertyType>::get_location () const
  {
    return location;
  }



  template <int dim, int spacedim, typename PropertyType>
  void
  Particle<dim,spacedim,PropertyType>::set_reference_location (const Point<dim> &new_loc)
  {
    reference_location = new_loc;
  }



  template <int dim, int spacedim, typename PropertyType>
  const Point<dim> &
  Particle<dim,spacedim,PropertyType>::get_reference_location () const
  {
    return reference_location;
  }



  template <int dim, int spacedim, typename PropertyType>
  types::particle_index
  Particle<dim,spacedim,PropertyType>::get_id () const
  {
    return id;
  }



  template <int dim, int spacedim, typename PropertyType>
  void
  Particle<dim,spacedim,PropertyType>::set_property_pool (PropertyPool<PropertyType> &new_property_pool)
  {
    property_pool = &new_property_pool;
  }



  template <int dim, int spacedim, typename PropertyType>
  bool
  Particle<dim,spacedim,PropertyType>::has_properties () const
  {
    return (property_pool != NULL)
           && (properties != PropertyPool<PropertyType>::invalid_handle);
  }



  template <int dim, int spacedim, typename PropertyType>
  void
  Particle<dim,spacedim,PropertyType>::set_properties (const ArrayView<const PropertyType> &new_properties)
  {
    Assert(property_pool != NULL,
           ExcMessage("A particle was asked to set its properties without an "
                      "associated PropertyPool."))

    if (properties == PropertyPool<PropertyType>::invalid_handle)
      properties = property_pool->allocate_properties_array();

    const ArrayView<PropertyType> old_properties = property_pool->get_properties(properties);

    std::copy(new_properties.begin(),new_properties.end(),&old_properties[0]);
  }



  template <int dim, int spacedim, typename PropertyType>
  const ArrayView<const PropertyType>
  Particle<dim,spacedim,PropertyType>::get_properties () const
  {
    Assert(has_properties(),
           ExcMessage("A particle without additional properties was asked for "
                      "its properties."));

    return property_pool->get_properties(properties);
  }



  template <int dim, int spacedim, typename PropertyType>
  const ArrayView<PropertyType>
  Particle<dim,spacedim,PropertyType>::get_properties ()
  {
    Assert(has_properties(),
           ExcMessage("A particle without additional properties was asked for "
                      "its properties."));

    return property_pool->get_properties(properties);
  }
}

DEAL_II_NAMESPACE_CLOSE

DEAL_II_NAMESPACE_OPEN

#include "particle.inst"

DEAL_II_NAMESPACE_CLOSE
