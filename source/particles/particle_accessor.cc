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

#include <deal.II/particles/particle_accessor.h>

DEAL_II_NAMESPACE_OPEN

namespace Particles
{
  template <int dim, int spacedim>
  ParticleAccessor<dim, spacedim>::ParticleAccessor()
    : map(nullptr)
    , particle()
  {}



  template <int dim, int spacedim>
  ParticleAccessor<dim, spacedim>::ParticleAccessor(
    const std::multimap<internal::LevelInd, Particle<dim, spacedim>> &map,
    const typename std::multimap<internal::LevelInd,
                                 Particle<dim, spacedim>>::iterator & particle)
    : map(const_cast<
          std::multimap<internal::LevelInd, Particle<dim, spacedim>> *>(&map))
    , particle(particle)
  {}



  template <int dim, int spacedim>
  void
  ParticleAccessor<dim, spacedim>::write_data(void *&data) const
  {
    Assert(particle != map->end(), ExcInternalError());

    particle->second.write_data(data);
  }



  template <int dim, int spacedim>
  void
  ParticleAccessor<dim, spacedim>::set_location(const Point<spacedim> &new_loc)
  {
    Assert(particle != map->end(), ExcInternalError());

    particle->second.set_location(new_loc);
  }



  template <int dim, int spacedim>
  const Point<spacedim> &
  ParticleAccessor<dim, spacedim>::get_location() const
  {
    Assert(particle != map->end(), ExcInternalError());

    return particle->second.get_location();
  }



  template <int dim, int spacedim>
  void
  ParticleAccessor<dim, spacedim>::set_reference_location(
    const Point<dim> &new_loc)
  {
    Assert(particle != map->end(), ExcInternalError());

    particle->second.set_reference_location(new_loc);
  }



  template <int dim, int spacedim>
  const Point<dim> &
  ParticleAccessor<dim, spacedim>::get_reference_location() const
  {
    Assert(particle != map->end(), ExcInternalError());

    return particle->second.get_reference_location();
  }



  template <int dim, int spacedim>
  types::particle_index
  ParticleAccessor<dim, spacedim>::get_id() const
  {
    Assert(particle != map->end(), ExcInternalError());

    return particle->second.get_id();
  }



  template <int dim, int spacedim>
  void
  ParticleAccessor<dim, spacedim>::set_property_pool(
    PropertyPool &new_property_pool)
  {
    Assert(particle != map->end(), ExcInternalError());

    particle->second.set_property_pool(new_property_pool);
  }



  template <int dim, int spacedim>
  bool
  ParticleAccessor<dim, spacedim>::has_properties() const
  {
    Assert(particle != map->end(), ExcInternalError());

    return particle->second.has_properties();
  }



  template <int dim, int spacedim>
  void
  ParticleAccessor<dim, spacedim>::set_properties(
    const std::vector<double> &new_properties)
  {
    Assert(particle != map->end(), ExcInternalError());

    particle->second.set_properties(new_properties);
  }



  template <int dim, int spacedim>
  void
  ParticleAccessor<dim, spacedim>::set_properties(
    const ArrayView<const double> &new_properties)
  {
    Assert(particle != map->end(), ExcInternalError());

    particle->second.set_properties(new_properties);
  }



  template <int dim, int spacedim>
  const ArrayView<const double>
  ParticleAccessor<dim, spacedim>::get_properties() const
  {
    Assert(particle != map->end(), ExcInternalError());

    return particle->second.get_properties();
  }



  template <int dim, int spacedim>
  typename Triangulation<dim, spacedim>::cell_iterator
  ParticleAccessor<dim, spacedim>::get_surrounding_cell(
    const Triangulation<dim, spacedim> &triangulation) const
  {
    Assert(particle != map->end(), ExcInternalError());

    const typename Triangulation<dim, spacedim>::cell_iterator cell(
      &triangulation, particle->first.first, particle->first.second);
    return cell;
  }



  template <int dim, int spacedim>
  const ArrayView<double>
  ParticleAccessor<dim, spacedim>::get_properties()
  {
    Assert(particle != map->end(), ExcInternalError());

    return particle->second.get_properties();
  }



  template <int dim, int spacedim>
  std::size_t
  ParticleAccessor<dim, spacedim>::serialized_size_in_bytes() const
  {
    Assert(particle != map->end(), ExcInternalError());

    return particle->second.serialized_size_in_bytes();
  }



  template <int dim, int spacedim>
  void
  ParticleAccessor<dim, spacedim>::next()
  {
    Assert(particle != map->end(), ExcInternalError());
    ++particle;
  }



  template <int dim, int spacedim>
  void
  ParticleAccessor<dim, spacedim>::prev()
  {
    Assert(particle != map->begin(), ExcInternalError());
    --particle;
  }



  template <int dim, int spacedim>
  bool
  ParticleAccessor<dim, spacedim>::
  operator!=(const ParticleAccessor<dim, spacedim> &other) const
  {
    return (map != other.map) || (particle != other.particle);
  }



  template <int dim, int spacedim>
  bool
  ParticleAccessor<dim, spacedim>::
  operator==(const ParticleAccessor<dim, spacedim> &other) const
  {
    return (map == other.map) && (particle == other.particle);
  }
} // namespace Particles

DEAL_II_NAMESPACE_CLOSE

DEAL_II_NAMESPACE_OPEN

#include "particle_accessor.inst"

DEAL_II_NAMESPACE_CLOSE
