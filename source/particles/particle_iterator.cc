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

#include <deal.II/particles/particle_iterator.h>

DEAL_II_NAMESPACE_OPEN

namespace Particles
{
  template <int dim, int spacedim>
  ParticleIterator<dim,spacedim>::ParticleIterator (const std::multimap<types::LevelInd, Particle<dim,spacedim> > &map,
                                                    const typename std::multimap<types::LevelInd, Particle<dim,spacedim> >::iterator &particle)
    :
    accessor (map, particle)
  {}



  template <int dim, int spacedim>
  ParticleAccessor<dim,spacedim> &
  ParticleIterator<dim,spacedim>::operator *()
  {
    return accessor;
  }



  template <int dim, int spacedim>
  ParticleAccessor<dim,spacedim> *
  ParticleIterator<dim,spacedim>::operator ->()
  {
    return &(this->operator* ());
  }



  template <int dim, int spacedim>
  const ParticleAccessor<dim,spacedim> &
  ParticleIterator<dim,spacedim>::operator *() const
  {
    return accessor;
  }



  template <int dim, int spacedim>
  const ParticleAccessor<dim,spacedim> *
  ParticleIterator<dim,spacedim>::operator ->() const
  {
    return &(this->operator* ());
  }


  template <int dim, int spacedim>
  ParticleIterator<dim,spacedim> &
  ParticleIterator<dim,spacedim>::operator =(const ParticleIterator &other)
  {
    accessor = other.accessor;
    return *this;
  }



  template <int dim, int spacedim>
  bool
  ParticleIterator<dim,spacedim>::operator != (const ParticleIterator<dim,spacedim> &other) const
  {
    return accessor != other.accessor;
  }



  template <int dim, int spacedim>
  bool
  ParticleIterator<dim,spacedim>::operator == (const ParticleIterator<dim,spacedim> &other) const
  {
    return accessor == other.accessor;
  }



  template <int dim, int spacedim>
  ParticleIterator<dim,spacedim> &
  ParticleIterator<dim,spacedim>::operator++()
  {
    accessor.next();
    return *this;
  }



  template <int dim, int spacedim>
  ParticleIterator<dim,spacedim>
  ParticleIterator<dim,spacedim>::operator++(int)
  {
    ParticleIterator tmp(*this);
    operator++ ();

    return tmp;
  }



  template <int dim, int spacedim>
  ParticleIterator<dim,spacedim> &
  ParticleIterator<dim,spacedim>::operator--()
  {
    accessor.prev();
    return *this;
  }



  template <int dim, int spacedim>
  ParticleIterator<dim,spacedim>
  ParticleIterator<dim,spacedim>::operator--(int)
  {
    ParticleIterator tmp(*this);
    operator-- ();

    return tmp;
  }
}


DEAL_II_NAMESPACE_CLOSE

DEAL_II_NAMESPACE_OPEN

#include "particle_iterator.inst"

DEAL_II_NAMESPACE_CLOSE
