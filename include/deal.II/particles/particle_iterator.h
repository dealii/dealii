// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_particles_particle_iterator_h
#define dealii_particles_particle_iterator_h

#include <deal.II/base/config.h>

#include <deal.II/particles/particle_accessor.h>

DEAL_II_NAMESPACE_OPEN

namespace Particles
{
  // Forward declaration
#ifndef DOXYGEN
  template <int, int>
  class ParticleHandler;
#endif

  /**
   * A class that is used to iterate over particles. Together with the
   * ParticleAccessor class this is used to hide the internal implementation
   * of the particle class and the particle container.
   */
  template <int dim, int spacedim = dim>
  class ParticleIterator
  {
  public:
    /**
     * A type for the storage container for particles.
     */
    using particle_container =
      typename ParticleAccessor<dim, spacedim>::particle_container;

    /**
     * Empty constructor. Such an object is not usable!
     */
    ParticleIterator() = default;

    /**
     * Constructor of the iterator. Takes a reference to the particle
     * container, an iterator to the cell, and the particle index within that
     * cell.
     */
    ParticleIterator(
      const typename particle_container::iterator particles_in_cell,
      const PropertyPool<dim, spacedim>          &property_pool,
      const unsigned int                          particle_index_within_cell);

    /**
     * Dereferencing operator, returns a reference to an accessor. Usage is thus
     * like <tt>(*i).get_id ();</tt>
     */
    const ParticleAccessor<dim, spacedim> &
    operator*() const;

    /**
     * Dereferencing operator, non-@p const version.
     */
    ParticleAccessor<dim, spacedim> &
    operator*();

    /**
     * Dereferencing operator, returns a pointer of the particle pointed to.
     * Usage is thus like <tt>i->get_id ();</tt>
     *
     * There is a @p const and a non-@p const version.
     */
    const ParticleAccessor<dim, spacedim> *
    operator->() const;

    /**
     * Dereferencing operator, non-@p const version.
     */
    ParticleAccessor<dim, spacedim> *
    operator->();

    /**
     * Compare for equality.
     */
    bool
    operator==(const ParticleIterator<dim, spacedim> &) const;

    /**
     * Compare for inequality.
     */
    bool
    operator!=(const ParticleIterator<dim, spacedim> &) const;

    /**
     * Prefix <tt>++</tt> operator: <tt>++iterator</tt>. This operator advances
     * the iterator to the next element and returns a reference to
     * <tt>*this</tt>.
     */
    ParticleIterator &
    operator++();

    /**
     * Postfix <tt>++</tt> operator: <tt>iterator++</tt>. This operator advances
     * the iterator to the next element, but returns an iterator to the element
     * previously pointed to.
     */
    ParticleIterator
    operator++(int);

    /**
     * Prefix <tt>\--</tt> operator: <tt>\--iterator</tt>. This operator moves
     * the iterator to the previous element and returns a reference to
     * <tt>*this</tt>.
     */
    ParticleIterator &
    operator--();

    /**
     * Postfix <tt>\--</tt> operator: <tt>iterator\--</tt>. This operator moves
     * the iterator to the previous element, but returns an iterator to the
     * element previously pointed to.
     */
    ParticleIterator
    operator--(int);

    /**
     * Return the state of the iterator.
     */
    IteratorState::IteratorStates
    state() const;

    /**
     * Mark the class as bidirectional iterator and declare some alias which
     * are standard for iterators and are used by algorithms to enquire about
     * the specifics of the iterators they work on.
     */
    using iterator_category = std::bidirectional_iterator_tag;
    using value_type        = ParticleAccessor<dim, spacedim>;
    using difference_type   = std::ptrdiff_t;
    using pointer           = ParticleAccessor<dim, spacedim> *;
    using reference         = ParticleAccessor<dim, spacedim> &;

  private:
    /**
     * The accessor to the actual particle.
     */
    ParticleAccessor<dim, spacedim> accessor;
  };



  // ------------------------------ inline functions -------------------------

  template <int dim, int spacedim>
  inline ParticleIterator<dim, spacedim>::ParticleIterator(
    const typename particle_container::iterator particles_in_cell,
    const PropertyPool<dim, spacedim>          &property_pool,
    const unsigned int                          particle_index_within_cell)
    : accessor(particles_in_cell, property_pool, particle_index_within_cell)
  {}



  template <int dim, int spacedim>
  inline ParticleAccessor<dim, spacedim> &
  ParticleIterator<dim, spacedim>::operator*()
  {
    return accessor;
  }



  template <int dim, int spacedim>
  inline ParticleAccessor<dim, spacedim> *
  ParticleIterator<dim, spacedim>::operator->()
  {
    return &(this->operator*());
  }



  template <int dim, int spacedim>
  inline const ParticleAccessor<dim, spacedim> &
  ParticleIterator<dim, spacedim>::operator*() const
  {
    return accessor;
  }



  template <int dim, int spacedim>
  inline const ParticleAccessor<dim, spacedim> *
  ParticleIterator<dim, spacedim>::operator->() const
  {
    return &(this->operator*());
  }



  template <int dim, int spacedim>
  inline bool
  ParticleIterator<dim, spacedim>::operator!=(
    const ParticleIterator<dim, spacedim> &other) const
  {
    return accessor != other.accessor;
  }



  template <int dim, int spacedim>
  inline bool
  ParticleIterator<dim, spacedim>::operator==(
    const ParticleIterator<dim, spacedim> &other) const
  {
    return accessor == other.accessor;
  }



  template <int dim, int spacedim>
  inline ParticleIterator<dim, spacedim> &
  ParticleIterator<dim, spacedim>::operator++()
  {
    accessor.next();
    return *this;
  }



  template <int dim, int spacedim>
  inline ParticleIterator<dim, spacedim>
  ParticleIterator<dim, spacedim>::operator++(int)
  {
    ParticleIterator tmp(*this);
    operator++();

    return tmp;
  }



  template <int dim, int spacedim>
  inline ParticleIterator<dim, spacedim> &
  ParticleIterator<dim, spacedim>::operator--()
  {
    accessor.prev();
    return *this;
  }



  template <int dim, int spacedim>
  inline ParticleIterator<dim, spacedim>
  ParticleIterator<dim, spacedim>::operator--(int)
  {
    ParticleIterator tmp(*this);
    operator--();

    return tmp;
  }



  template <int dim, int spacedim>
  inline IteratorState::IteratorStates
  ParticleIterator<dim, spacedim>::state() const
  {
    return accessor.state();
  }

} // namespace Particles

DEAL_II_NAMESPACE_CLOSE

#endif
