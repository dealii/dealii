// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2021 by the deal.II authors
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

#ifndef dealii_particles_particle_accessor_h
#define dealii_particles_particle_accessor_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator_base.h>

#include <deal.II/particles/particle.h>

DEAL_II_NAMESPACE_OPEN

namespace Particles
{
  // Forward declarations
#ifndef DOXYGEN
  template <int, int>
  class ParticleIterator;
  template <int, int>
  class ParticleHandler;
#endif

  /**
   * Accessor class used by ParticleIterator to access particle data.
   */
  template <int dim, int spacedim = dim>
  class ParticleAccessor
  {
  public:
    /**
     * A type for the storage container for particles.
     */
    using particle_container =
      std::vector<std::vector<Particle<dim, spacedim>>>;

    /**
     * @copydoc Particle::write_particle_data_to_memory
     */
    void *
    write_particle_data_to_memory(void *data) const;


    /**
     * @copydoc Particle::read_particle_data_from_memory
     */
    const void *
    read_particle_data_from_memory(const void *data);

    /**
     * Set the location of this particle. Note that this does not check
     * whether this is a valid location in the simulation domain.
     *
     * @param [in] new_location The new location for this particle.
     *
     * @note In parallel programs, the ParticleHandler class stores particles
     *   on both the locally owned cells, as well as on ghost cells. The
     *   particles on the latter are *copies* of particles owned on other
     *   processors, and should therefore be treated in the same way as
     *   ghost entries in @ref GlossGhostedVector "vectors with ghost elements"
     *   or @ref GlossGhostCell "ghost cells": In both cases, one should
     *   treat the ghost elements or cells as `const` objects that shouldn't
     *   be modified even if the objects allow for calls that modify
     *   properties. Rather, properties should only be modified on processors
     *   that actually *own* the particle.
     */
    void
    set_location(const Point<spacedim> &new_location);

    /**
     * Get the location of this particle.
     *
     * @return The location of this particle.
     */
    const Point<spacedim> &
    get_location() const;

    /**
     * Set the reference location of this particle.
     *
     * @param [in] new_reference_location The new reference location for
     * this particle.
     *
     * @note In parallel programs, the ParticleHandler class stores particles
     *   on both the locally owned cells, as well as on ghost cells. The
     *   particles on the latter are *copies* of particles owned on other
     *   processors, and should therefore be treated in the same way as
     *   ghost entries in @ref GlossGhostedVector "vectors with ghost elements"
     *   or @ref GlossGhostCell "ghost cells": In both cases, one should
     *   treat the ghost elements or cells as `const` objects that shouldn't
     *   be modified even if the objects allow for calls that modify
     *   properties. Rather, properties should only be modified on processors
     *   that actually *own* the particle.
     */
    void
    set_reference_location(const Point<dim> &new_reference_location);

    /**
     * Return the reference location of this particle in its current cell.
     */
    const Point<dim> &
    get_reference_location() const;

    /**
     * Return the ID number of this particle.
     */
    types::particle_index
    get_id() const;

    /**
     * Tell the particle where to store its properties (even if it does not
     * own properties). Usually this is only done once per particle, but
     * since the particle generator does not know about the properties
     * we want to do it not at construction time. Another use for this
     * function is after particle transfer to a new process.
     */
    void
    set_property_pool(PropertyPool<dim, spacedim> &property_pool);

    /**
     * Return whether this particle has a valid property pool and a valid
     * handle to properties.
     */
    bool
    has_properties() const;

    /**
     * Set the properties of this particle.
     *
     * @param [in] new_properties A vector containing the
     * new properties for this particle.
     *
     * @note In parallel programs, the ParticleHandler class stores particles
     *   on both the locally owned cells, as well as on ghost cells. The
     *   particles on the latter are *copies* of particles owned on other
     *   processors, and should therefore be treated in the same way as
     *   ghost entries in @ref GlossGhostedVector "vectors with ghost elements"
     *   or @ref GlossGhostCell "ghost cells": In both cases, one should
     *   treat the ghost elements or cells as `const` objects that shouldn't
     *   be modified even if the objects allow for calls that modify
     *   properties. Rather, properties should only be modified on processors
     *   that actually *own* the particle.
     */
    void
    set_properties(const std::vector<double> &new_properties);

    /**
     * Set the properties of this particle.
     *
     * @param [in] new_properties An ArrayView pointing to memory locations
     * containing the new properties for this particle.
     *
     * @note In parallel programs, the ParticleHandler class stores particles
     *   on both the locally owned cells, as well as on ghost cells. The
     *   particles on the latter are *copies* of particles owned on other
     *   processors, and should therefore be treated in the same way as
     *   ghost entries in @ref GlossGhostedVector "vectors with ghost elements"
     *   or @ref GlossGhostCell "ghost cells": In both cases, one should
     *   treat the ghost elements or cells as `const` objects that shouldn't
     *   be modified even if the objects allow for calls that modify
     *   properties. Rather, properties should only be modified on processors
     *   that actually *own* the particle.
     */
    void
    set_properties(const ArrayView<const double> &new_properties);

    /**
     * Get write-access to properties of this particle.
     *
     * @return An ArrayView of the properties of this particle.
     */
    const ArrayView<double>
    get_properties();

    /**
     * Get read-access to properties of this particle.
     *
     * @return An ArrayView of the properties of this particle.
     */
    const ArrayView<const double>
    get_properties() const;

    /**
     * Return the size in bytes this particle occupies if all of its data is
     * serialized (i.e. the number of bytes that is written by the write_data
     * function of this class).
     */
    std::size_t
    serialized_size_in_bytes() const;

    /**
     * Get a cell iterator to the cell surrounding the current particle.
     * As particles are organized in the structure of a triangulation,
     * but the triangulation itself is not stored in the particle this
     * operation requires a reference to the triangulation.
     */
    const typename Triangulation<dim, spacedim>::cell_iterator &
    get_surrounding_cell() const;

    /**
     * @deprecated: Deprecated version of the function with the same
     * name above.
     */
    DEAL_II_DEPRECATED
    const typename Triangulation<dim, spacedim>::cell_iterator &
    get_surrounding_cell(
      const Triangulation<dim, spacedim> &triangulation) const;

    /**
     * Serialize the contents of this class using the [BOOST serialization
     * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
     */
    template <class Archive>
    void
    serialize(Archive &ar, const unsigned int version);

    /**
     * Advance the ParticleAccessor to the next particle.
     */
    void
    next();

    /**
     * Move the ParticleAccessor to the previous particle.
     */
    void
    prev();

    /**
     * Inequality operator.
     */
    bool
    operator!=(const ParticleAccessor<dim, spacedim> &other) const;

    /**
     * Equality operator.
     */
    bool
    operator==(const ParticleAccessor<dim, spacedim> &other) const;

    /**
     * Return the state of the accessor.
     */
    IteratorState::IteratorStates
    state() const;

  private:
    /**
     * Construct an invalid accessor. Such an object is not usable.
     */
    ParticleAccessor();

    /**
     * Construct an accessor from a reference to a container, an iterator to the
     * current cell, and the particle index within that cell.
     * This constructor is `private` so that it can only be accessed by
     * friend classes.
     */
    ParticleAccessor(
      const particle_container &particles,
      const typename Triangulation<dim, spacedim>::active_cell_iterator &cell,
      const unsigned int particle_index_within_cell);

    /**
     * Returns a reference to the current Particle. Because the
     * internal structure may change this is not intended for public use and
     * only a convenience function for internal purposes.
     */
    const Particle<dim, spacedim> &
    get_particle() const;

    /**
     * Non-@p const version of the function with the same name above.
     */
    Particle<dim, spacedim> &
    get_particle();

    /**
     * A pointer to the container that stores the particles. Obviously,
     * this accessor is invalidated if the container changes.
     */
    particle_container *particles;

    /**
     * Cell iterator to the cell of the current particle.
     */
    typename Triangulation<dim, spacedim>::active_cell_iterator cell;

    /**
     * Index to the cell this particle is stored in at the moment. This could be
     * read from the member variable cell, but is used in many performance
     * critical functions and is therefore stored and updated separately.
     */
    unsigned int active_cell_index;

    /**
     * Local index of the particle within its current cell.
     */
    unsigned int particle_index_within_cell;

    // Make ParticleIterator a friend to allow it constructing
    // ParticleAccessors.
    template <int, int>
    friend class ParticleIterator;
    template <int, int>
    friend class ParticleHandler;
  };



  template <int dim, int spacedim>
  template <class Archive>
  void
  ParticleAccessor<dim, spacedim>::serialize(Archive &          ar,
                                             const unsigned int version)
  {
    return get_particle().serialize(ar, version);
  }


  // ------------------------- inline functions ------------------------------

  template <int dim, int spacedim>
  inline ParticleAccessor<dim, spacedim>::ParticleAccessor()
    : particles(nullptr)
    , cell()
    , active_cell_index(numbers::invalid_unsigned_int)
    , particle_index_within_cell(numbers::invalid_unsigned_int)
  {}



  template <int dim, int spacedim>
  inline ParticleAccessor<dim, spacedim>::ParticleAccessor(
    const particle_container &particles,
    const typename Triangulation<dim, spacedim>::active_cell_iterator &cell,
    const unsigned int particle_index_within_cell)
    : particles(const_cast<particle_container *>(&particles))
    , cell(cell)
    , particle_index_within_cell(particle_index_within_cell)
  {
    if (cell.state() == IteratorState::valid)
      active_cell_index = cell->active_cell_index();
    else if (cell.state() == IteratorState::past_the_end)
      active_cell_index = particles.size();
    else
      active_cell_index = numbers::invalid_unsigned_int;
  }



  template <int dim, int spacedim>
  inline const void *
  ParticleAccessor<dim, spacedim>::read_particle_data_from_memory(
    const void *data)
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    return get_particle().read_particle_data_from_memory(data);
  }



  template <int dim, int spacedim>
  inline void *
  ParticleAccessor<dim, spacedim>::write_particle_data_to_memory(
    void *data) const
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    return get_particle().write_particle_data_to_memory(data);
  }



  template <int dim, int spacedim>
  inline void
  ParticleAccessor<dim, spacedim>::set_location(const Point<spacedim> &new_loc)
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    get_particle().set_location(new_loc);
  }



  template <int dim, int spacedim>
  inline const Point<spacedim> &
  ParticleAccessor<dim, spacedim>::get_location() const
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    return get_particle().get_location();
  }



  template <int dim, int spacedim>
  inline void
  ParticleAccessor<dim, spacedim>::set_reference_location(
    const Point<dim> &new_loc)
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    get_particle().set_reference_location(new_loc);
  }



  template <int dim, int spacedim>
  inline const Point<dim> &
  ParticleAccessor<dim, spacedim>::get_reference_location() const
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    return get_particle().get_reference_location();
  }



  template <int dim, int spacedim>
  inline types::particle_index
  ParticleAccessor<dim, spacedim>::get_id() const
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    return get_particle().get_id();
  }



  template <int dim, int spacedim>
  inline void
  ParticleAccessor<dim, spacedim>::set_property_pool(
    PropertyPool<dim, spacedim> &new_property_pool)
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    get_particle().set_property_pool(new_property_pool);
  }



  template <int dim, int spacedim>
  inline bool
  ParticleAccessor<dim, spacedim>::has_properties() const
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    return get_particle().has_properties();
  }



  template <int dim, int spacedim>
  inline void
  ParticleAccessor<dim, spacedim>::set_properties(
    const std::vector<double> &new_properties)
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    get_particle().set_properties(new_properties);
  }



  template <int dim, int spacedim>
  inline void
  ParticleAccessor<dim, spacedim>::set_properties(
    const ArrayView<const double> &new_properties)
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    get_particle().set_properties(new_properties);
  }



  template <int dim, int spacedim>
  inline const ArrayView<const double>
  ParticleAccessor<dim, spacedim>::get_properties() const
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    return get_particle().get_properties();
  }



  template <int dim, int spacedim>
  inline const typename Triangulation<dim, spacedim>::cell_iterator &
  ParticleAccessor<dim, spacedim>::get_surrounding_cell() const
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    return cell;
  }



  template <int dim, int spacedim>
  inline const typename Triangulation<dim, spacedim>::cell_iterator &
  ParticleAccessor<dim, spacedim>::get_surrounding_cell(
    const Triangulation<dim, spacedim> & /*triangulation*/) const
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    return cell;
  }



  template <int dim, int spacedim>
  inline const ArrayView<double>
  ParticleAccessor<dim, spacedim>::get_properties()
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    return get_particle().get_properties();
  }



  template <int dim, int spacedim>
  inline std::size_t
  ParticleAccessor<dim, spacedim>::serialized_size_in_bytes() const
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    return get_particle().serialized_size_in_bytes();
  }



  template <int dim, int spacedim>
  inline void
  ParticleAccessor<dim, spacedim>::next()
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    ++particle_index_within_cell;

    if (particle_index_within_cell > (*particles)[active_cell_index].size() - 1)
      {
        particle_index_within_cell = 0;

        do
          {
            ++cell;
            ++active_cell_index;
          }
        while (cell.state() == IteratorState::valid &&
               (*particles)[active_cell_index].size() == 0);
      }
  }



  template <int dim, int spacedim>
  inline void
  ParticleAccessor<dim, spacedim>::prev()
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    if (particle_index_within_cell > 0)
      --particle_index_within_cell;
    else
      {
        do
          {
            --cell;
            --active_cell_index;

            if (cell.state() != IteratorState::valid)
              {
                active_cell_index = particles->size();
                break;
              }

            if ((*particles)[active_cell_index].size() > 0)
              {
                particle_index_within_cell =
                  (*particles)[active_cell_index].size() - 1;
                break;
              }
          }
        while ((*particles)[active_cell_index].size() == 0);
      }
  }



  template <int dim, int spacedim>
  inline bool
  ParticleAccessor<dim, spacedim>::
  operator!=(const ParticleAccessor<dim, spacedim> &other) const
  {
    return !(*this == other);
  }



  template <int dim, int spacedim>
  inline bool
  ParticleAccessor<dim, spacedim>::
  operator==(const ParticleAccessor<dim, spacedim> &other) const
  {
    return (particles == other.particles) && (cell == other.cell) &&
           (particle_index_within_cell == other.particle_index_within_cell);
  }



  template <int dim, int spacedim>
  inline IteratorState::IteratorStates
  ParticleAccessor<dim, spacedim>::state() const
  {
    if (particles != nullptr && cell.state() == IteratorState::valid &&
        particle_index_within_cell < (*particles)[active_cell_index].size())
      return IteratorState::valid;
    else if (particles != nullptr &&
             cell.state() == IteratorState::past_the_end &&
             particle_index_within_cell == 0)
      return IteratorState::past_the_end;
    else
      return IteratorState::invalid;

    return IteratorState::invalid;
  }



  template <int dim, int spacedim>
  inline Particle<dim, spacedim> &
  ParticleAccessor<dim, spacedim>::get_particle()
  {
    return (*particles)[active_cell_index][particle_index_within_cell];
  }



  template <int dim, int spacedim>
  inline const Particle<dim, spacedim> &
  ParticleAccessor<dim, spacedim>::get_particle() const
  {
    return (*particles)[active_cell_index][particle_index_within_cell];
  }

} // namespace Particles

DEAL_II_NAMESPACE_CLOSE

namespace boost
{
  namespace geometry
  {
    namespace index
    {
      // Forward declaration of bgi::indexable
      template <class T>
      struct indexable;

      /**
       * Make sure we can construct an RTree from Particles::ParticleAccessor
       * objects.
       */
      template <int dim, int spacedim>
      struct indexable<dealii::Particles::ParticleAccessor<dim, spacedim>>
      {
        /**
         * boost::rtree expects a const reference to an indexable object. For
         * a Particles::Particle object, this is its reference location.
         */
        using result_type = const dealii::Point<spacedim> &;

        result_type
        operator()(const dealii::Particles::ParticleAccessor<dim, spacedim>
                     &accessor) const
        {
          return accessor.get_location();
        }
      };
    } // namespace index
  }   // namespace geometry
} // namespace boost

#endif
