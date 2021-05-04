// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2019 by the deal.II authors
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
    typename Triangulation<dim, spacedim>::cell_iterator
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

  private:
    /**
     * Construct an invalid accessor. Such an object is not usable.
     */
    ParticleAccessor();

    /**
     * Construct an accessor from a reference to a map and an iterator to the
     * map. This constructor is `private` so that it can only be accessed by
     * friend classes.
     */
    ParticleAccessor(
      const std::multimap<internal::LevelInd, Particle<dim, spacedim>> &map,
      const typename std::multimap<internal::LevelInd,
                                   Particle<dim, spacedim>>::iterator
        &particle);

  private:
    /**
     * A pointer to the container that stores the particles. Obviously,
     * this accessor is invalidated if the container changes.
     */
    std::multimap<internal::LevelInd, Particle<dim, spacedim>> *map;

    /**
     * An iterator into the container of particles. Obviously,
     * this accessor is invalidated if the container changes.
     */
    typename std::multimap<internal::LevelInd,
                           Particle<dim, spacedim>>::iterator particle;

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
    return particle->second.serialize(ar, version);
  }


  // ------------------------- inline functions ------------------------------

  template <int dim, int spacedim>
  inline ParticleAccessor<dim, spacedim>::ParticleAccessor()
    : map(nullptr)
    , particle()
  {}



  template <int dim, int spacedim>
  inline ParticleAccessor<dim, spacedim>::ParticleAccessor(
    const std::multimap<internal::LevelInd, Particle<dim, spacedim>> &map,
    const typename std::multimap<internal::LevelInd,
                                 Particle<dim, spacedim>>::iterator & particle)
    : map(const_cast<
          std::multimap<internal::LevelInd, Particle<dim, spacedim>> *>(&map))
    , particle(particle)
  {}



  template <int dim, int spacedim>
  inline const void *
  ParticleAccessor<dim, spacedim>::read_particle_data_from_memory(
    const void *data)
  {
    Assert(particle != map->end(), ExcInternalError());

    return particle->second.read_particle_data_from_memory(data);
  }



  template <int dim, int spacedim>
  inline void *
  ParticleAccessor<dim, spacedim>::write_particle_data_to_memory(
    void *data) const
  {
    Assert(particle != map->end(), ExcInternalError());

    return particle->second.write_particle_data_to_memory(data);
  }



  template <int dim, int spacedim>
  inline void
  ParticleAccessor<dim, spacedim>::set_location(const Point<spacedim> &new_loc)
  {
    Assert(particle != map->end(), ExcInternalError());

    particle->second.set_location(new_loc);
  }



  template <int dim, int spacedim>
  inline const Point<spacedim> &
  ParticleAccessor<dim, spacedim>::get_location() const
  {
    Assert(particle != map->end(), ExcInternalError());

    return particle->second.get_location();
  }



  template <int dim, int spacedim>
  inline void
  ParticleAccessor<dim, spacedim>::set_reference_location(
    const Point<dim> &new_loc)
  {
    Assert(particle != map->end(), ExcInternalError());

    particle->second.set_reference_location(new_loc);
  }



  template <int dim, int spacedim>
  inline const Point<dim> &
  ParticleAccessor<dim, spacedim>::get_reference_location() const
  {
    Assert(particle != map->end(), ExcInternalError());

    return particle->second.get_reference_location();
  }



  template <int dim, int spacedim>
  inline types::particle_index
  ParticleAccessor<dim, spacedim>::get_id() const
  {
    Assert(particle != map->end(), ExcInternalError());

    return particle->second.get_id();
  }



  template <int dim, int spacedim>
  inline void
  ParticleAccessor<dim, spacedim>::set_property_pool(
    PropertyPool<dim, spacedim> &new_property_pool)
  {
    Assert(particle != map->end(), ExcInternalError());

    particle->second.set_property_pool(new_property_pool);
  }



  template <int dim, int spacedim>
  inline bool
  ParticleAccessor<dim, spacedim>::has_properties() const
  {
    Assert(particle != map->end(), ExcInternalError());

    return particle->second.has_properties();
  }



  template <int dim, int spacedim>
  inline void
  ParticleAccessor<dim, spacedim>::set_properties(
    const std::vector<double> &new_properties)
  {
    Assert(particle != map->end(), ExcInternalError());

    particle->second.set_properties(new_properties);
  }



  template <int dim, int spacedim>
  inline void
  ParticleAccessor<dim, spacedim>::set_properties(
    const ArrayView<const double> &new_properties)
  {
    Assert(particle != map->end(), ExcInternalError());

    particle->second.set_properties(new_properties);
  }



  template <int dim, int spacedim>
  inline const ArrayView<const double>
  ParticleAccessor<dim, spacedim>::get_properties() const
  {
    Assert(particle != map->end(), ExcInternalError());

    return particle->second.get_properties();
  }



  template <int dim, int spacedim>
  inline typename Triangulation<dim, spacedim>::cell_iterator
  ParticleAccessor<dim, spacedim>::get_surrounding_cell(
    const Triangulation<dim, spacedim> &triangulation) const
  {
    Assert(particle != map->end(), ExcInternalError());

    const typename Triangulation<dim, spacedim>::cell_iterator cell(
      &triangulation, particle->first.first, particle->first.second);
    return cell;
  }



  template <int dim, int spacedim>
  inline const ArrayView<double>
  ParticleAccessor<dim, spacedim>::get_properties()
  {
    Assert(particle != map->end(), ExcInternalError());

    return particle->second.get_properties();
  }



  template <int dim, int spacedim>
  inline std::size_t
  ParticleAccessor<dim, spacedim>::serialized_size_in_bytes() const
  {
    Assert(particle != map->end(), ExcInternalError());

    return particle->second.serialized_size_in_bytes();
  }



  template <int dim, int spacedim>
  inline void
  ParticleAccessor<dim, spacedim>::next()
  {
    Assert(particle != map->end(), ExcInternalError());
    ++particle;
  }



  template <int dim, int spacedim>
  inline void
  ParticleAccessor<dim, spacedim>::prev()
  {
    Assert(particle != map->begin(), ExcInternalError());
    --particle;
  }



  template <int dim, int spacedim>
  inline bool
  ParticleAccessor<dim, spacedim>::
  operator!=(const ParticleAccessor<dim, spacedim> &other) const
  {
    return (map != other.map) || (particle != other.particle);
  }



  template <int dim, int spacedim>
  inline bool
  ParticleAccessor<dim, spacedim>::
  operator==(const ParticleAccessor<dim, spacedim> &other) const
  {
    return (map == other.map) && (particle == other.particle);
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
