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

#ifndef dealii_particles_particle_accessor_h
#define dealii_particles_particle_accessor_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator_base.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/property_pool.h>

#include <boost/geometry/index/indexable.hpp>
#include <boost/serialization/array_wrapper.hpp>

#include <list>


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
     * Data structure to describe the particles in a given cell. This is used
     * inside an std::list in `particle_container`. The storage of this field
     * is typically handled by ParticleHandler, using an std::list of this
     * structure.
     *
     * There are four main reasons for the present design:
     * <ul>
     * <li>The list represents a compact storage for all particles on a cell,
     * handling e.g. the case where only a few cells hold particles without
     * storage overhead.</li>
     * <li>The particles on a cell are identified by their handle object,
     * i.e., a single integer, again ensuring only 4 byte storage cost for the
     * handling of the particles for the case of sufficiently many particles
     * per cell.</li>
     * <li>Combined with a cache object that holds an iterator to the
     * `std::list` of particles (8 bytes per cell), or the `std::list::end` in
     * case no particles are present on a cell, there is a fast access from a
     * ParticleAccessor -> surrounding cell (by access to
     * ParticlesInCell::cell_iterator), as well as fast access from cell ->
     * all particles (through the cache). It also allows for a fast iteration
     * through all particles, by incrementing either the index of particles
     * within a cell, or, if at the end of the cell, to the next element in
     * the outer list. The cache is simple to keep consistent because the
     * iterators into `std::list` remain valid upon insertion or deletion of
     * entries in the list, as specified by `std::list`'s API.
     * <li>In order to detect the start and end of all cells with particles,
     * the std::list contains some dummy entry with a `cell_iterator` past the
     * end of valid cells, which is used as a criterion to terminate the loops
     * of ParticleAccessor, again minimizing the computational cost of handling
     * the loop over particles.</li>
     * </ul>
     */
    struct ParticlesInCell
    {
      /**
       * Default constructor.
       */
      ParticlesInCell() = default;

      /**
       * Construct from a vector of particles and a cell iterator.
       */
      ParticlesInCell(
        const std::vector<typename PropertyPool<dim, spacedim>::Handle>
          &particles,
        const typename Triangulation<dim, spacedim>::active_cell_iterator &cell)
        : particles(particles)
        , cell(cell)
      {}

      /**
       * A vector of particles on a cell.
       */
      std::vector<typename PropertyPool<dim, spacedim>::Handle> particles;

      /**
       * The underlying cell.
       */
      typename Triangulation<dim, spacedim>::active_cell_iterator cell;
    };

    /**
     * A type for the storage container for particles.
     */
    using particle_container = std::list<ParticlesInCell>;

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
     *   ghost entries in
     *   @ref GlossGhostedVector "vectors with ghost elements"
     *   or
     *   @ref GlossGhostCell "ghost cells":
     *   In both cases, one should
     *   treat the ghost elements or cells as `const` objects that shouldn't
     *   be modified even if the objects allow for calls that modify
     *   properties. Rather, properties should only be modified on processors
     *   that actually *own* the particle.
     */
    void
    set_location(const Point<spacedim> &new_location);

    /**
     * Get read-access to the location of this particle.
     *
     * @return The location of this particle.
     */
    const Point<spacedim> &
    get_location() const;

    /**
     * Get read- and write-access to the location of this particle.
     * Note that changing the location does not check
     * whether this is a valid location in the simulation domain.
     *
     * @note In parallel programs, the ParticleHandler class stores particles
     *   on both the locally owned cells, as well as on ghost cells. The
     *   particles on the latter are *copies* of particles owned on other
     *   processors, and should therefore be treated in the same way as
     *   ghost entries in
     *   @ref GlossGhostedVector "vectors with ghost elements"
     *   or
     *   @ref GlossGhostCell "ghost cells":
     *   In both cases, one should
     *   treat the ghost elements or cells as `const` objects that shouldn't
     *   be modified even if the objects allow for calls that modify
     *   properties. Rather, properties should only be modified on processors
     *   that actually *own* the particle.
     *
     * @return The location of this particle.
     */
    Point<spacedim> &
    get_location();

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
     *   ghost entries in
     *   @ref GlossGhostedVector "vectors with ghost elements"
     *   or
     *   @ref GlossGhostCell "ghost cells":
     *   In both cases, one should
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
     * Set the ID number of this particle.
     */
    void
    set_id(const types::particle_index &new_id);

    /**
     * Return the ID number of this particle.
     */
    types::particle_index
    get_id() const;

    /**
     * Return a particle ID number local to each MPI process. This number
     * enables the direct array access (similar to
     * LinearAlgebra::distributed::Vector::local_element()) to quantities used
     * for local computations. Use
     * ParticleHandler::get_max_local_particle_index() to query suitable array
     * sizes.
     *
     * @note The number returned by this function is not stable between calls
     * of ParticleHandler::sort_particles_into_subdomains_and_cells() and
     * therefore it cannot be used to track quantities across those
     * calls. Furthermore, the numbers returned by this function are typically
     * not refreshed in case individual particles are removed from the
     * underlying ParticleHandler object, which means that the returned
     * indices of all particles on an MPI process typically do not form a
     * contiguous interval of numbers. Use particle properties for storing
     * persistent information and checkpointing the actual data.
     */
    types::particle_index
    get_local_index() const;

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
     *   ghost entries in
     *   @ref GlossGhostedVector "vectors with ghost elements"
     *   or
     *   @ref GlossGhostCell "ghost cells":
     *    In both cases, one should
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
     *   ghost entries in
     *   @ref GlossGhostedVector "vectors with ghost elements"
     *   or
     *   @ref GlossGhostCell "ghost cells":
     *   In both cases, one should
     *   treat the ghost elements or cells as `const` objects that shouldn't
     *   be modified even if the objects allow for calls that modify
     *   properties. Rather, properties should only be modified on processors
     *   that actually *own* the particle.
     */
    void
    set_properties(const ArrayView<const double> &new_properties);

    /**
     * Set the properties of this particle, assuming that the properties
     * stored on this particle correspond to a rank-1 Tensor object. In
     * particular, this means that the number of properties stored on
     * the particle must equal `dim`.
     *
     * @param [in] new_properties A Tensor containing the new properties for
     *   this particle.
     *
     * @note In parallel programs, the ParticleHandler class stores particles
     *   on both the locally owned cells, as well as on ghost cells. The
     *   particles on the latter are *copies* of particles owned on other
     *   processors, and should therefore be treated in the same way as
     *   ghost entries in
     *   @ref GlossGhostedVector "vectors with ghost elements"
     *   or
     *   @ref GlossGhostCell "ghost cells":
     *   In both cases, one should
     *   treat the ghost elements or cells as `const` objects that shouldn't
     *   be modified even if the objects allow for calls that modify
     *   properties. Rather, properties should only be modified on processors
     *   that actually *own* the particle.
     */
    void
    set_properties(const Tensor<1, dim> &new_properties);

    /**
     * Get write-access to properties of this particle.
     *
     * @return An ArrayView of the properties of this particle.
     */
    ArrayView<double>
    get_properties()
    {
      // The implementation is up here inside the class declaration because
      // NVCC (at least in 12.5 and 12.6) otherwise produce a compile error:
      //
      // error: no declaration matches ‘dealii::ArrayView<__remove_cv(const
      // double)> dealii::Particles::ParticleAccessor<dim,
      // spacedim>::get_properties()’
      //
      // See https://github.com/dealii/dealii/issues/17148
      Assert(state() == IteratorState::valid, ExcInternalError());

      return property_pool->get_properties(get_handle());
    }


    /**
     * Get read-access to properties of this particle.
     *
     * @return An ArrayView of the properties of this particle.
     */
    ArrayView<const double>
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
     */
    const typename Triangulation<dim, spacedim>::cell_iterator &
    get_surrounding_cell() const;

    /**
     * Write the data of this object to a stream for the purpose of
     * serialization using the [BOOST serialization
     * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
     */
    template <class Archive>
    void
    save(Archive &ar, const unsigned int version) const;

    /**
     * Read the data of this object from a stream for the purpose of
     * serialization using the [BOOST serialization
     * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
     */
    template <class Archive>
    void
    load(Archive &ar, const unsigned int version);

#ifdef DOXYGEN
    /**
     * Write and read the data of this object from a stream for the purpose
     * of serialization using the [BOOST serialization
     * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
     */
    template <class Archive>
    void
    serialize(Archive &archive, const unsigned int version);
#else
    // This macro defines the serialize() method that is compatible with
    // the templated save() and load() method that have been implemented.
    BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif

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
      const typename particle_container::iterator particles_in_cell,
      const PropertyPool<dim, spacedim>          &property_pool,
      const unsigned int                          particle_index_within_cell);

    /**
     * Returns a reference to the current Particle. Because the
     * internal structure may change this is not intended for public use and
     * only a convenience function for internal purposes.
     */
    const typename PropertyPool<dim, spacedim>::Handle &
    get_handle() const;

    /**
     * Non-@p const version of the function with the same name above.
     */
    typename PropertyPool<dim, spacedim>::Handle &
    get_handle();

    /**
     * An iterator to the particles in the current cell within the
     * particle_container object.
     */
    typename particle_container::iterator particles_in_cell;

    /**
     * A pointer to the property pool that stores the actual particle data.
     */
    PropertyPool<dim, spacedim> *property_pool;

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
  inline void
  ParticleAccessor<dim, spacedim>::load(Archive &ar, const unsigned int)
  {
    unsigned int n_properties = 0;

    Point<spacedim>                       location;
    Point<dim>                            reference_location;
    types::particle_index                 id;
    ar &location &reference_location &id &n_properties;

    set_location(location);
    set_reference_location(reference_location);
    set_id(id);

    if (n_properties > 0)
      {
        ArrayView<double> properties(get_properties());
        Assert(
          properties.size() == n_properties,
          ExcMessage(
            "This particle was serialized with " +
            std::to_string(n_properties) +
            " properties, but the new property handler provides space for " +
            std::to_string(properties.size()) +
            " properties. Deserializing a particle only works for matching property sizes."));

        ar &boost::serialization::make_array(properties.data(), n_properties);
      }
  }



  template <int dim, int spacedim>
  template <class Archive>
  inline void
  ParticleAccessor<dim, spacedim>::save(Archive &ar, const unsigned int) const
  {
    unsigned int n_properties = 0;
    if ((property_pool != nullptr) &&
        (get_handle() != PropertyPool<dim, spacedim>::invalid_handle))
      n_properties = get_properties().size();

    Point<spacedim>       location           = get_location();
    Point<dim>            reference_location = get_reference_location();
    types::particle_index id                 = get_id();

    ar &location &reference_location &id &n_properties;

    if (n_properties > 0)
      ar &boost::serialization::make_array(get_properties().data(),
                                           n_properties);
  }


  // ------------------------- inline functions ------------------------------

  template <int dim, int spacedim>
  inline ParticleAccessor<dim, spacedim>::ParticleAccessor()
    : particles_in_cell(typename particle_container::iterator())
    , property_pool(nullptr)
    , particle_index_within_cell(numbers::invalid_unsigned_int)
  {}



  template <int dim, int spacedim>
  inline ParticleAccessor<dim, spacedim>::ParticleAccessor(
    const typename particle_container::iterator particles_in_cell,
    const PropertyPool<dim, spacedim>          &property_pool,
    const unsigned int                          particle_index_within_cell)
    : particles_in_cell(particles_in_cell)
    , property_pool(const_cast<PropertyPool<dim, spacedim> *>(&property_pool))
    , particle_index_within_cell(particle_index_within_cell)
  {}



  template <int dim, int spacedim>
  inline const void *
  ParticleAccessor<dim, spacedim>::read_particle_data_from_memory(
    const void *data)
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

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
        const ArrayView<double> particle_properties =
          property_pool->get_properties(get_handle());
        const unsigned int size = particle_properties.size();
        for (unsigned int i = 0; i < size; ++i)
          particle_properties[i] = *pdata++;
      }

    return static_cast<const void *>(pdata);
  }



  template <int dim, int spacedim>
  inline void *
  ParticleAccessor<dim, spacedim>::write_particle_data_to_memory(
    void *data) const
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    types::particle_index *id_data = static_cast<types::particle_index *>(data);
    *id_data                       = get_id();
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
          property_pool->get_properties(get_handle());
        for (unsigned int i = 0; i < particle_properties.size(); ++i, ++pdata)
          *pdata = particle_properties[i];
      }

    return static_cast<void *>(pdata);
  }



  template <int dim, int spacedim>
  inline void
  ParticleAccessor<dim, spacedim>::set_location(const Point<spacedim> &new_loc)
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    property_pool->set_location(get_handle(), new_loc);
  }



  template <int dim, int spacedim>
  inline const Point<spacedim> &
  ParticleAccessor<dim, spacedim>::get_location() const
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    return property_pool->get_location(get_handle());
  }



  template <int dim, int spacedim>
  inline Point<spacedim> &
  ParticleAccessor<dim, spacedim>::get_location()
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    return property_pool->get_location(get_handle());
  }



  template <int dim, int spacedim>
  inline void
  ParticleAccessor<dim, spacedim>::set_reference_location(
    const Point<dim> &new_loc)
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    property_pool->set_reference_location(get_handle(), new_loc);
  }



  template <int dim, int spacedim>
  inline const Point<dim> &
  ParticleAccessor<dim, spacedim>::get_reference_location() const
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    return property_pool->get_reference_location(get_handle());
  }



  template <int dim, int spacedim>
  inline types::particle_index
  ParticleAccessor<dim, spacedim>::get_id() const
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    return property_pool->get_id(get_handle());
  }



  template <int dim, int spacedim>
  inline types::particle_index
  ParticleAccessor<dim, spacedim>::get_local_index() const
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    return get_handle();
  }



  template <int dim, int spacedim>
  inline void
  ParticleAccessor<dim, spacedim>::set_id(const types::particle_index &new_id)
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    property_pool->set_id(get_handle(), new_id);
  }



  template <int dim, int spacedim>
  inline bool
  ParticleAccessor<dim, spacedim>::has_properties() const
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    // Particles always have a property pool associated with them,
    // but we can access properties only if there is a valid handle.
    // The only way a particle can have no valid handle if it has
    // been moved-from -- but that leaves an object in an invalid
    // state, and so we can just assert that that can't be the case.
    Assert((get_handle() != PropertyPool<dim, spacedim>::invalid_handle),
           ExcInternalError());
    return (property_pool->n_properties_per_slot() > 0);
  }



  template <int dim, int spacedim>
  inline void
  ParticleAccessor<dim, spacedim>::set_properties(
    const std::vector<double> &new_properties)
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    set_properties(
      ArrayView<const double>(new_properties.data(), new_properties.size()));
  }



  template <int dim, int spacedim>
  inline void
  ParticleAccessor<dim, spacedim>::set_properties(
    const ArrayView<const double> &new_properties)
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    const ArrayView<double> property_values =
      property_pool->get_properties(get_handle());

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
  inline void
  ParticleAccessor<dim, spacedim>::set_properties(
    const Tensor<1, dim> &new_properties)
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    // A Tensor object is not an array, so we cannot just create an
    // ArrayView object for it. Rather, copy the data into a true
    // array and make the ArrayView from that.
    double array[dim];
    for (unsigned int d = 0; d < dim; ++d)
      array[d] = new_properties[d];

    set_properties(make_array_view(array));
  }



  template <int dim, int spacedim>
  inline ArrayView<const double>
  ParticleAccessor<dim, spacedim>::get_properties() const
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    return property_pool->get_properties(get_handle());
  }



  template <int dim, int spacedim>
  inline const typename Triangulation<dim, spacedim>::cell_iterator &
  ParticleAccessor<dim, spacedim>::get_surrounding_cell() const
  {
    Assert(state() == IteratorState::valid, ExcInternalError());
    Assert(particles_in_cell->cell.state() == IteratorState::valid,
           ExcInternalError());

    return particles_in_cell->cell;
  }



  // template <int dim, int spacedim>
  // inline ArrayView<double>
  //   ParticleAccessor<dim, spacedim>::get_properties()



  template <int dim, int spacedim>
  inline std::size_t
  ParticleAccessor<dim, spacedim>::serialized_size_in_bytes() const
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    std::size_t size = sizeof(get_id()) +
                       sizeof(double) * spacedim + // get_location()
                       sizeof(double) * dim;       // get_reference_location()

    if (has_properties())
      {
        size += sizeof(double) * get_properties().size();
      }
    return size;
  }



  template <int dim, int spacedim>
  inline void
  ParticleAccessor<dim, spacedim>::next()
  {
    Assert(state() == IteratorState::valid, ExcInternalError());

    ++particle_index_within_cell;

    if (particle_index_within_cell >= particles_in_cell->particles.size())
      {
        particle_index_within_cell = 0;
        ++particles_in_cell;
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
        --particles_in_cell;
        particle_index_within_cell = particles_in_cell->particles.empty() ?
                                       0 :
                                       particles_in_cell->particles.size() - 1;
      }
  }



  template <int dim, int spacedim>
  inline bool
  ParticleAccessor<dim, spacedim>::operator!=(
    const ParticleAccessor<dim, spacedim> &other) const
  {
    return !(*this == other);
  }



  template <int dim, int spacedim>
  inline bool
  ParticleAccessor<dim, spacedim>::operator==(
    const ParticleAccessor<dim, spacedim> &other) const
  {
    return (property_pool == other.property_pool) &&
           (particles_in_cell == other.particles_in_cell) &&
           (particle_index_within_cell == other.particle_index_within_cell);
  }



  template <int dim, int spacedim>
  inline IteratorState::IteratorStates
  ParticleAccessor<dim, spacedim>::state() const
  {
    if (property_pool != nullptr &&
        particles_in_cell->cell.state() == IteratorState::valid &&
        particle_index_within_cell < particles_in_cell->particles.size())
      return IteratorState::valid;
    else if (property_pool != nullptr &&
             particles_in_cell->cell.state() == IteratorState::past_the_end &&
             particle_index_within_cell == 0)
      return IteratorState::past_the_end;
    else
      return IteratorState::invalid;

    return IteratorState::invalid;
  }



  template <int dim, int spacedim>
  inline typename PropertyPool<dim, spacedim>::Handle &
  ParticleAccessor<dim, spacedim>::get_handle()
  {
    return particles_in_cell->particles[particle_index_within_cell];
  }



  template <int dim, int spacedim>
  inline const typename PropertyPool<dim, spacedim>::Handle &
  ParticleAccessor<dim, spacedim>::get_handle() const
  {
    return particles_in_cell->particles[particle_index_within_cell];
  }

} // namespace Particles

DEAL_II_NAMESPACE_CLOSE

namespace boost
{
  namespace geometry
  {
    namespace index
    {
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
