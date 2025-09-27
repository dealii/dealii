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

#ifndef dealii_particles_particle_h
#define dealii_particles_particle_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/point.h>
#include <deal.II/base/types.h>

#include <deal.II/particles/property_pool.h>

#include <boost/geometry/index/indexable.hpp>
#include <boost/serialization/array.hpp>

#include <cstdint>


DEAL_II_NAMESPACE_OPEN

/**
 * A namespace that contains all classes that are related to the particle
 * implementation, in particular the fundamental Particle class.
 */
namespace Particles
{
  namespace internal
  {
    /**
     * Internal alias of cell level/index pair.
     */
    using LevelInd = std::pair<int, int>;
  } // namespace internal

  /**
   * A class that represents a particle in a domain that is meshed by
   * a triangulation of some kind. The data this class stores is the
   * position of the particle in the overall space, the position of
   * the particle in the reference coordinate system of the cell it is
   * currently in, an ID number that is unique among all particles,
   * and a variable number of "properties".
   *
   * The "properties" attached to each object of this class are
   * stored by a PropertyPool object. These properties are
   * stored as an array of `double` variables that can be accessed
   * via an ArrayView object. For example, if one wanted to equip
   * each particle with a "temperature" and "chemical composition"
   * property that is advected along with the particle (and may change
   * from time step to time step based on some differential equation,
   * for example), then one would allocate two properties per particle
   * in the PropertyPool object.
   *
   * In practice, however, one often wants to associate properties
   * with particles that are not just independent numbers as in the
   * situation above. An example would be if one wanted to track the
   * stress or strain that a particle is subjected to -- a tensor-valued
   * quantity. In these cases, one would <i>interpret</i> these scalar
   * properties as the <i>components of the stress or strain</i>. In
   * other words, one would first tell the PropertyPool to allocate
   * as many properties per particle as there are components in the
   * tensor one wants to track, and then write small conversion functions that
   * take the ArrayView of scalar properties returned by the
   * get_properties() function and convert it to a tensor of the
   * appropriate type. This can then be evaluated and evolved in each
   * time step. A second conversion function would convert back from a
   * tensor to an ArrayView object to store the updated data back in the
   * particle via the set_properties() function.
   *
   * There are of course cases where the properties one cares about are
   * not real (or, in computers, floating point) numbers but rather
   * categorical: For example, one may want to mark some particles
   * as "red", "blue", or "green". The property might then either be
   * represented as an integer, or as an element of an `enum`. In these
   * cases, one would need to come up with a way to <i>represent</i>
   * these sorts of categorical fields in terms of floating point
   * numbers. For example, one could map "red" to the floating point number
   * 1.0, "blue" to 2.0, and "green" to 3.0. The conversion functions
   * to translate between these two representations should then not be very
   * difficult to write either.
   */
  template <int dim, int spacedim = dim>
  class Particle
  {
  public:
    /**
     * Empty constructor for Particle, creates a particle at the
     * origin.
     */
    Particle();

    /**
     * Constructor for Particle. This function creates a particle with the
     * specified ID at the specified location. Note that there is no check for
     * duplicate particle IDs so the user must make sure the IDs are unique over
     * all processes. Data is stored in a global PropertyPool object
     * (corresponding to the global "heap") but can later be transferred to
     * another property pool by calling set_property_pool().
     *
     * @param[in] location Initial location of particle.
     * @param[in] reference_location Initial location of the particle
     * in the coordinate system of the reference cell.
     * @param[in] id Globally unique ID number of particle.
     */
    Particle(const Point<spacedim>      &location,
             const Point<dim>           &reference_location,
             const types::particle_index id);

    /**
     * Copy-constructor for Particle. This function creates a particle with
     * exactly the state of the input argument. The copied data is stored in a
     * global PropertyPool object (corresponding to the global "heap") but can
     * later be transferred to another property pool by calling
     * set_property_pool().
     */
    Particle(const Particle<dim, spacedim> &particle);

    /**
     * Constructor for Particle. This function creates a particle from a data
     * vector. Data is stored in a global PropertyPool object (corresponding to
     * the global "heap") but can later be transferred to another property pool
     * by calling set_property_pool(). This constructor is usually called after
     * serializing a particle by calling the write_data() function.
     *
     * @param[in,out] begin_data A pointer to a memory location from which
     * to read the information that completely describes a particle. This
     * class then de-serializes its data from this memory location and
     * advances the pointer beyond the data that has been read to initialize
     * the particle information.
     *
     * @param[in,out] property_pool An optional pointer to a property pool
     * that is used to manage the property data used by this particle. If this
     * argument is not provided, then a global property pool is used; on the
     * other hand,
     * if a non-null pointer is provided, this constructor assumes @p begin_data
     * contains serialized data of the same length and type that is allocated
     * by @p property_pool. If the data pointer provided here corresponds
     * to data for a particle that has properties, then this function will only
     * succeed if a property pool is provided as second argument that is able to
     * store the correct number of properties per particle.
     */
    Particle(const void                       *&begin_data,
             PropertyPool<dim, spacedim> *const property_pool = nullptr);

    /**
     * Move constructor for Particle, creates a particle from an existing
     * one by stealing its state.
     */
    Particle(Particle<dim, spacedim> &&particle) noexcept;

    /**
     * Copy assignment operator.
     */
    Particle<dim, spacedim> &
    operator=(const Particle<dim, spacedim> &particle);

    /**
     * Move assignment operator.
     */
    Particle<dim, spacedim> &
    operator=(Particle<dim, spacedim> &&particle) noexcept;

    /**
     * Destructor. Releases the property handle if it is valid, and
     * therefore frees that memory space for other particles. (Note:
     * the memory is managed by the property pool, and the pool is responsible
     * for what happens to the memory.
     */
    ~Particle();

    /**
     * Write particle data into a data array. The array is expected
     * to be large enough to take the data, and the void pointer should
     * point to the first entry of the array to which the data should be
     * written. This function is meant for serializing all particle properties
     * and later de-serializing the properties by calling the appropriate
     * constructor Particle(void *&data, PropertyPool *property_pool = nullptr);
     *
     * @param [in] data The memory location to write particle data
     *   into.
     *
     * @return A pointer to the next byte after the array to which data has
     *   been written.
     */
    void *
    write_particle_data_to_memory(void *data) const;


    /**
     * Update all of the data associated with a particle: id,
     * location, reference location and, if any, properties by using a
     * data array. The array is expected to be large enough to take the data,
     * and the void pointer should point to the first entry of the array to
     * which the data should be written. This function is meant for
     * de-serializing the particle data without requiring that a new Particle
     * class be built. This is used in the ParticleHandler to update the
     * ghost particles without de-allocating and re-allocating memory.
     *
     * @param[in] data A pointer to a memory location from which
     * to read the information that completely describes a particle. This
     * class then de-serializes its data from this memory location.
     *
     * @return A pointer to the next byte after the array from which data has
     *   been read.
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
     * Get the location of this particle.
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
     * Return the ID number of this particle. The ID of a particle is intended
     * to be a property that is globally unique even in parallel computations
     * and is transferred along with other properties of a particle if it
     * moves from a cell owned by the current processor to a cell owned by
     * a different processor, or if ownership of the cell it is on is
     * transferred to a different processor.
     */
    types::particle_index
    get_id() const;

    /**
     * Set the ID number of this particle. The ID of a particle is intended
     * to be a property that is globally unique even in parallel computations
     * and is transferred along with other properties of a particle if it
     * moves from a cell owned by the current processor to a cell owned by
     * a different processor, or if ownership of the cell it is on is
     * transferred to a different processor. As a consequence, when setting
     * the ID of a particle, care needs to be taken to ensure that particles
     * have globally unique IDs. (The ParticleHandler does not itself check
     * whether particle IDs so set are globally unique in a parallel setting
     * since this would be a very expensive operation.)
     *
     * @param[in] new_id The new ID number for this particle.
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
    set_id(const types::particle_index &new_id);

    /**
     * Tell the particle where to store its properties (even if it does not
     * own properties). Usually this is only done once per particle, but
     * since the particle does not know about the properties,
     * we want to do it not at construction time. Another use for this
     * function is after particle transfer to a new process.
     *
     * If a particle already stores properties in a property pool, then
     * their values are saved, the memory is released in the previous
     * property pool, and a copy of the particle's properties will be
     * allocated in the new property pool.
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
     * @param [in] new_properties An ArrayView containing the
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
     *   In both cases, one should
     *   treat the ghost elements or cells as `const` objects that shouldn't
     *   be modified even if the objects allow for calls that modify
     *   properties. Rather, properties should only be modified on processors
     *   that actually *own* the particle.
     */
    void
    set_properties(const ArrayView<const double> &new_properties);

    /**
     * Get write-access to properties of this particle. If the
     * particle has no properties yet, but has access to a
     * PropertyPool object it will allocate properties to
     * allow writing into them. If it has no properties and
     * has no access to a PropertyPool this function will
     * throw an exception.
     *
     * @return An ArrayView of the properties of this particle.
     */
    ArrayView<double>
    get_properties();

    /**
     * Get read-access to properties of this particle. If the particle
     * has no properties this function throws an exception.
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
     * Note that in order to store the properties correctly, the property pool
     * of this particle has to be known at the time of reading, i.e.
     * set_property_pool() has to have been called, before this function is
     * called.
     */
    template <class Archive>
    void
    load(Archive &ar, const unsigned int version);

    /**
     * Free the memory of the property pool
     */
    void
    free_properties();

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

  private:
    /**
     * A global property pool used when a particle is not associated with
     * a property pool that belongs to, for example, a ParticleHandler.
     */
    static PropertyPool<dim, spacedim> global_property_pool;

    /**
     * A pointer to the property pool. Necessary to translate from the
     * handle to the actual memory locations.
     */
    PropertyPool<dim, spacedim> *property_pool;

    /**
     * A handle to all particle properties
     */
    typename PropertyPool<dim, spacedim>::Handle property_pool_handle;
  };



  /* ---------------------- inline and template functions ------------------ */

  template <int dim, int spacedim>
  template <class Archive>
  inline void
  Particle<dim, spacedim>::load(Archive &ar, const unsigned int)
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
  Particle<dim, spacedim>::save(Archive &ar, const unsigned int) const
  {
    unsigned int n_properties = 0;
    if ((property_pool != nullptr) &&
        (property_pool_handle != PropertyPool<dim, spacedim>::invalid_handle))
      n_properties = get_properties().size();

    Point<spacedim>       location           = get_location();
    Point<dim>            reference_location = get_reference_location();
    types::particle_index id                 = get_id();

    ar &location &reference_location &id &n_properties;

    if (n_properties > 0)
      ar &boost::serialization::make_array(get_properties().data(),
                                           n_properties);
  }



  template <int dim, int spacedim>
  inline void
  Particle<dim, spacedim>::set_location(const Point<spacedim> &new_loc)
  {
    property_pool->set_location(property_pool_handle, new_loc);
  }



  template <int dim, int spacedim>
  inline const Point<spacedim> &
  Particle<dim, spacedim>::get_location() const
  {
    return property_pool->get_location(property_pool_handle);
  }



  template <int dim, int spacedim>
  inline Point<spacedim> &
  Particle<dim, spacedim>::get_location()
  {
    return property_pool->get_location(property_pool_handle);
  }



  template <int dim, int spacedim>
  inline void
  Particle<dim, spacedim>::set_reference_location(const Point<dim> &new_loc)
  {
    property_pool->set_reference_location(property_pool_handle, new_loc);
  }



  template <int dim, int spacedim>
  inline const Point<dim> &
  Particle<dim, spacedim>::get_reference_location() const
  {
    return property_pool->get_reference_location(property_pool_handle);
  }



  template <int dim, int spacedim>
  inline types::particle_index
  Particle<dim, spacedim>::get_id() const
  {
    return property_pool->get_id(property_pool_handle);
  }



  template <int dim, int spacedim>
  inline void
  Particle<dim, spacedim>::set_id(const types::particle_index &new_id)
  {
    property_pool->set_id(property_pool_handle, new_id);
  }



  template <int dim, int spacedim>
  inline void
  Particle<dim, spacedim>::set_property_pool(
    PropertyPool<dim, spacedim> &new_property_pool)
  {
    // First, we do want to save any properties that may
    // have previously been set, and copy them over to the memory allocated
    // on the new pool.
    //
    // It is possible that a particle currently has no properties -- for
    // example if it has been created without an associated property
    // pool (i.e., uses the default global pool which does not store any
    // properties) but that the new pool has properties. In that case,
    // there is simply nothing to transfer -- but the register_particle()
    // call here will make sure that the newly allocated properties are
    // zero-initialized.
    const typename PropertyPool<dim, spacedim>::Handle new_handle =
      new_property_pool.register_particle();

    const Point<spacedim>       location           = get_location();
    const Point<dim>            reference_location = get_reference_location();
    const types::particle_index id                 = get_id();

    if (/* old pool */ has_properties())
      {
        ArrayView<const double> old_properties = this->get_properties();
        ArrayView<double>       new_properties =
          new_property_pool.get_properties(new_handle);
        std::copy(old_properties.cbegin(),
                  old_properties.cend(),
                  new_properties.begin());
      }

    // Now release the old memory handle
    property_pool->deregister_particle(property_pool_handle);


    // Then set the pointer to the property pool we want to use. Also set the
    // handle to any properties.
    property_pool        = &new_property_pool;
    property_pool_handle = new_handle;

    // Now also store the saved locations
    set_location(location);
    set_reference_location(reference_location);
    set_id(id);
  }



  template <int dim, int spacedim>
  inline ArrayView<const double>
  Particle<dim, spacedim>::get_properties() const
  {
    if (has_properties() == false)
      return {};
    else
      return property_pool->get_properties(property_pool_handle);
  }



  template <int dim, int spacedim>
  inline bool
  Particle<dim, spacedim>::has_properties() const
  {
    // Particles always have a property pool associated with them,
    // but we can access properties only if there is a valid handle.
    // The only way a particle can have no valid handle if it has
    // been moved-from -- but that leaves an object in an invalid
    // state, and so we can just assert that that can't be the case.
    Assert((property_pool_handle !=
            PropertyPool<dim, spacedim>::invalid_handle),
           ExcInternalError());
    return (property_pool->n_properties_per_slot() > 0);
  }

} // namespace Particles

DEAL_II_NAMESPACE_CLOSE // Do not convert for module purposes


  namespace boost
{
  namespace geometry
  {
    namespace index
    {
      /**
       * Make sure we can construct an RTree of Particles::Particle objects.
       */
      template <int dim, int spacedim>
      struct indexable<dealii::Particles::Particle<dim, spacedim>>
      {
        /**
         * boost::rtree expects a const reference to an indexable object. For
         * a Particles::Particle object, this is its reference location.
         */
        using result_type = const dealii::Point<spacedim> &;

        result_type
        operator()(
          const dealii::Particles::Particle<dim, spacedim> &particle) const
        {
          return particle.get_location();
        }
      };

    } // namespace index
  }   // namespace geometry
} // namespace boost


DEAL_II_NAMESPACE_OPEN // Do not convert for module purposes
  DEAL_II_NAMESPACE_CLOSE

#endif
