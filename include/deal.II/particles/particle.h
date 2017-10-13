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

#ifndef dealii_particles_particle_h
#define dealii_particles_particle_h

#include <deal.II/particles/property_pool.h>

#include <deal.II/base/point.h>
#include <deal.II/base/types.h>
#include <deal.II/base/array_view.h>

DEAL_II_NAMESPACE_OPEN

/**
 * A namespace for all type definitions related to particles.
 */
namespace types
{
  /**
   * Typedef of cell level/index pair. TODO: replace this by the
   * active_cell_index.
   */
  typedef std::pair<int, int> LevelInd;

  /* Type definitions */

#ifdef DEAL_II_WITH_64BIT_INDICES
  /**
   * The type used for indices of particles. While in
   * sequential computations the 4 billion indices of 32-bit unsigned integers
   * is plenty, parallel computations using hundreds of processes can overflow
   * this number and we need a bigger index space. We here utilize the same
   * build variable that controls the dof indices because the number
   * of degrees of freedom and the number of particles are typically on the same
   * order of magnitude.
   *
   * The data type always indicates an unsigned integer type.
   */
  typedef unsigned long long int particle_index;

  /**
   * An identifier that denotes the MPI type associated with
   * types::global_dof_index.
   */
#  define PARTICLE_INDEX_MPI_TYPE MPI_UNSIGNED_LONG_LONG
#else
  /**
   * The type used for indices of particles. While in
   * sequential computations the 4 billion indices of 32-bit unsigned integers
   * is plenty, parallel computations using hundreds of processes can overflow
   * this number and we need a bigger index space. We here utilize the same
   * build variable that controls the dof indices because the number
   * of degrees of freedom and the number of particles are typically on the same
   * order of magnitude.
   *
   * The data type always indicates an unsigned integer type.
   */
  typedef unsigned int particle_index;

  /**
   * An identifier that denotes the MPI type associated with
   * types::global_dof_index.
   */
#  define PARTICLE_INDEX_MPI_TYPE MPI_UNSIGNED
#endif
}

/**
 * A namespace that contains all classes that are related to the particle
 * implementation, in particular the fundamental Particle class.
 */
namespace Particles
{
  /**
   * Base class of particles - represents a particle with position,
   * an ID number and a variable number of properties. This class
   * can be extended to include data related to a particle by the property
   * manager.
   *
   * @ingroup Particle
   * @author Rene Gassmoeller, 2017
   *
   */
  template <int dim, int spacedim=dim>
  class Particle
  {
  public:
    /**
     * Empty constructor for Particle, creates a particle at the
     * origin.
     */
    Particle ();

    /**
     * Constructor for Particle, creates a particle with the specified
     * ID at the specified location. Note that there is no
     * check for duplicate particle IDs so the user must
     * make sure the IDs are unique over all processes.
     *
     * @param[in] location Initial location of particle.
     * @param[in] reference_location Initial location of the particle
     * in the coordinate system of the reference cell.
     * @param[in] id Globally unique ID number of particle.
     */
    Particle (const Point<spacedim> &location,
              const Point<dim> &reference_location,
              const types::particle_index id);

    /**
     * Copy-Constructor for Particle, creates a particle with exactly the
     * state of the input argument. Note that since each particle has a
     * handle for a certain piece of the property memory, and is responsible
     * for registering and freeing this memory in the property pool this
     * constructor registers a new chunk, and copies the properties.
     */
    Particle (const Particle<dim,spacedim> &particle);

    /**
     * Constructor for Particle, creates a particle from a data vector.
     * This constructor is usually called after serializing a particle by
     * calling the write_data function.
     *
     * @param[in,out] begin_data A pointer to a memory location from which
     * to read the information that completely describes a particle. This
     * class then de-serializes its data from this memory location and
     * advance the pointer accordingly.
     *
     * @param[in,out] property_pool An optional pointer to a property pool
     * that is used to manage the property data used by this particle. Note that
     * if a non-null pointer is handed over this constructor assumes @p begin_data
     * contains serialized data of the same length and type that is allocated
     * by @p property_pool.
     */
    Particle (const void *&begin_data,
              PropertyPool &property_pool);

    /**
     * Move constructor for Particle, creates a particle from an existing
     * one by stealing its state.
     */
    Particle (Particle<dim,spacedim> &&particle);

    /**
     * Copy assignment operator.
     */
    Particle<dim,spacedim> &operator=(const Particle<dim,spacedim> &particle);

    /**
     * Move assignment operator.
     */
    Particle<dim,spacedim> &operator=(Particle<dim,spacedim> &&particle);

    /**
     * Destructor. Releases the property handle if it is valid, and
     * therefore frees that memory space for other particles. (Note:
     * the memory is managed by the property pool, and the pool is responsible
     * for what happens to the memory.
     */
    ~Particle ();

    /**
     * Write particle data into a data array. The array is expected
     * to be large enough to take the data, and the void pointer should
     * point to the first element in which the data should be written. This
     * function is meant for serializing all particle properties and
     * afterwards de-serializing the properties by calling the appropriate
     * constructor Particle(void *&data, PropertyPool *property_pool = NULL);
     *
     * @param [in,out] data The memory location to write particle data
     * into. This pointer points to the begin of the memory, in which the
     * data will be written and it will be advanced by the serialized size
     * of this particle.
     */
    void
    write_data(void *&data) const;

    /**
      * Set the location of this particle. Note that this does not check
      * whether this is a valid location in the simulation domain.
      *
      * @param [in] new_location The new location for this particle.
      */
    void
    set_location (const Point<spacedim> &new_location);

    /**
     * Get the location of this particle.
     *
     * @return The location of this particle.
     */
    const Point<spacedim> &
    get_location () const;

    /**
     * Set the reference location of this particle.
     *
     * @param [in] new_reference_location The new reference location for
     * this particle.
     */
    void
    set_reference_location (const Point<dim> &new_reference_location);

    /**
     * Return the reference location of this particle in its current cell.
     */
    const Point<dim> &
    get_reference_location () const;

    /**
     * Return the ID number of this particle.
     */
    types::particle_index
    get_id () const;

    /**
     * Tell the particle where to store its properties (even if it does not
     * own properties). Usually this is only done once per particle, but
     * since the particle does not know about the properties,
     * we want to do it not at construction time. Another use for this
     * function is after particle transfer to a new process.
     */
    void
    set_property_pool(PropertyPool &property_pool);

    /**
      * Returns whether this particle has a valid property pool and a valid
      * handle to properties.
      */
    bool
    has_properties () const;

    /**
     * Set the properties of this particle.
     *
     * @param [in] new_properties An ArrayView containing the
     * new properties for this particle.
     */
    void
    set_properties (const ArrayView<const double> &new_properties);

    /**
     * Get write-access to properties of this particle.
     *
     * @return An ArrayView of the properties of this particle.
     */
    const ArrayView<double>
    get_properties ();

    /**
     * Get read-access to properties of this particle.
     *
     * @return An ArrayView of the properties of this particle.
     */
    const ArrayView<const double>
    get_properties () const;

    /**
     * Returns the size in bytes this particle occupies if all of its data is
     * serialized (i.e. the number of bytes that is written by the write_data
     * function of this class).
     */
    std::size_t
    serialized_size_in_bytes() const;

    /**
     * Write the data of this object to a stream for the purpose of
     * serialization.
     */
    template <class Archive>
    void save (Archive &ar, const unsigned int version) const;

    /**
     * Read the data of this object from a stream for the purpose of
     * serialization.
     */
    template <class Archive>
    void load (Archive &ar, const unsigned int version);

    BOOST_SERIALIZATION_SPLIT_MEMBER()

  private:
    /**
     * Current particle location.
     */
    Point<spacedim>             location;

    /**
     * Current particle location in the reference cell.
     */
    Point<dim>             reference_location;

    /**
     * Globally unique ID of particle.
     */
    types::particle_index  id;

    /**
     * A pointer to the property pool. Necessary to translate from the
     * handle to the actual memory locations.
     */
    PropertyPool *property_pool;

    /**
     * A handle to all particle properties
     */
    PropertyPool::Handle properties;
  };

  /* -------------------------- inline and template functions ---------------------- */

  template <int dim, int spacedim>
  template <class Archive>
  void Particle<dim,spacedim>::load (Archive &ar, const unsigned int)
  {
    unsigned int n_properties = 0;

    ar &location
    & reference_location
    & id
    & n_properties;

    if (n_properties > 0)
      {
        properties = new double[n_properties];
        ar &boost::serialization::make_array(properties, n_properties);
      }
  }

  template <int dim, int spacedim>
  template <class Archive>
  void Particle<dim,spacedim>::save (Archive &ar, const unsigned int) const
  {
    unsigned int n_properties = 0;
    if ((property_pool != NULL) && (properties != PropertyPool::invalid_handle))
      n_properties = get_properties().size();

    ar &location
    & reference_location
    & id
    & n_properties;

    if (n_properties > 0)
      ar &boost::serialization::make_array(properties, n_properties);

  }
}

DEAL_II_NAMESPACE_CLOSE

#endif

