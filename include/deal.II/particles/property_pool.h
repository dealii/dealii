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

#ifndef dealii_particles_property_pool_h
#define dealii_particles_property_pool_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/point.h>


DEAL_II_NAMESPACE_OPEN

namespace types
{
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
  using particle_index = uint64_t;

#  ifdef DEAL_II_WITH_MPI
  /**
   * An identifier that denotes the MPI type associated with
   * types::global_dof_index.
   */
#    define DEAL_II_PARTICLE_INDEX_MPI_TYPE MPI_UINT64_T
#  endif

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
  using particle_index = unsigned int;

#  ifdef DEAL_II_WITH_MPI
  /**
   * An identifier that denotes the MPI type associated with
   * types::global_dof_index.
   */
#    define DEAL_II_PARTICLE_INDEX_MPI_TYPE MPI_UNSIGNED
#  endif
#endif
} // namespace types

namespace Particles
{
  /**
   * This class manages a memory space in which all particles associated with
   * a ParticleHandler store their properties. It also stores the locations
   * and reference locations of particles.
   *
   * The rationale for this class is
   * that because typically every particle stores the same number of
   * properties, and because algorithms generally traverse over all particles
   * doing the same operation on all particles' properties, it is more efficient
   * to let the memory used for properties be handled by a central manager.
   * Particles then do not store a pointer to a memory area in which they store
   * their properties, but instead a "handle" that the PropertyPool class then
   * translates into a pointer to concrete memory.
   *
   * All this said, the current implementation only provides this kind of
   * interface, but still uses simple new/delete allocation for every
   * set of properties requested by a particle. Additionally, the current
   * implementation assumes the same number of properties per particle, but of
   * course the PropertyType could contain a pointer to dynamically allocated
   * memory with varying sizes per particle (this memory would not be managed by
   * this class).
   */
  template <int dim, int spacedim = dim>
  class PropertyPool
  {
  public:
    /**
     * Typedef for the handle that is returned to the particles, and that
     * uniquely identifies the slot of memory that is reserved for this
     * particle.
     */
    using Handle = unsigned int;

    /**
     * Define a default (invalid) value for handles.
     */
    static const Handle invalid_handle;

    /**
     * Constructor. Stores the number of properties per reserved slot.
     */
    PropertyPool(const unsigned int n_properties_per_slot);

    /**
     * Destructor. This function ensures that all memory that had
     * previously been allocated using allocate_properties_array()
     * has also been returned via deallocate_properties_array().
     */
    ~PropertyPool();

    /**
     * Clear the dynamic memory allocated by this class. This function
     * ensures that all memory that had previously been allocated using
     * allocate_properties_array() has also been returned via
     * deallocate_properties_array().
     */
    void
    clear();

    /**
     * Return a new handle that allows a particle to store information such as
     * properties and locations. This also allocated memory in this PropertyPool
     * variable.
     */
    Handle
    register_particle();

    /**
     * Return a handle obtained by register_particle() and mark the memory
     * allocated for storing the particle's data as free for re-use.
     */
    void
    deregister_particle(Handle &handle);

    /**
     * Return the location of a particle identified by the given `handle`.
     */
    const Point<spacedim> &
    get_location(const Handle handle) const;

    /**
     * Set the location of a particle identified by the given `handle`.
     */
    void
    set_location(const Handle handle, const Point<spacedim> &new_location);

    /**
     * Return the reference_location of a particle identified by the given
     * `handle`.
     */
    const Point<dim> &
    get_reference_location(const Handle handle) const;

    /**
     * Set the reference location of a particle identified by the given
     * `handle`.
     */
    void
    set_reference_location(const Handle      handle,
                           const Point<dim> &new_reference_location);

    /**
     * Return the ID number of this particle identified by the given
     * `handle`.
     */
    types::particle_index
    get_id(const Handle handle) const;

    /**
     * Set the ID number of this particle identified by the given
     * `handle`.
     */
    void
    set_id(const Handle handle, const types::particle_index &new_id);

    /**
     * Return an ArrayView to the properties that correspond to the given
     * handle @p handle.
     */
    ArrayView<double>
    get_properties(const Handle handle);

    /**
     * Reserve the dynamic memory needed for storing the properties of
     * @p size particles.
     */
    void
    reserve(const std::size_t size);

    /**
     * Return how many properties are stored per slot in the pool.
     */
    unsigned int
    n_properties_per_slot() const;

  private:
    /**
     * The number of properties that are reserved per particle.
     */
    const unsigned int n_properties;

    /**
     * A vector that stores the locations of particles. It is indexed in the
     * same way as the `reference_locations` and `properties` arrays, i.e., via
     * handles.
     */
    std::vector<Point<spacedim>> locations;

    /**
     * A vector that stores the reference locations of particles. It is indexed
     * in the same way as the `locations` and `properties` arrays, i.e., via
     * handles.
     */
    std::vector<Point<dim>> reference_locations;

    /**
     * A vector that stores the unique identifiers of particles. It is indexed
     * in the same way as the `locations` and `properties` arrays, i.e., via
     * handles.
     */
    std::vector<types::particle_index> ids;

    /**
     * The currently allocated properties (whether assigned to
     * a particle or available for assignment). It is indexed the same way as
     * the `locations` and `reference_locations` arrays via handles.
     */
    std::vector<double> properties;

    /**
     * A collection of handles that have been created by
     * allocate_properties_array() and have been destroyed by
     * deallocate_properties_array(). Since the memory is still
     * allocated these handles can be reused for new particles
     * to avoid memory allocation.
     */
    std::vector<Handle> currently_available_handles;
  };



  /* ---------------------- inline and template functions ------------------ */

  template <int dim, int spacedim>
  inline const Point<spacedim> &
  PropertyPool<dim, spacedim>::get_location(const Handle handle) const
  {
    const std::vector<double>::size_type data_index =
      (handle != invalid_handle) ? handle : 0;

    // Ideally we would need to assert that 'handle' has not been deallocated
    // by searching through 'currently_available_handles'. However, this
    // is expensive and this function is performance critical, so instead
    // just check against the array range, and rely on the fact
    // that handles are invalidated when handed over to
    // deallocate_properties_array().
    Assert(data_index <= locations.size() - 1,
           ExcMessage("Invalid location handle. This can happen if the "
                      "handle was duplicated and then one copy was deallocated "
                      "before trying to access the properties."));

    return locations[data_index];
  }



  template <int dim, int spacedim>
  inline void
  PropertyPool<dim, spacedim>::set_location(const Handle           handle,
                                            const Point<spacedim> &new_location)
  {
    const std::vector<double>::size_type data_index =
      (handle != invalid_handle) ? handle : 0;

    // Ideally we would need to assert that 'handle' has not been deallocated
    // by searching through 'currently_available_handles'. However, this
    // is expensive and this function is performance critical, so instead
    // just check against the array range, and rely on the fact
    // that handles are invalidated when handed over to
    // deallocate_properties_array().
    Assert(data_index <= locations.size() - 1,
           ExcMessage("Invalid location handle. This can happen if the "
                      "handle was duplicated and then one copy was deallocated "
                      "before trying to access the properties."));

    locations[data_index] = new_location;
  }



  template <int dim, int spacedim>
  inline const Point<dim> &
  PropertyPool<dim, spacedim>::get_reference_location(const Handle handle) const
  {
    const std::vector<double>::size_type data_index =
      (handle != invalid_handle) ? handle : 0;

    // Ideally we would need to assert that 'handle' has not been deallocated
    // by searching through 'currently_available_handles'. However, this
    // is expensive and this function is performance critical, so instead
    // just check against the array range, and rely on the fact
    // that handles are invalidated when handed over to
    // deallocate_properties_array().
    Assert(data_index <= reference_locations.size() - 1,
           ExcMessage("Invalid location handle. This can happen if the "
                      "handle was duplicated and then one copy was deallocated "
                      "before trying to access the properties."));

    return reference_locations[data_index];
  }



  template <int dim, int spacedim>
  inline void
  PropertyPool<dim, spacedim>::set_reference_location(
    const Handle      handle,
    const Point<dim> &new_reference_location)
  {
    const std::vector<double>::size_type data_index =
      (handle != invalid_handle) ? handle : 0;

    // Ideally we would need to assert that 'handle' has not been deallocated
    // by searching through 'currently_available_handles'. However, this
    // is expensive and this function is performance critical, so instead
    // just check against the array range, and rely on the fact
    // that handles are invalidated when handed over to
    // deallocate_properties_array().
    Assert(data_index <= reference_locations.size() - 1,
           ExcMessage("Invalid location handle. This can happen if the "
                      "handle was duplicated and then one copy was deallocated "
                      "before trying to access the properties."));

    reference_locations[data_index] = new_reference_location;
  }



  template <int dim, int spacedim>
  inline types::particle_index
  PropertyPool<dim, spacedim>::get_id(const Handle handle) const
  {
    const std::vector<double>::size_type data_index =
      (handle != invalid_handle) ? handle : 0;

    // Ideally we would need to assert that 'handle' has not been deallocated
    // by searching through 'currently_available_handles'. However, this
    // is expensive and this function is performance critical, so instead
    // just check against the array range, and rely on the fact
    // that handles are invalidated when handed over to
    // deallocate_properties_array().
    Assert(data_index <= ids.size() - 1,
           ExcMessage("Invalid location handle. This can happen if the "
                      "handle was duplicated and then one copy was deallocated "
                      "before trying to access the properties."));

    return ids[data_index];
  }



  template <int dim, int spacedim>
  inline void
  PropertyPool<dim, spacedim>::set_id(const Handle                 handle,
                                      const types::particle_index &new_id)
  {
    const std::vector<double>::size_type data_index =
      (handle != invalid_handle) ? handle : 0;

    // Ideally we would need to assert that 'handle' has not been deallocated
    // by searching through 'currently_available_handles'. However, this
    // is expensive and this function is performance critical, so instead
    // just check against the array range, and rely on the fact
    // that handles are invalidated when handed over to
    // deallocate_properties_array().
    Assert(data_index <= ids.size() - 1,
           ExcMessage("Invalid location handle. This can happen if the "
                      "handle was duplicated and then one copy was deallocated "
                      "before trying to access the properties."));

    ids[data_index] = new_id;
  }



  template <int dim, int spacedim>
  inline ArrayView<double>
  PropertyPool<dim, spacedim>::get_properties(const Handle handle)
  {
    const std::vector<double>::size_type data_index =
      (handle != invalid_handle) ? handle * n_properties : 0;

    // Ideally we would need to assert that 'handle' has not been deallocated
    // by searching through 'currently_available_handles'. However, this
    // is expensive and this function is performance critical, so instead
    // just check against the array range, and rely on the fact
    // that handles are invalidated when handed over to
    // deallocate_properties_array().
    Assert(data_index <= properties.size() - n_properties,
           ExcMessage("Invalid property handle. This can happen if the "
                      "handle was duplicated and then one copy was deallocated "
                      "before trying to access the properties."));

    return ArrayView<double>(properties.data() + data_index, n_properties);
  }


} // namespace Particles

DEAL_II_NAMESPACE_CLOSE

#endif
