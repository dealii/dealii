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


DEAL_II_NAMESPACE_OPEN

namespace Particles
{
  /**
   * This class manages a memory space in which all particles associated with
   * a ParticleHandler store their properties. The rationale for this class is
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
     * The currently allocated properties (whether assigned to
     * a particle or available for assignment).
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
  inline ArrayView<double>
  PropertyPool<dim, spacedim>::get_properties(const Handle handle)
  {
    const std::vector<double>::size_type data_index =
      (handle != invalid_handle) ? handle * n_properties : 0;

    // Ideally we would need to assert that 'handle' has not been deallocated
    // by searching through 'currently_available_handles'. However, this
    // is expensive and this function is performance critical, so instead
    // just check against the properties range, and rely on the fact
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
